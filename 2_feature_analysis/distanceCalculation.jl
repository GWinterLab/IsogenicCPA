#############################################
# !! WORKING VERSION !!
# Author(s): Amanda Ng R.H.
# Created on: 27 Mar 2023
# Last updated on: 18 Apr 2023
# Documentation status: Good
#
# Function for calculating the dissimilarity between treatments and a reference compound as Robust Hellinger Distance using the BioProfiling.jl package. (Parameters following Loan Vuillard's parameters in the example notebooks for BioProfiling.jl.)
# I wrote this script for calculating the distance in command line using an SBATCH script.
# Example use case: execute_distanceCalculation.sh for GW015_006 full-scale CPA
#
# Note: This script has been tested and is ready for use.
#############################################

#######################
# Precompile package(s)
#######################
using ArgParse
using CSV, StatsBase, Statistics, DataFrames, UMAP, RCall, FreqTables
using MultipleTesting, Random, MultivariateStats, Distributed
using BioProfiling
using Distances

#############
# Function(s)
#############
"""
    _print(message)

Apparently, Julia stalls the printing of messages till after a function has been carried out.
This stalling is quite annoying when I want to check on real time progress of an internal function.
To overcome this stalling, you need to flush the buffer that is added to the output messages of the
internal function.

This function allows you to print messages when implemented in functions.

# Argument(s)
 - `message::String`: the message to print out.

# Example(s)
```julia-repl
julia> _print("Hello world")
Hello world
```
"""
function _print(message::AbstractString)
    println(message)
    flush(stdout)
end

"""
    makeDirectory(directory)

Make a directory if it does not exist already.

# Argument(s)
 - `directory::String`: the full path of the directory to make or check if it exists already.
"""
function makeDirectory(directory)
    if isdir(directory) == false
        mkdir(directory)
        _print("The following directory has been made:\n$(directory)")
    else
        _print("The following directory already exists:\n$(directory)")
    end
end

"""
    exportDataFrame(
        output_path::String,
        dataframe::DataFrame;
        dataframe_description = ""
    )

Export a dataframe as a CSV file. If the user has provided a description for the dataframe being
exported, print the description as well.

# Argument(s)
 - `output_path::String`: The full path to export the dataframe as.
 - `dataframe::DataFrame`: Dataframe to export

# Optional argument(s)
 - `dataframe_description::String`: A description of the dataframe being exported.
"""
function exportDataFrame(
        output_path::String,
        dataframe::DataFrame;
        dataframe_description = ""
    )
    _print("Dataframe being exported: $(dataframe_description)")
    CSV.write(output_path, dataframe)
    _print("EXPORTED: Dataframe exported at\n$(output_path)")
end

"""
    retrieveFeatures(
        selectedFeatures_path::String,
        selectionStrategy::String,
        treatment::String
    )

Function for retrieving features using a specific selection strategy/treatment.

# Argument(s)
 - `selectedFeatures_path::String`: The full path to the CSV file with the features selected using different feature
    selection strategies
 - `selectionStrategy::String`: Name of the selection strategy of interest (e.g. "global", "treatment-ccentric")
 - `treatment::String`: Name of the treatment if selectionStrategy = "treatment-ccentric". Else, treatment is not
    used. Deafaults to "all".

# Return(s)
 - `features::Array`: List of features retrieved for the specific selection strategy/treatment.
"""
function retrieveFeatures(
        selectedFeatures_path::String,
        selectionStrategy::String,
        treatment = "all"
    )
    _print("Features selected by $(selectionStrategy) for $(treatment) treatment")
    selectedFeatures = CSV.read(selectedFeatures_path, DataFrame)
    if selectionStrategy == "treatment-ccentric"
        features = filter(:Treatment => x -> x ==(treatment), selectedFeatures)[1, :Features_Selected]
    else
        features = filter(:Selection_Strategy => x -> x ==(selectionStrategy), selectedFeatures)[1, :Features_Selected]
    end
    features = split(features, " ")
    return features
end

"""
    umapReduction(
        profile_path::String,
        cell::String,
        features::Array;
        n_components = 4,
        metric = CosineDist(),
        min_dist = 2
    )

Function for reducing selected features of a profile to a user-defined number of components by UMAP. The default optional arguments used here are from Loan Vuillard's implementation (https://pubmed.ncbi.nlm.nih.gov/34935929/).

# Argument(s)
 - `profile_path::String`: Full path to the morphology profile. The profile used should be AFTER robust Z
    standardization and baseline feature selection. Otherwise, batch/plate effects will be very obvious.
 - `cell::String`: Name of the cell line to use data from (e.g. "c662_rko_wt").
 - `features::Array`: List of features to use for UMAP reduction.

# Optional argument(s)
 - `n_components::Integer`: Number of UMAP components to reduce the profile to. Defaults to 4.
 - `metric::{SemiMetric, Symbol}`: The metric to use for calculating the distance in the input space. See
    https://github.com/dillondaudert/UMAP.jl/blob/master/src/umap_.jl for further information. Defaults to
    CosineDist().
 - `min_dist::Real`: The minimum spacing of points in the output embedding. Defaults to 2.
"""
function umapReduction(
        profile_path::String,
        cell::String,
        features::Array;
        n_components = 4,
        metric = CosineDist(),
        min_dist = 2
    )
    _print("Starting UMAP reduction...")
    
    # Define the column names to use for the UMAP-reduced profile
    columns = []
    n = 0
    while n < n_components
        n += 1
        columns = push!(columns, "UMAP$(n)")
    end
    columns = [Symbol(column) for column in columns]
    
    # Prep the profile for UMAP reduction
    profile = CSV.read(profile_path, DataFrame)
    if cell != "all"
        profile = filter(:Metadata_Cell => x -> x ==(cell), profile)
    end
    
    # Trim the profile to the features selected and the relevant metadata columns
    metadataToKeep = ["Metadata_Treatment", "Metadata_Cell", "Metadata_Well", "Metadata_Plate"]
    featuresToKeep = vcat(
        metadataToKeep,
        features
    )
    profile = profile[!, featuresToKeep]
    
    # Reduce the profile by UMAP to n_components
    matrix = convert(Matrix, profile[!, features])
    reducedProfile = umap(matrix', n_components; metric = metric, min_dist = min_dist)
    reducedProfile = convert(DataFrame, reducedProfile')
    rename!(reducedProfile, columns)
    _print("COMPLETED: UMAP reduction")
    
    # Add the metadata columns
    for metadata in metadataToKeep
        reducedProfile[!, metadata] = profile[!, metadata]
    end
    _print("First five rows of the UMAP reduced profile:")
    display(first(reducedProfile, 5))
    
    # Return the UMAP reduced profile
    return reducedProfile
end

"""
    calculateDistance(
        reducedProfile::DataFrame,
        referenceCompound::String
    )

Function for calculating the pairwise distance between all treatments and a reference compound. Default parameters for robust_morphological_perturbation_value used in Loan Vuillard's implementation (https://pubmed.ncbi.nlm.nih.gov/34935929/).

# Argument(s)
 - `reducedProfile::DataFrame`: UMAP-reduced profile.
 - `referenceCompound::String`: Name of the compound to use as a reference compound in distance calculation.
"""
function calculateDistance(
        reducedProfile::DataFrame,
        referenceCompound::String
    )
    _print("Starting distance calculation...")
    
    # Prep the reducedProfile for distance calculation
    expUMAP = Experiment(reducedProfile, description = "UMAP-reduced profile")
    filters = Array{BioProfiling.AbstractReduce, 1}()
    metadataToKeep = ["Metadata_Treatment", "Metadata_Cell", "Metadata_Well", "Metadata_Plate"]
    push!(filters, NameSelector(x -> !any(occursin.(metadataToKeep, String(x)))))
    filter!(expUMAP, filters)
    
    # Calculate the distance using Robust Hellinger Distance
    distance = robust_morphological_perturbation_value(
        expUMAP, :Metadata_Treatment, referenceCompound;
        nb_rep = 5000, dist = :RobHellinger
    )
    @assert !any(ismissing.(distance.RMPV))
    _print("COMPLETED: Distance calculation for $(referenceCompound)")
    return distance
end

####################################################################
# Execution of distance calculation by interfacing with command line
####################################################################
# Interfacing with command line
function parse_commandline()
    
    # Retrieve the arguments from command line
    s = ArgParseSettings()
    
    # Information to help the user use this julia function in command line
    @add_arg_table s begin
        "selectedFeatures_path"
            help = "Path to the CSV file containing the features selected with different feature selection strategies."
            arg_type = String
            required = true
        "selectionStrategy"
            help = "Name of the selection strategy."
            arg_type = String
            required = true
        "treatment"
            help = "Name of the treatment or 'all'."
            arg_type = String
            required = true
        "profile_path"
            help = "Path to the CSV file containing the feature measurements after robust Z standardization and baseline feature selection."
            arg_type = String
            required = true
        "cell"
            help = "Name of cell line to use (e.g. c662_rko_wt)."
            arg_type = String
            required = true
        "referenceCompound"
            help = "Name of the compound to use as a reference for distance calculation."
            arg_type = String
            required = true
        "output_dir"
            help = "Path to the directory for exporting the outputs to."
            arg_type = String
            required = true
    end
    
    return parse_args(s)
end

# Function for executing the functions for UMAP reduction and distance calculation
function main()
    
    # Retrieve the arguments from command line
    parsed_args = parse_commandline()
    selectedFeatures_path = parsed_args["selectedFeatures_path"]
    selectionStrategy = parsed_args["selectionStrategy"]
    treatment = parsed_args["treatment"]
    profile_path = parsed_args["profile_path"]
    cell = parsed_args["cell"]
    referenceCompound = parsed_args["referenceCompound"]
    output_dir = parsed_args["output_dir"]
    
    # Make the output_dir if it does not exist already
    # Note: This step might not work if the output directory is part of another non-existing directory,
    # which is why I prefer to make the output directory beforehand (see execute_distanceCalculation.sh)
    makeDirectory(output_dir)
    
    # Retrieve the features
    features = retrieveFeatures(selectedFeatures_path, selectionStrategy, treatment)
    
    _print("\n")
    
    # Carry out the UMAP reduction on default settings if the UMAP does not already exist
    # and export the UMAP
    output_path = "$(output_dir)/umapReducedProfile_$(cell)_$(selectionStrategy).csv"
    if isfile(output_path) == false
        reducedProfile = umapReduction(profile_path, cell, features)
        dataframe_description = "UMAP-reduced profile"
        exportDataFrame(output_path, reducedProfile; dataframe_description = dataframe_description)
    else
        _print("UMAP-reduced profile has been generated before. Using previously generated UMAP-reduced profile from:\n$(output_path)")
        reducedProfile = CSV.read(output_path, DataFrame)
    end
    
    _print("\n")
    
    # Carry out the distance calculation
    distance = calculateDistance(reducedProfile, referenceCompound)
    
    # Export the distance calculated
    output_path = "$(output_dir)/distance_$(cell)_$(selectionStrategy)_$(referenceCompound).csv"
    dataframe_description = "Distance calculated for vs. $(referenceCompound)"
    exportDataFrame(output_path, distance; dataframe_description = dataframe_description)
    
end

# # TEST VERSION for checking the mapping of parameters from command line
# function main()
#     parsed_args = parse_commandline()
#     println("Parsed args:")
#     for (arg,val) in parsed_args
#         println("  $arg  =>  $val")
#     end
# end

main()