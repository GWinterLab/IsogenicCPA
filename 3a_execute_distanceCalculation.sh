#!/bin/bash
#SBATCH --output=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_logs/distanceCalculation_%j.log     #%j stands for unique job ID
#SBATCH --error=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_errors/distanceCalculation_%j.err    #%j stands for unique job ID
#SBATCH --job-name="distanceCalculation"
#SBATCH --partition=shortq    # job queue where the job is submitted to
#SBATCH --qos=shortq    # qos must match the partition
#SBATCH --nodes=1     # number of physical nodes
#SBATCH --ntasks=1    # 1 task
#SBATCH --cpus-per-task=10    # 1 task on 10 CPU
#SBATCH --time=12:00:00   # Job time is max 12 hours (it's cancelled if it's still executing after 12 hours)
#SBATCH --mem=8000    # using 8 Gb of memory
#SBATCH --mail-type=end   # send an email when this job ends
#SBATCH --mail-user=ang@cemm.at   # email your CeMM account

#############################################
# Author(s): Amanda Ng R.H.
# Created on: 18 Apr 2023
# Last updated on: 19 Apr 2023 (adjusted the sbatch settings)
# Overall status: Tested and used on 19 Apr 2023
#
# Script for executing distance caculation using distance_calculation.jl.
#############################################

/bin/hostname

date
start=$(date +%s)

##################
# Assign variables
##################
# User inputs
parent_dir="/research/lab_winter/users/ang/projects/GW015_SuFEX_IMiDs/GW015_006__full_cpa"
selectedFeatures_path="${parent_dir}/cleanPostFeatureExtraction/output_dir/2_ProfileAssemblyAndFeatureSelection_output/featuresSelected_output.csv"
selectionStrategy="global_c662_rko_wt"
treatment="all"
profile_path="${parent_dir}/cleanPostFeatureExtraction/output_dir/2_ProfileAssemblyAndFeatureSelection_output/concatenated__profile.csv"
cell="c662_rko_wt"
referenceCompound="DMSO"
output_dir="${parent_dir}/cleanPostFeatureExtraction/output_dir/3_ProfileInterpretation_output/BioProfiling_output"

# Location of the julia script for distance calculation
dependency_script=/research/lab_winter/users/ang/cell_painting_assay_pipeline/02_feature_analysis/distanceCalculation.jl

# Make the directory for the outputs from the distance calculation
if [[ ! -e ${output_dir} ]]; then
    mkdir -p ${output_dir}
fi

##################################
# Execute the distance calculation
##################################
# Establish source for virtual environments
source ~/miniconda3/etc/profile.d/conda.sh

# Activate conda environment required for distance calculation
# This same environment has the necessary packages for executing the distance calculation
conda activate juliaEnvironment

julia ${dependency_script} ${selectedFeatures_path} ${selectionStrategy} ${treatment} ${profile_path} ${cell} ${referenceCompound} ${output_dir}

# Close the conda environment
conda deactivate

#################
# Update the user
#################
# Name of output files expected from the distanceCalculation.jl
umapFileName="${output_dir}/umapReducedProfile_${cell}_${selectionStrategy}.csv"
distanceFileName="${output_dir}/distance_${cell}_${selectionStrategy}.csv"
output_files=(umapFileName distanceFileName)

# Check if the files have been exported
for file in "${output_files[@}]"
do
    if [ -f "${file}" ]; then
        printf "PASS/Expected output file from distanceCalculation.jl found:\n${file}\n\n"
    else
        printf "FAIL: The following file was NOT found:\n${file}\n\n"
        exit 1
    fi
done

# Update the user that the distance calculation and UMAP reduction has been completed
date
end=$(date +%s)
runtime=$((${end}-${start}))
printf "\nTime taken: ${runtime} s"

### END OF SCRIPT ###