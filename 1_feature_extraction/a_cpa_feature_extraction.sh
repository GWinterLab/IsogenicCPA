#!/bin/bash
#SBATCH --output=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_logs/a_cpa_feature_extraction_%j.log     #%j stands for unique job ID
#SBATCH --error=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_errors/a_cpa_feature_extraction_%j.err    #%j stands for unique job ID
#SBATCH --job-name="A_FeatureExtraction"
#SBATCH --partition=shortq    # job queue where the job is submitted to
#SBATCH --qos=shortq    # qos must match the partition
#SBATCH --nodes=1     # number of physical nodes
#SBATCH --ntasks=1    # 1 task
#SBATCH --cpus-per-task=10    # 1 task on 10 CPU
#SBATCH --time=3:00:00   # Job time is max 2 hours (it's cancelled if it's still executing after 3 hours)
#SBATCH --mem-per-cpu=8G    # using 8 Gb of memory
#SBATCH --mail-type=end   # send an email when this job ends
#SBATCH --mail-user=ang@cemm.at   # email your CeMM account

#############################################
# Author(s): Amanda Ng R.H.
# Created on: 10 Feb 2022
# Last updated on: 06 Sep 2023
# Batch processing status: Ready
# Overall status: Tested
# (using DEV_metadata_modules.py which handles the treatment assignment better)
#
# Master script A of the Cell Painting Assay (CPA) feature extraction pipeline which:
# 1. Generates the metadata file for cellprofiler-based feature extraction and illumination correction.
# 2. Calculate the illumination matrix per plate.
#############################################

/bin/hostname

###############################################################################
# Check if the user has provided the necessary variables for running the script
###############################################################################
if [ $# -ne 6 ]; then
  printf "\nERROR: User has not provided the necessary variables for running the script."
  printf "\nUsage: %s mode raw_data_dir channel_annotation transferlist output_parent_dir" ${0}
  printf "\nmode: Specify if the pipeline is being used in test or non-test mode"
  printf "\nraw_data_dir: Path to the parent directory containing the subdirectories of raw images. Each subdirectory is a plate/imaging run."
  printf "\ntransferlist: .csv file containing information on the name and concentration of the compound used for each well."
  printf "\nchannel_annotation: Path to the .csv file containing information on the channels."
  printf "\nplate_annotation: Path to the .csv file containing the information on the plates."
  printf "\noutput_parent_dir: Path to the parent directory for all the outputs from the feature extraction pipeline of the Cell Painting Assay."

  printf "\n\n~~~~~~~~~~ Inputs recieved ~~~~~~~~~~"
  printf "\nmode:\n${1}"
  printf "\nraw_data_dir:\n${2}"
  printf "\ntransferlist:\n${3}"
  printf "\nchannel_annotation:\n${4}"
  printf "\nplate_annotation:\n${5}"
  printf "\noutput_parent_dir:\n${6}"
  exit 1
fi

##################
# Assign variables
##################
# User inputs
mode="${1}"
raw_data_dir="${2}"
transferlist="${3}"
channel_annotation="${4}"
plate_annotation="${5}"
output_parent_dir="${6}"

# Batching variable
number_of_batches=20

# Parent directory for all scripts for the feature extraction
script_parent_dir=/research/lab_winter/users/ang/isogenicCPA_repo/1_feature_extraction

# Location of dependency scripts
dependency_scripts_dir=${script_parent_dir}/dependency_scripts

# Script/cppipe paths
metadata_modules=${dependency_scripts_dir}/metadata_modules.py
module_2_cppipe=${dependency_scripts_dir}/module_2.cppipe

# Module output parent directories
module_1_output=${output_parent_dir}/module_1_output
module_2_output=${output_parent_dir}/module_2_output
module_3_output=${output_parent_dir}/module_3_output

# Module 3 output file paths
module_3_metadata_path=${module_3_output}/module_3_metadata.csv
module_3subset_output_path=${module_3_output}/module_3subset_output.csv

# Array of directories that may not exist yet
make_dir_array=(${output_parent_dir} ${module_1_output} ${module_2_output} ${module_3_output})

date
start=$(date +%s)

# Record the user input/check the variable assignment
printf "\n### Running part A of the CPA feature extraction pipeline ###"
printf "\n~~~~~~~~~~ User input ~~~~~~~~~~"
printf "\mode: %s" ${mode}
printf "\nraw_data_dir: %s" ${raw_data_dir}
printf "\nchannel_annotation: %s" ${pretrained_model}
printf "\ntransferlist: %s" ${transferlist}
printf "\noutput_parent_dir: %s ${output_parent_dir}"
printf "\n"

###########################################
# Establish source for virtual environments
###########################################
source ~/miniconda3/etc/profile.d/conda.sh

#############################################################################
# Check if the channel and treatment annotation sheets are in the .csv format
#############################################################################
for path in ${channel_annotation} ${transferlist}
do
  if ! [[ "${path}" == *".csv" ]]; then
    printf "\nERROR: The following path provided should be a .csv file but is not:\n%s\nPlease check the path and/or convert the file into the .csv format.\n" ${path}
    exit 1;
  fi
done

#################################################
# Make any directory if they do not exist already
#################################################
for dir in "${make_dir_array[@]}"
do
  mkdir -p ${dir}
done

################################################################################
# MODULE 1 execution: generate metadata file containing the image paths and their associated metadata
################################################################################
# Activate the metadata_modules environment
conda activate metadata_modules

printf "\n### Module 1 output:"
python3 ${metadata_modules} --make_metadata1 ${raw_data_dir} ${channel_annotation} ${plate_annotation} ${transferlist} ${module_1_output}

# Check if module 1 has managed to generate the metadata sheet
# otherwise, exit the pipeline
compiled_file=${module_1_output}/module_1_output_compiled.csv
if ! [[ -f "${compiled_file}" ]] ; then
  printf "\nERROR in module 1 execution.\n"
  exit 1;
fi

# Deactivate the environment
conda deactivate

################################################################################
# MODULE 2 execution: generate the illumination correction matrix for each channel PER PLATE
################################################################################
# Activate the virtual environment specific for the use of CellProfiler v4.2.1
conda activate cellprofiler

printf "\n### Running module 2...\n"

for file in ${module_1_output}/*
  do
    if [[ ${file} != ${compiled_file} ]] ; then
      printf "${file}\n"
      cellprofiler -c -r -p ${module_2_cppipe} -o ${module_2_output} --data-file ${file}
    fi
  done

# Check if module 2 has managed to generate the metadata sheet
# otherwise, exit the pipeline
# if the output from ls ${module_2_output}/*.npy 1> results in an error (i.e. no npy files can be found, redirect the error message from standard error (2) to 1 which would result in the "else" situation)
if ls ${module_2_output}/*.npy 1> /dev/null 2>&1; then
  printf "\nIllumination correction matrix for each channel has been generated.\n"
else
  printf "\nERROR in module 2 execution.\n"
  exit 1;
fi

# Deactivate virtual environment
conda deactivate

################################################################################
# MODULE 3 execution: update the metadata file with the paths to the illumination correction matrices
################################################################################
# Activate metadata_modules environment
conda activate metadata_modules

printf "\n### Module 3 output:"
python3 ${metadata_modules} --update_metadata1 ${module_2_output} ${channel_annotation} ${compiled_file} ${module_3_metadata_path}

# Check if module 3 has managed to generate the metadata sheet
# otherwise, exit the pipeline
if [[ -f "${module_3_metadata_path}" ]] ; then
  printf "\nThe metadata sheet has been updated with the illumination correction matrix paths and file names.\n"
else
  printf "\nERROR in module 3 execution.\n"
  exit 1;
fi

# If test mode is triggered, make a subset of the metadata sheet for downstream testing of the pipeline
# Note: I only used the test mode when I was setting up the feature extraction pipeline.
if [[ "${mode}" == "test" ]]; then
  printf "\n### Test mode output:\nUser has indicated that the pipeline should be run in testing mode.\nMaking a subset of the metadata sheet for downstream usage in the pipeline...\n"
  python3 ${metadata_modules} --make_metadata1_subset ${module_3_metadata_path} ${module_3subset_output_path}

# Otherwise, split the metadata sheet into X equal parts
# where X is controlled by the number_of_batches
else
  printf "\n### Splitting the metadata file from module 3 into %s mini-metadata files for running the subsequent scripts in separate batches...\n" ${number_of_batches}
  python3 ${metadata_modules} --split_metadata1 ${module_3_metadata_path} ${module_3_output} ${number_of_batches}

  # Check if the metadata was split by checking if te first mini CSV file was generated
  if [[ -f "${module_3_output}/1_mini.csv" ]] ; then
    printf "\nThe metadata sheet has been split.\n"
  else
    printf "\nERROR in split_metadata1 execution.\n"
    exit 1;
  fi

fi

# Deactivate the environment
conda deactivate

#################
# Update the user
#################
printf "\n### COMPLETED: a_cpa_feature extraction.sh ###\n"
date
end=$(date +%s)
runtime=$((${end}-${start}))
printf "\nTime taken: ${runtime} s"
printf "\nPlease proceed to running b_cpa_feature_extraction.sh\n"

### END OF MASTER SCRIPT A OF THE CPA FEATURE EXTRACTION PIPELINE ###
