#!/bin/bash
#SBATCH --output=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_logs/c_cpa_feature_extraction_%j.log     #%j stands for unique job ID
#SBATCH --error=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_errors/c_cpa_feature_extraction_%j.err    #%j stands for unique job ID
#SBATCH --job-name="C_FeatureExtraction"
#SBATCH --partition=shortq    # job queue where the job is submitted to
#SBATCH --qos=shortq    # qos must match the partition
#SBATCH --nodes=1     # number of physical nodes
#SBATCH --ntasks=1    # 1 task
#SBATCH --cpus-per-task=5    # 1 task on 5 CPU
#SBATCH --time=1:00:00   # Job time is max 1 hour (it's cancelled if it's still executing after 1 hour)
#SBATCH --mem=8000    # using 8 Gb of memory
#SBATCH --mail-type=end   # send an email when this job ends
#SBATCH --mail-user=ang@cemm.at   # email your CeMM account

#############################################
# Author(s): Amanda Ng R.H.
# Created on: 14 Feb 2022
# Last updated on: 06 Sep 2023
# Batch processing status: Ready
# Overall status: Tested and ready for use
#
# Master script C of the Cell Painting Assay (CPA) feature extraction pipeline which:
# 1. Submits the scripts for updating the metadata file with the filename and path details of the segmentation masks, and runs the feature extraction cppipe.
#############################################

/bin/hostname

###############################################################################
# Check if the user has provided the necessary variables for running the script
###############################################################################
if [ $# -ne 1 ]; then
  printf "\nERROR: User has not provided the necessary variables for running the script."
  printf "\nUsage: %s output_parent_dir" ${0}
  printf "\noutput_parent_dir: Path to the parent directory for all the outputs from the feature extraction pipeline of the Cell Painting Assay."
  exit 1
fi

###########################################
# Establish source for virtual environments
###########################################
source ~/miniconda3/etc/profile.d/conda.sh

##################
# Assign variables
##################
# User input
output_parent_dir="${1}"

# Parent directory for all scripts for the feature extraction
script_parent_dir=/research/lab_winter/users/ang/isogenicCPA_repo/1_feature_extraction

# Location of dependency scripts
dependency_scripts_dir=${script_parent_dir}/dependency_scripts

# Script/cppipe paths
metadata_modules=${dependency_scripts_dir}/metadata_modules.py
module_6b_7b_batch=${dependency_scripts_dir}/module_6b_7b_batch.sh

# Module output/relevant input directories
module_3_output=${output_parent_dir}/module_3_output
module_5b_output=${output_parent_dir}/module_5b_output
module_6b_output=${output_parent_dir}/module_6b_output
module_7b_output=${output_parent_dir}/module_7b_output

# Temporary directory used for the inputs in the whole cell segmentation in the previous master script (b_cpa_feature_extraction)
temp_whole_cell_segmentation_input=${output_parent_dir}/temp_whole_cell_segmentation_input

# Module 6b input file
module_3_metadata_path=${module_3_output}/module_3_metadata.csv

# Array of directories that may not exist yet
make_dir_array=(${module_6b_output} ${module_7b_output})

date
start=$(date +%s)

# Record the user input/check the variable assignment
printf "\n### Running part C of the CPA feature extraction pipeline ###"
printf "\n~~~~~~~~~~ User input ~~~~~~~~~~"
printf "\noutput_parent_dir: %s" ${output_parent_dir}
printf "\n"

################################################################################
# Check if the expected module 5b outputs are present before proceeding with c_cpa_feature_extraction
################################################################################
# Activate metadata_modules environment
conda activate metadata_modules

status=$(python3 ${metadata_modules} --check_script_b_completion ${module_3_metadata_path} ${module_5b_output})
echo ${status}
if [[ ${status} == *"Overall: Fail"* ]] ; then
  printf "\nERROR: b_cpa_feature_extraction/module 5b (nuclei and whole cell segmentation) has not run to completion or there is an error. Please check.\n"
  exit 1;
else
  printf "\n### b_cpa_feature_extraction has run to completion. Continuing with running c_cpa_feature_extraction.\n"
  if [ -d ${temp_whole_cell_segmentation_input} ] ; then
    printf "\nDeleting the temporary folder used for whole cell segmentation in b_cpa_feature_extraction..."
    rm -r ${temp_whole_cell_segmentation_input}
    multiple_models=true
  else
    multiple_models=false
  fi
fi

# Deactivate the virtual environment
conda deactivate

#################################################
# Make any directory if they do not exist already
#################################################
for dir in "${make_dir_array[@]}"
do
  mkdir -p ${dir}
done

########################################
# Modules 6b and 7b execution in batches
########################################
# Make the module_6b_7b_batch executable
dos2unix ${module_6b_7b_batch}
chmod +x ${module_6b_7b_batch}

printf "\n### SUBMITTING: Scripts for metadata update (module 6b) and feature extraction (module 7b)\n"

#==============================================#
# Function for executing modules 6b and 7b #
function execute_modules6b7b {

  # Arguments for the function
  module_6b_output=${1}
  module_7b_output=${2}
  whole_cell_segmentation_dir=${3}

  # Print the agruments recieved for checking purposes
  printf "\n~~~~~~~~~~ User input ~~~~~~~~~~"
  printf "\nmodule_6b_output: %s" ${module_6b_output}
  printf "\nmodule_7b_output: %s" ${module_7b_output}
  printf "\nwhole_cell_segmentation_dir: %s" ${whole_cell_segmentation_dir}

  # Counter for tracking and naming the batches
  counter=0

  # For each mini metadata file i.e. each batch
  for mini_metadata in ${module_3_output}/*mini.csv
  do
    counter=$((${counter}+1))

    # Full path to the metadata files generated by module 6b for use in module 7b
    updated_mini_metadata=${module_6b_output}/${counter}_mini.csv

    # Make the folder for storing the information and output from module 7b for each batch
    temp_module_7b_output=${module_7b_output}/${counter}
    mkdir -p ${temp_module_7b_output}

    # Run modules 6b and 7b
    sbatch ${module_6b_7b_batch} ${module_5b_output} ${mini_metadata} ${updated_mini_metadata} ${temp_module_7b_output} ${whole_cell_segmentation_dir}

  done
}
# End of function for executing modules 6b and 7b
#==============================================#

# Run modules 6b and 7b in different manners depending on whether multiple whole cell segmentation models were used or not
if [ "${multiple_models}" = "true" ] ; then
  printf "DETECTED: Multiple whole cell segmentation models were used.\n"

  for subdir in $(ls ${module_5b_output}/whole_cell_segmentation/)
  do

    # Make subdirectories in modules 6b and 7b for each whole cell segmentation model used
    module_6b_subdir=${module_6b_output}/${subdir}
    mkdir -p ${module_6b_subdir}
    module_7b_subdir=${module_7b_output}/${subdir}
    mkdir -p ${module_7b_subdir}

    # Define the whole cell segmentation model-specific directory
    whole_cell_segmentation_dir=${module_5b_output}/whole_cell_segmentation/${subdir}

    # Run modules 6b and 7b
    execute_modules6b7b ${module_6b_subdir} ${module_7b_subdir} ${whole_cell_segmentation_dir}
  done

else
  printf "DETECTED: Only one whole cell segmentation model was used used.\n"
  whole_cell_segmentation_dir=${module_5b_output}/whole_cell_segmentation
  execute_modules6b7b ${module_6b_output} ${module_7b_output} ${whole_cell_segmentation_dir}
fi

printf "\n### SUBMITTED: Scripts for metadata update (module 6b) and feature extraction (module 7b)\n"

#################
# Update the user
#################
printf "\n### COMPLETED: c_cpa_feature_extraction.sh ###\nPlease wait for the scripts submitted above to be completed before proceeding to d_cpa_feature_extraction.sh."
date
end=$(date +%s)
runtime=$((${end}-${start}))
printf "\nTime taken: ${runtime} s"

### END OF MASTER SCRIPT C OF THE CPA FEATURE EXTRACTION PIPELINE ###
