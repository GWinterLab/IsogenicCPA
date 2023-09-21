#!/bin/bash
#SBATCH --output=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_logs/module_6b_7b_batch_%j.log    #%j stands for unique job ID
#SBATCH --error=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_errors/module_6b_7b_batch_%j.err    #%j stands for unique job ID
#SBATCH --job-name="module_6b_7b_batch"
#SBATCH --partition=shortq    # job queue where the job is submitted to
#SBATCH --qos=shortq    # qos must match the partition
#SBATCH --nodes=1     # number of physical nodes
#SBATCH --ntasks=1    # 1 task
#SBATCH --cpus-per-task=5    # 1 task on 5 CPU
#SBATCH --time=10:00:00   # Job time is max 10 hours (it's cancelled if it's still executing after 10 hours)
#SBATCH --mem=8000        # using 8 Gb of memory
#SBATCH --mail-type=end   # send an email when this job ends
#SBATCH --mail-user=ang@cemm.at   # email your CeMM account

#############################################
# Author(s): Amanda Ng R.H.
# Created on: 14 Feb 2022
# Last updated on: 06 Sep 2023
#
# Script for running modules 6b and 7b sequentially.
# Used by c_cpa_feature_extraction.sh to batch process feature extraction.
#############################################

/bin/hostname

###########################################
# Establish source for virtual environments
###########################################
source ~/miniconda3/etc/profile.d/conda.sh

#####################
# Variable assignment
#####################
# User input
module_5b_output="${1}"
mini_metadata="${2}"
updated_mini_metdata="${3}"
module_7b_output="${4}"
whole_cell_segmentation_dir="${5}"

# Parent directory for all scripts for the feature extraction
script_parent_dir=/research/lab_winter/users/ang/isogenicCPA_repo/1_feature_extraction

# Location of dependency scripts
dependency_scripts_dir=${script_parent_dir}/dependency_scripts

# Additional variable(s)
nuclei_segmentation_dir=${module_5b_output}/nuclei_segmentation

# Script/cppipe paths
metadata_modules=${dependency_scripts_dir}/metadata_modules.py
module_7b_cppipe=${dependency_scripts_dir}/module_7b.cppipe

#################
# Update the user
#################
printf "\n### STARTING: Modules 6b and 7b ###\n"
printf "\n~~~~~~~~~~ User input ~~~~~~~~~~"
printf "\nmodule_5b_output: %s" ${module_5b_output}
printf "\nmini_metadata: %s" ${mini_metadata}
printf "\nupdated_mini_metdata: %s" ${updated_mini_metdata}
printf "\nmodule_7b_output: %s" ${module_7b_output}
printf "\nwhole_cell_segmentation_dir: %s" ${whole_cell_segmentation_dir}
printf "\n"
date
start=$(date +%s)

################################################################################
# MODULE 6b execution: update the metadata file with the paths to the nuclei and whole cell segmentation masks
################################################################################
# Activate metadata_modules environment
conda activate metadata_modules

printf "\n### Module 6b output:"
python3 ${metadata_modules} --update_metadata2 ${nuclei_segmentation_dir} ${whole_cell_segmentation_dir} ${mini_metadata} ${updated_mini_metdata}

# Check if module 6 was executed successfully
if ! [[ -f "${updated_mini_metdata}" ]] ; then
  printf "\nERROR in module 6 execution.\n"
  exit 1;
fi

# Deactivate the environment
conda deactivate

################################################################################
# MODULE 7b execution: illumination correction of all images and feature extraction
################################################################################
# Activate the cellprofiler environment
conda activate cellprofiler

printf "\n### Running module 7b...\n"
cellprofiler -c -r -p ${module_7b_cppipe} -o ${module_7b_output} --data-file ${updated_mini_metdata}
# 15 Feb 2022 dev note: I am running into issues with the FilterObjects module, where an error keeps popping up and I can't figure out why. I've turned my attention to prepping Master Script B for batching temporarily and then to improving the segmentation before returning to working on this issue.

# Check if module 7b was executed successfully
# by checking if there are any .png files in the module_7b_output, which corresponds to the outlines of nuclei and/or whole cells
outlines=$(ls ${module_7b_output}/*.png 2> /dev/null | wc -l)
if [[ "${outlines}" != "0" ]]; then
  printf "\nFeature extraction has been completed.\n### COMPLETED: module_6b_7b_batch has run to completion successfully. ###\n"
else
  printf "\nERROR in module 7b execuion.\n"
  exit 1;
fi

# Deactivate the environment
conda deactivate

#################
# Update the user
#################
date
end=$(date +%s)
runtime=$((${end}-${start}))
printf "\nTime taken: ${runtime} s"
