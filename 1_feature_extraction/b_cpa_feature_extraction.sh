#!/bin/bash
#SBATCH --output=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_logs/b_cpa_feature_extraction_%j.log     #%j stands for unique job ID
#SBATCH --error=/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_errors/b_cpa_feature_extraction_%j.err    #%j stands for unique job ID
#SBATCH --job-name="B_FeatureExtraction"
#SBATCH --partition=shortq    # job queue where the job is submitted to
#SBATCH --qos=shortq    # qos must match the partition
#SBATCH --nodes=1     # number of physical nodes
#SBATCH --ntasks=1    # 1 task
#SBATCH --cpus-per-task=10    # 1 task on 10 CPU
#SBATCH --time=3:00:00   # Job time is max 3 hours (it's cancelled if it's still executing after 3 hours)
#SBATCH --mem-per-cpu=8G    # using 8 Gb of memory
#SBATCH --mail-type=end   # send an email when this job ends
#SBATCH --mail-user=ang@cemm.at   # email your CeMM account

#############################################
# Author(s): Amanda Ng R.H.
# Created on: 14 Feb 2022
# Last updated on: 06 Sep 2023
# Batch processing status: Ready
# Overall status: Tested and ready for use
#
# Master script B of the Cell Paintinf Assay (CPA) pipeline which:
# 1. Prepares the images for nuclei and whole cell segmentation by cellpose.
# 2. Submits the scripts (module_5b.sh) for nuclei and whole cell segmentation by cellpose.
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

# Script/cppipe path(s)
module_4b_cppipe=${dependency_scripts_dir}/module_4b.cppipe
module_5b_script=${dependency_scripts_dir}/module_5b.sh

# Relevant module input and output directories
module_3_output=${output_parent_dir}/module_3_output
module_4b_output=${output_parent_dir}/module_4b_output
module_5b_output=${output_parent_dir}/module_5b_output

# Module 5b output sub-directories
nuclei_output=${module_5b_output}/nuclei_segmentation
whole_cell_output=${module_5b_output}/whole_cell_segmentation

# Cellpose model(s) to use for whole cell segmentation
whole_cell_model=cyto
whole_cell_model_array=(${whole_cell_model})
# Note: You can other models here if you have trained them.
# e.g.:
# whole_cell_model=cyto
# whole_cell_model_2="${script_parent_dir}/trained_whole_cell_model/20220313_1051__generalist_cyto_regular"
# whole_cell_model_array=(${whole_cell_model} ${whole_cell_model_2})

date
start=$(date +%s)

# Record the user input/check the variable assignment
printf "\n### Running part B of the CPA feature extraction pipeline ###"
printf "\n~~~~~~~~~~ User input ~~~~~~~~~~"
printf "\noutput_parent_dir: %s" ${output_parent_dir}
printf "\n"

#################################################
# Make any directory if they do not exist already
#################################################
make_dir_array=(${module_4b_output} ${module_5b_output} ${nuclei_output} ${whole_cell_output} ${whole_cell_output_2})
for dir in "${make_dir_array[@]}"
do
  mkdir -p ${dir}
done

#############################################################################
# MODULE 4b execution: prep the images for nuclei and whole cell segmentation
#############################################################################
# Activate cellprofiler environment
conda activate cellprofiler

printf "\n### Module 4b output:"
counter=0
for mini_metadata in ${module_3_output}/*_mini.csv
do
  counter=$((${counter}+1))
  output_dir=${module_4b_output}/${counter}
  mkdir -p ${output_dir}
  cellprofiler -c -r -p ${module_4b_cppipe} -o ${output_dir} --data-file ${mini_metadata}

  # Check if module 4b was executed successfully
  # else, exit the pipeline
  if ls ${output_dir}/*.tiff 1> /dev/null 2>&1; then

    # Move the tiff images prepped for segmentation
    for string in "dna" "dna_agp"
    do
      output_subdir=${output_dir}_${string}
      mkdir -p ${output_subdir}
      mv ${output_dir}/*${string}.tiff ${output_subdir}
    done

    # Delete the temporary output directory for module 4
    rmdir ${output_dir}

    printf "\nThe images for segmentation have been prepared and shifted to the intended locations."

  else
    printf "\nERROR in module 4b execution.\n"
    exit 1;

  fi

done

# deactivate the virtual environment
conda deactivate

#####################################################################
# MODULE 5b execution: nuclei and whole cell segmentation via CellPose
#####################################################################
# Make the module 5b scripts executable
dos2unix ${module_5b_script}
chmod +x ${module_5b_script}

# Run the nuclei segmentation for each subset of prepped nuclei segmentation images
printf "\n### SUBMITTING: Scripts for nuclei segmentation\n"
for input_dir in ${module_4b_output}/*_dna
do
  sbatch ${module_5b_script} ${input_dir} ${nuclei_output} nuclei
done
printf "\n### SUBMITTED: Scripts for nuclei segmentation\n"


printf "\n### SUBMITTING: Scripts for whole cell segmentation\n"

### Single segmentation model ###
if [ "${#whole_cell_model_array[@]}" -eq "1" ] ; then
  for input_dir in ${module_4b_output}/*_dna_agp
  do
    sbatch ${module_5b_script} ${input_dir} ${whole_cell_output} ${whole_cell_model_array[0]}
  done

### Multiple segmentation models ###
else
  for model in "${whole_cell_model_array[@]}"
  do
    if [[ -f ${model} ]]; then
      model_shortname="$(basename -- ${model})"
    else
      model_shortname=${model}
    fi

    # Make a copy of the composite images prepped for segmentation temporarily so that whole cell segmentation using multiple models can be done simultaneously without any conflicts
    temp_composite_dir=${output_parent_dir}/temp_whole_cell_segmentation_input/${model_shortname}
    mkdir -p ${temp_composite_dir}
    cp -r -p ${module_4b_output}/*_dna_agp ${temp_composite_dir}

    # Make the directory for the segmentation masks for each model
    output_dir=${whole_cell_output}/${model_shortname}

    # Carry out the whole cell segmentation in batches
    for input_dir in ${temp_composite_dir}/*_dna_agp
    do
      sbatch ${module_5b_script} ${input_dir} ${output_dir} ${model}
    done
  done


fi

printf "\n### SUBMITTED: Scripts for whole cell segmentation\n"

# #################
# # DEV RUN/TESTING
# #################
# input_dir=${module_4b_output}/1_dna
# output_dir=${output_parent_dir}/test_module_4b_output_2
# mkdir -p ${output_dir}
# sbatch ${module_5b_script} ${input_dir} ${output_dir} nuclei

#################
# Update the user
#################
printf "\n### COMPLETED: b_cpa_feature extraction.sh ###\n"
date
end=$(date +%s)
runtime=$((${end}-${start}))
printf "\nTime taken: ${runtime} s"
printf "\nPlease wait for the segmentation to run to completion before running c_cpa_feature_extraction.sh\n"

### END OF MASTER SCRIPT B OF THE CPA FEATURE EXTRACTION PIPELINE ###
