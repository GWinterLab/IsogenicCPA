#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --job-name='segmentation_of_cell_images'
#SBATCH --output='/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_logs/segmentation_of_cell_images_%j.log'
#SBATCH --error='/research/lab_winter/users/ang/cell_painting_assay_pipeline/temp_errors/segmentation_of_cell_images_%j.err'
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH -m block
#SBATCH --mail-type=end   # send an email when this job ends
#SBATCH --mail-user=ang@cemm.at   # email your CeMM account

#############################################
# Author(s): Amanda Ng R.H.
# Created on: 14 Feb 2022
# Last updated on: 06 Sep 2023
#
# Master script C connecting modules 6b and 8b in the cpa_feature_extraction pipeline, which is dependent on the modules/scripts in the dependency_scripts folder.
#############################################

/bin/hostname

################################################################################
# Activate the virtual environment specific for the segmentation carried out here with cellpose
################################################################################
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cellpose

#####################
# Variable assignment
#####################
# user input
input_dir="${1}"
output_dir="${2}"
pretrained_model="${3}"

#############################
# Segmentation of cell images
#############################
# check if the user input was mapped correctly
printf "\n### STARTING: Segmentation of cell images ###"
printf "\n~~~~~~~~~~ User input ~~~~~~~~~~"
printf "\ninput_dir: %s" ${input_dir}
printf "\noutput_dir: %s" ${output_dir}
printf "\npretrained_model: %s" ${pretrained_model}
printf "\n"

# timestamp and noting start time of script execution
date
start=$(date +%s)

# make output directory if it does not exist already
mkdir -p ${output_dir}


if [[ "${pretrained_model}" == "nuclei" ]]; then
  # segment nuclei from the illumination corrected DNA channel images (which is greyscale)
  python -m cellpose --dir ${input_dir} --pretrained_model ${pretrained_model} --diameter 0. --save_png --chan 0 --exclude_on_edges

else
  # segment images where the cytoplasm is green and the nucleus is blue
  # these images are composites comprising the DNA (blue) and AGP (green) channels which
  # have been illumination corrected in module 4b
  # target to produce: x.png and x.npy
  python -m cellpose --dir ${input_dir} --pretrained_model ${pretrained_model} --diameter 0. --save_png --chan 2 --chan2 3 --exclude_on_edges

fi

# delete the .npy files
rm ${input_dir}/*.npy

############################################################
# Transfer the masks output to the intended output directory
############################################################
printf '\n Transferring the segmentation masks to the specified output location...'
mv ${input_dir}/*_cp_masks.png ${output_dir}
printf "\n### COMPLETED: Segmentation of cell images ###\n"

# timestamp and noting the end time and run time for the script
date
end=$(date +%s)
runtime=$((${end}-${start}))
printf "\nTime taken: ${runtime} s"

####################################
# Deactivate the virtual environment
####################################
conda deactivate
