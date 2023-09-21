#######################################
# Commands to run for specific use case
# Use case: GW015_006__full_cpa
#
# Written by: Amanda Ng R.H.
# Last updated on: 06 Sep 2023
#######################################

# Part A --> completed on 06 April 2023
# Note: It takes about 60 min per plate to complete part A (illumination correction).
parent_dir=/research/lab_winter/users/ang/projects/GW015_SuFEX_IMiDs/GW015_006__full_cpa
cell_line_dir_array=(c662_rko_wt c1141_rko_ko c1327_rko_oe)
a_path=/research/lab_winter/users/ang/cell_painting_assay_pipeline/01_feature_extraction/DEV_a_cpa_feature_extraction.sh
mode=non-test
transferlist=${parent_dir}/transferlist.csv
channel_annotation=${parent_dir}/channel_annotation.csv
plate_annotation=${parent_dir}/plate_annotation.csv
for dir in "${cell_line_dir_array[@]}"
do
  printf "${dir}\n"
  raw_data_dir=${parent_dir}/${dir}/raw_data_dir
  output_dir=${parent_dir}/${dir}/output_dir2
  sbatch ${a_path} ${mode} ${raw_data_dir} ${transferlist} ${channel_annotation} ${plate_annotation} ${output_dir}
done

# Part B --> started at 1544, 06 April 2023
# Note: It took about 2.0 h to prepare the images for segmentation and a further 1.5 hours to complete the segmentation for 20 batches of images (231 images per batch)
parent_dir=/research/lab_winter/users/ang/projects/GW015_SuFEX_IMiDs/GW015_006__full_cpa
cell_line_dir_array=(c662_rko_wt c1141_rko_ko c1327_rko_oe)
b_path=/research/lab_winter/users/ang/cell_painting_assay_pipeline/01_feature_extraction/b_cpa_feature_extraction.sh
for dir in "${cell_line_dir_array[@]}"
do
  printf "${dir}\n"
  output_dir=${parent_dir}/${dir}/output_dir2
  sbatch ${b_path} ${output_dir}
done

# Part C --> completed on 21 Mar 2023
# Note: The processing time took between 4.5 to 8.0 hours per batch of 231 sets of images
parent_dir=/research/lab_winter/users/ang/projects/GW015_SuFEX_IMiDs/GW015_006__full_cpa
cell_line_dir_array=(c662_rko_wt c1141_rko_ko c1327_rko_oe)
c_path=/research/lab_winter/users/ang/cell_painting_assay_pipeline/01_feature_extraction/c_cpa_feature_extraction.sh
for dir in "${cell_line_dir_array[@]}"
do
  printf "${dir}\n"
  output_dir=${parent_dir}/${dir}/output_dir2
  sbatch ${c_path} ${output_dir}
done

# Proceed to post-feature extraction pipeline by starting the jupyter lad IDE
input_dir=/research/lab_winter/users/ang
jupyterlab_script=/research/lab_winter/users/ang/jupyterlab_cluster_setup/jupyterlab_setup.sh
sbatch ${jupyterlab_script} ${input_dir}
log=/research/lab_winter/users/ang/jupyterlab_cluster_setup/logs/jupyter-notebook_"<JOB ID>".log
less ${log}
# Copy and paste the URL the log or you can set up an SSH tunnel (the instructions are also printed in the log)