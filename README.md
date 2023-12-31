# Isogenic CPA
The isogenic Cell Painting Assay pipeline uses cell morphological features from isogenic cell lines treated with the same compounds. The computational pipeline and tools provided here can be used to predict and characterize multi-specific drugs by comparing the morphological features changed across isogenic cell lines expressing different levels of a presenter protein. This pipeline was used in the manuscript titled "Discovery of molecular glue degraders via isogenic morphological profiling" for the presenter protein, CRBN. [Check out the paper](https://doi.org/10.1021/acschembio.3c00598) 😊

<img src="https://github.com/GWinterLab/IsogenicCPA/assets/91612461/15eaf164-3c57-478f-97b3-38b43e3c7a77" width="500">

## 0| Virtual environment requirements
The files here use multiple virtual environments. The requirements for setting up these environments can be found in [0_virtual_environment_requirements](https://github.com/GWinterLab/IsogenicCPA/tree/3f7198ac1248695ab301744906474c780cae933f/0_virtual_environments_required).

## 1| Feature extraction
The feature extraction procedure used here is adapted from the conventional feature extraction pipeline used by the [JUMP-CPA Consortium](https://github.com/broadinstitute/imaging-platform-pipelines/tree/master/JUMP_production). The key changes are:
  - Exclusion of the RNA channel
  - Usage of [cellpose](https://github.com/MouseLand/cellpose) (deep learning-based segmentation)

To run the feature extraction procedure, you can find the instructions and a brief explanation in [isogenicCPApipeline.pdf](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/isogenicCPApipeline.pdf). You can also find an example of the commands used in [1_commands_to_run.sh](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/1_commands_to_run.sh).

## 2| Data pre-processing (QC, Assembly of morphological profiles and Feature selection)
The instructions and details of the data pre-processing can be found in [isogenicCPApipeline.pdf](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/isogenicCPApipeline.pdf). The functions used are from custom modules in [2_feature_analysis](https://github.com/GWinterLab/IsogenicCPA/tree/3f7198ac1248695ab301744906474c780cae933f/2_feature_analysis) and the asociated notebook ([2_QC_ProfileAssembly_FeatureSelection.ipynb](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/2_QC_ProfileAssembly_FeatureSelection.ipynb)) provides an example of the function usage .

## 3| Profile calculations
In the manuscript, calculations were done with the morphological profiles in two ways:
 - Using one isogenic cell line
   <br>Calculating the morphological perturbation strength of compounds compared to DMSO (basal morphology) using the commands in [3a_execute_distanceCalculation.sh](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/3a_execute_distanceCalculation.sh) which is based off of [BioProfiling.jl](https://github.com/menchelab/BioProfiling.jl).
 - Using all three isogenic cell lines
   <br>Calculating the likelihood that a compound is dependent on the presenter protein (expressed at varying levels across the isogenic cell lines). The relevant functions used are from custom modules in [2_feature_analysis](https://github.com/GWinterLab/IsogenicCPA/tree/3f7198ac1248695ab301744906474c780cae933f/2_feature_analysis) and an example of their usage can be found in [3b_ProfileCalculations.ipynb](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/3b_ProfileCalculations.ipynb).

## 4| Profile interpretations
Finally, once the profile calculations are completed, interpretations can be made. Examples of the interpretations (with example outputs) can be found in [4_ProfileInterpretation.ipynb](https://github.com/GWinterLab/IsogenicCPA/blob/3f7198ac1248695ab301744906474c780cae933f/4_ProfileInterpretation.ipynb).
