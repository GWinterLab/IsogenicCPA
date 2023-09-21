#############################################
# !! WORKING VERSION !!
# Author(s): Amanda Ng R.H.
# Created on: 31 Mar 2023
# Last updated on: 21 Sep 2023 (added the functions for comparing treatment-centric features)
# Documentation status:
# In progress
# See https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html
# for more information on the Sphinx documentation style.
#
# Functions for profile interpretation.
#############################################
# Import packages
import os
import warnings

import pandas as pd

import seaborn as sns
from matplotlib import pylab as plt
from matplotlib.patches import Patch
import distinctipy

from numpy.linalg import norm
import statistics
import numpy as np
from scipy.stats import mannwhitneyu
import math

# Function for calculating the induction score
def calculateInductionScores(
    selectedFeatures,
    profile,
    treatmentsToUse = "all"
    
):
    """
    Function for calculating the induction score for each image per treatment using the associated
    treatment-centric features. These features were found to correlate with the CRBN status of cell lines
    for the given treatment. The induction score aggregates the robust Z scores of all features per image
    per treatment, effectively representing the magnitude of change in morphology toward a specific
    state of cell morphology.
    
    Mathematical(?) considerations:
     1. The aggregate score of all features would be high just by having more features. To counteract this
        situation, the aggregate score is averaged over all features to give the induction score.
     2. The features could correlate with CRBN status directly OR inversely, which means I should NOT
        aggregate the measurements directly. If I carried out direct aggregation, directly and inversely
        correlated measurements could end up negating each other. In other words, I should aggregate
        the absolute measurements.
        
    :param selectedFeatures: Dataframe with the features selected using different strategies.
    :type selectedFeatures: pd.DataFrame
    :param profile: Dataframe with all the measurements made per feature per treatment per cell line.
    :type profile: pd.DataFrame
    :param treatmentsToUse: If set to "all", the induction score per image will be calculated for all
        treatments in the profile and with at least one treatment-centric feature. Else, the user can
        also provide a list of treatments to restrict the induction score calculation to. Defaults to
        "all".
    :type treatmentsToUse: list or str
    :return inductionScores: Dataframe with the induction scores calculated along with some metadata
        i.e. "Metadata_Treatment", "Metadata_Cell", "Metadata_Plate", "Metadata_Well", "UpdatedImageNumber"
    :rtype inductionScores: pd.DataFrame
    """
    # Generate a dictionary for mapping the treatment to its features selected
    treatment2features_dict = dict()
    for i in selectedFeatures.index:
        treatment = selectedFeatures.loc[i, "Treatment"]
        features = selectedFeatures.loc[i, "Features_Selected"]
        try:
            features = features.split(" ")
            treatment2features_dict[treatment] = features
        except:
            print(f"No features selected: {treatment}\n")
    
    if treatmentsToUse == "all":
        treatmentsToUse = list(treatment2features_dict.keys())
        
    # Define the metadata columns to keep
    metadataToKeep = [
        "Metadata_Treatment",
        "Metadata_Cell",
        "Metadata_Plate",
        "Metadata_Well",
        "UpdatedImageNumber"
    ]
    
    # For each treatment
    inductionScores = []
    for treatment in treatmentsToUse:
        if treatment in treatment2features_dict:
        
            # Retrieve the features selected for the treatment
            features = treatment2features_dict[treatment]

            # Trim down the profile to just the treatment
            treatmentProfile = profile[profile["Metadata_Treatment"] == treatment]
            if len(treatmentProfile) > 0:

                # Copy the metadata from the profile to the inductionScore
                inductionScore = treatmentProfile[metadataToKeep].copy()

                # Trim down the treatmentProfile further to just the features selected
                treatmentProfile = treatmentProfile[features]

                # Calculate the induction score for each image
                # Note: The way the induction score is calculated assumes that the direction of feature change
                # is true for each observation i.e. the direction for all observations is NOT considered.
                # In theory, you could consider the overall direction of change (positive or negative).
                # I did not have enough control compounds to check if this consideration would improve the
                # prediction of CRBN dependency, but a future user could take this consideration into account.
                scores = []
                for i in treatmentProfile.index:
                    score = 0
                    for feature in features:
                        measurement = treatmentProfile.loc[i, feature]
                        score += abs(measurement)
                    score = score/len(features)
                    scores.append(score)

                # Append the induction score to the treatmentProfile
                inductionScore["InductionScore"] = scores

                # Calculate the normalized induction score (min-max normalization)
                # so that the scores can be compared across treatments
                minScore = min(scores)
                maxScore = max(scores)
                normalizedScores = []
                for score in scores:
                    normalizedScore = (score - minScore)/(maxScore - minScore)
                    normalizedScores.append(normalizedScore)
                inductionScore["NormalizedInductionScore"] = normalizedScores

                # Store the induction score calculated as separate dataframes
                inductionScores.append(inductionScore)
            
        else:
            print(f"Treatment not in profile or has no features selected: {treatment}\n")
        
    # Concatenate all the dataframes with the induction score
    inductionScores = pd.concat(inductionScores)
    
    return(inductionScores)


# Function for comparing induction scores between KO and OE
def compareInduction(
    inductionScores,
    ko = "c1141_rko_ko",
    oe = "c1327_rko_oe"
):
    """
    Calculate the difference in the distribution of induction scores observed in KO and OE cells using
    the Mann-Whiteny U test. In this case, I expect the OE to have a distrbution that has higher induction
    scores in general than the KO, so I use a right tail test.
    
    The Mann-Whitney U test returns a p-value and the U statistic. The U statistic is dependent on the
    sample sizes of the distributions being compared, where the:
        maximum U = number of samples from KO * number of samples from OE
    I, thus, converted the U statistic into a fraction of the maximum U, which I call the corrected U
    that ranges from 0 to 1.
    
    :param inductionScores: Dataframe with the induction scores calculated along with some metadata
        i.e. "Metadata_Treatment", "Metadata_Cell", "Metadata_Plate", "Metadata_Well", "UpdatedImageNumber"
    :type inductionScores: pd.DataFrame
    :param ko: Name of the cell line in the inductionScores that corresponds to the KO
        (e.g. "c1141_rko_ko").
    :type ko: str
    :param oe: Name of the cell line in the inductionScores that corresponds to the OE
        (e.g. "c1327_rko_oe").
    :type oe: str
    :return output_df: Dataframe with the p-value and corrected U values calculated per treatment.
    :rtype output_df: pd.DataFrame
    """
    # Tirm down the induction scores to just those from RKO CRBN KO and RKO CRBN OE
    inductionScores = inductionScores[inductionScores["Metadata_Cell"].isin([ko, oe])]
    
    # Retrieve the set of treatments in the inductionScores
    treatments = set(inductionScores["Metadata_Treatment"].tolist())
    
    # Make a dictionary for storing the calculated values
    output_dict = dict()
    for key in ["Metadata_Treatment", "Corrected_u", "p-value"]:
        output_dict[key] = []
    
    # For each treatment
    for treatment in treatments:
        
        # Trim down the induction score to just those pertaining to the treatment
        treatmentScore = inductionScores[inductionScores["Metadata_Treatment"] == treatment]
        
        # Retrieve the scores corresponding to RKO CRBN KO
        ko_scores = treatmentScore[treatmentScore["Metadata_Cell"] == ko]["InductionScore"]
        
        # Retrieve the scores corresponding to RKO CRBN OE
        oe_scores = treatmentScore[treatmentScore["Metadata_Cell"] == oe]["InductionScore"]
        
        # Quantify the difference between the induction scores calculated for each cell line
        u, p = mannwhitneyu(oe_scores, ko_scores, alternative = "greater")
        
        # Correct the u as a proprtion of the maximum possible u
        corrected_u = u / (len(ko_scores) * len(oe_scores))
        
        # Update the dictionary
        output_dict["Metadata_Treatment"].append(treatment)
        output_dict["Corrected_u"].append(corrected_u)
        output_dict["p-value"].append(p)
        
    # Convert the output_dict into a dataframe
    output_df = pd.DataFrame.from_dict(
        data = output_dict,
        orient = "columns"
    )
    return(output_df)

# Function for calculating the similarity in tau values
# AFTER de-selection of potentially redundant features
def comparePairwise_tauValues(
    selectedFeatures,
    tauValues
):
    """
    Calculate the similarity in tau values of treatment-centric features.

    :param selectedFeatures: Dataframe with the features selected using different strategies.
    :type selectedFeatures: pd.DataFrame
    :param tauValues: Dataframe with the tau values calculated during treatment-centric feature
        selection. The tau values loosely reflect the strength and type of correlation (positive or
        negative) between the change in a given feature and the expression of CRBN (or another protein)
        across the isogenic cell lines used.
    :type tauValues: pd.DataFrame
    :returns df: DataFrame with the cosine similarity scores calculated for each pair of compounds
        compared.
    :rtype df: pd.DataFrame
    """
    # Retrieve the features that would be selected for the controls
    # AFTER de-selection of redundant features
    treatment2features_dict = dict()
    for i in selectedFeatures.index:
        treatment = selectedFeatures.loc[i, "Treatment"]
        features = selectedFeatures.loc[i, "Features_Selected"].split(" ")
        treatment2features_dict[treatment] = features
    
    # Create a dataframe for storing the cosine similarity results
    df = pd.DataFrame(
        columns = list(treatment2features_dict.keys()),
        index = list(treatment2features_dict.keys())
    )
    
    # List for tracking which treatments have already had the overlap calculations done
    calculatedAlready = []
    
    # Calculate the overlap in treatment-centric features between each pair of treatments
    for treatment_a, features_a in treatment2features_dict.items():
        calculatedAlready.append(treatment_a)
        df.loc[treatment_a, treatment_a] = 1.
        for treatment_b, features_b in treatment2features_dict.items():
            if treatment_b not in calculatedAlready:
                featuresToCompare = list(set(features_a + features_b))
                A = tauValues[tauValues["Metadata_Treatment"] == treatment_a][featuresToCompare].to_numpy()
                B = tauValues[tauValues["Metadata_Treatment"] == treatment_b][featuresToCompare].to_numpy().T
                cosine = np.dot(A, B)/(norm(A) * norm(B))
                df.loc[treatment_a, treatment_b] = cosine[0][0]
                df.loc[treatment_b, treatment_a] = cosine[0][0]
    
    return(df)

# Function that calculates the pairwise cosine similarities between tau values of compounds
# if they are predicted to have CRBN-dependent bioactivity
def calculate_cosineSimilarity(
    selectedFeatures,
    tauValues,
    compoundData,
    output_path,
    correctedU_threshold = 0.7,
    featureNumber_threshold = 5
):
    """
    Calculate the pairwise cosine similarity scores between the tau values of the treatment-centric
    features of compounds that are predicted to have bioactivity that is dependent on the effector
    protein of interest like CRBN.
    
    :param selectedFeatures: Dataframe with the features selected using different strategies.
    :type selectedFeatures: pd.DataFrame
    :param tauValues: Dataframe with the tau values calculated during treatment-centric feature
        selection. The tau values loosely reflect the strength and type of correlation (positive or
        negative) between the change in a given feature and the expression of CRBN (or another protein)
        across the isogenic cell lines used.
    :type tauValues: pd.DataFrame
    :param compoundData: Dataframe with calculations/data for each compound (Corrected U scores and
        number of treatment-centric features required)
    :type compoundData: pd.DataFrame
    :param output_path: Path to export the dataframe of cosine similarity values to.
    :type output_path: str
    :param correctedU_threshold: Minimum corrected U score for a compound to be predicted as having
        bioactivity dependent on an effector protein of interest. Defaults to 0.7.
    :type correctedU_threshold: float (can be from 0 to 1)
    :param featureNumber_threshold: Minimum number of treatment-centric features for a compound to have
        a reliable corrected U score. In theory, a higher number should reflect higher reliability.
        Defaults to 5.
    :type featureNumber_threshold: int
    :returns cosine_df: Dataframe of cosine similarity values calculated.
    :rtype cosine_df: pd.DataFrame
    """
    ##################################################################################################
    # Prep the datasets for calculating the cosine similarities between the treatment-centric features
    # of compounds predicted to have CRBN-dependent bioactivity
    #################################################################################################
    # Retrieve the compounds that pass the criteria for CRBN-dependent bioactivity
    predictedCompounds = compoundData[
        (compoundData["Corrected_u"] >= correctedU_threshold) &
        (compoundData["Number_Features"] >= featureNumber_threshold)
    ]["Compound"].tolist()
    
    # Trim the selectedFeatures and tauValues to just the predictedCompounds
    selectedFeatures = selectedFeatures[selectedFeatures["Treatment"].isin(predictedCompounds)].reset_index(drop = True)
    tauValues = tauValues[tauValues["Metadata_Treatment"].isin(predictedCompounds)].reset_index(drop = True)
    
    #######################################################
    # Execute and export the cosine similarity calculations
    #######################################################
    # Calculate the pairwise cosine similarity between the tau values of the treatment-centric features
    # of the pair of compounds being compared
    cosine_df = comparePairwise_tauValues(
        selectedFeatures,
        tauValues
    )
    
    # Export the cosine similarities calculated 
    cosine_df.to_csv(output_path, index = True)
    print(f"EXPORTED: Pairwise cosine similarities between tau values of treatment-centric features has been calculated and exported at:/n{output_path}")
    return(cosine_df)