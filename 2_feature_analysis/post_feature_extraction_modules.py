#############################################
# !! WORKING VERSION !!
# Author(s): Amanda Ng R.H.
# Created on: 29 Nov 2022
# Last updated on: 21 Jul 2023 (edited the modify_featureExtractionCSV to make it compatible for running in bash)
# Documentation status:
# Good --- up till line 1776
# See https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html
# for more information on the Sphinx documentation style.
#
# General functions for use in the post-feature extraction pipleine of the Cell Painting Assay.
# Note: Not all functions have been fully tested.
#############################################

############
# Package(s)
############
# Warnings
import warnings

# File management
import os

# Dataframe management
import pandas as pd

# Plotting
import seaborn as sns
from matplotlib import pylab as plt

# Math
import statistics
import numpy as np
from scipy import stats

# Isolation Forest
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest

# Viewing images
import imageio
import cv2

# For interfacing with command line
import argparse

######################
# Internal function(s)
######################
def _vprint(message, verbose = True):
    """
    Internal function that prints a message depending on verbose option set by the user.
    
    :param message: Message the user would like to print out.
    :type message: str
    :param verbose: Option for whether the message should be printed out or not,
        defaults to True.
    :type verbose: bool
    :return:
    """
    if verbose:
        print(message)
    return

##############################################################################
# Function for adding a unique identifier for each site within the cell
# line-specific profile known as "UpdatedImageNumber" and merge the Image.csv
# files across all batches as {cell}__Image.csv
#############################################################################
def modify_featureExtractionCSV(
    cell,
    module_7b_output_dir,
    output_dir
):
    """
    Function that adds "UpdatedImageNumber" column to individual item CSV files from
    the feature extraction pipeline. The "UpdatedImageNumber" is a unique identifier for
    sets of images (one for each channel) associated to a unique site in a unique plate
    for a cell line.
    
    It also merges the Image.csv files from each batch into a single {cell}__Image.csv.
    
    :param cell: Name of the directory where all the cell line data is stored. It is also
        used for naming the merged Image.csv file. (e.g. "rko_wt")
    :type cell: str
    :param parent_dir: Path to the directory containing all the data associated to the
        Cell Painting Assay run.
    :type parent_dir: str
    :param module_7b_output_dir: Path to the directory for outputing the merged Image.csv
        ({cell}_Image.csv) to.
    :type module_7b_output_dir: str
    :return image_merged_path: Path to the {cell}_Image.csv.
    :rtype image_merged_path: str
    """    
    # List of feature extraction items except "Experiment"
    items = ["Cells", "Cytoplasm", "Nuclei", "Image"]
    
    # Assemble a dictionary for mapping the item to the paths
    # to their CSV files across batches
    item2paths_dict = dict()
    for item in items:
        item2paths_dict[item] = []
        for batch_dir in os.listdir(module_7b_output_dir):
            batch_dir = f"{module_7b_output_dir}/{batch_dir}"
            if os.path.isdir(batch_dir):
                path = f"{batch_dir}/{item}.csv"
                item2paths_dict[item].append(path)
    
    # Add the UpdateImageNumber column for all items
    # This UpdateImageNumber is a unique identifier for each image in the dataset
    for item, paths in item2paths_dict.items():
        batch = 0
        for path in paths:
            batch += 1
            
            # DEV: Check what the path looks like for the first batch
            if batch == 1:
                print("Example of path to the CSV files (for checking purposes):")
                print(path)
            
            # Load the dataframe and retrieve the last row of the dataframe
            df = pd.read_csv(path)
            last_row_number = len(df) - 1
            
            # Add the UpdatedImageNumber to the first CSV file
            # and get the maximum ImageNumber value
            if batch == 1:
                df["UpdatedImageNumber"] = df["ImageNumber"]
                max_updatedImageNumber = df.loc[last_row_number, "UpdatedImageNumber"]
            
            # Add the UpdatedImageNumber to the subsequent CSV files
            # and add the UpdatedImageNumber which is the ImageNumber + max_UpdatedImageNumber
            else:
                df["UpdatedImageNumber"] = df["ImageNumber"] + max_updatedImageNumber
                max_updatedImageNumber = df.loc[last_row_number, "UpdatedImageNumber"]
                
            # Export the dataframes with the UpdatedImageNumber added
            df.to_csv(path, index = False)
            
        print(f"{item} --- UpdatedImageNumber added")
        
        # Merge the image CSV files across all batches
        if item == "Image":
            df_list = []
            for path in item2paths_dict[item]:
                df = pd.read_csv(path)
                df_list.append(df)
            df = pd.concat(df_list)
            image_merged_path = f"{output_dir}/{cell}__{item}.csv"
            df.to_csv(image_merged_path, index = False)
            print(f"{item} --- merged CSV files")
            
    print(f"COMPLETED: UpdatedImageNumber has been added to the batch CSV files of the following items:\n{items}")
    print(f"COMPLETED: Merger of all Image.csv files across batches as\n{image_merged_path}")
    return(image_merged_path)

###########################################################################
# Function for trimming out flagged images from any of the object CSV files
###########################################################################
def remove_flaggedImages(cell, df, flagged_df):
    """
    Function for trimming out the flagged images from any object CSV file.
    
    :param cell: Name of cell line in use (e.g. "rko_wt").
    :type cell: str
    :param df: Dataframe of the object CSV file.
    :type df: pd.DataFrame
    :param flagged_df: Dataframe containing information on what images to trim out.
    :type flagged_df: pd.DataFrame
    :return:
    """
    initial_image_number = len(df)
    for ref in set(flagged_df["Reference_Column"].tolist()):
        valuesToRemove = flagged_df[flagged_df["Reference_Column"] == ref]["Value"].tolist()
        df = df[~df[ref].isin(valuesToRemove)]
    print(f"{cell} | Number of images removed: {initial_image_number - len(df)}")
    return(df)

#################################################################
# Functions for flagging problematic images with Isolation Forest
#################################################################
def get_imageQuality_features(image_df, verbose = True):
    """
    Function for retrieving a list of features that correspponds to measures of image quality.
    
    :param image_df: Dataframe containing general measurements of Image object.
    :type image_df: pd.DataFrame
    :param verbose: Option for whether the progress message should be printed out or not,
        defaults to True.
    :type verbose: bool, optional.
    :return features: List of features that correspond to image quality.
    :rtype features: list
    """
    all_features = image_df.columns.tolist()
    i = 0
    features = []
    for f in all_features:
        if "ImageQuality" in f:
            features.append(f)
            _vprint(f"{str(i)}\t{f}\tQuality", verbose)
            i += 1
        if "Image" not in f:
            _vprint(f"{str(i)}\t{f}\tMeasurement", verbose)
            i += 1
    features.sort()
    return(features)

def flag_lowQuality_images(image_df, featuresToComputeOn):
    """
    Function that flages low quality images using Isolation Forest.
    
    :param image_df: Dataframe containing general measurements of Image object.
    :type image_df: pd.DataFrame
    :param featuresToComputeOn: List of features pertaining to image quality.
    :type featuresToComputeOn: list
    :return principal_df: Dataframe containing 2D PCA values calculated from the image quality measurements.
    :rtype principal_df: pd.DataFrame
    """
    ### Identify outliers using IsolationForest ###
    image_df = image_df.dropna(
        subset = featuresToComputeOn
    )
    
    x = image_df.loc[:, featuresToComputeOn].values
    y = image_df.loc[:, ["UpdatedImageNumber"]].values

    x = StandardScaler().fit_transform(x)

    print("Train IsolationForest")
    clf = IsolationForest(
        n_jobs = 1,
        n_estimators = 500,
        max_samples = "auto",
        contamination = "auto",
        random_state = 1000
                         )
    clf.fit(x)
    y_pred_train = clf.predict(x)
    outliers = y_pred_train
    print("Finished training")
    
    image_df["IsolationForest"] = outliers
    
    ### Get the number of bad images and OK images ###
    ### and prep the dataframe for visualizing in 2D PCA ###

    vectorlength = np.linalg.norm(x, axis = 1)

    pca = PCA(n_components = 2)
    principalComponents = pca.fit_transform(x)

    principal_df = pd.DataFrame(
        data = principalComponents,
        columns = ["PC1", "PC2"]
                               )
    principal_df = pd.concat(
        [principal_df, image_df["UpdatedImageNumber"], image_df["IsolationForest"]],
        axis = 1
    )
    
    del image_df

    print(pca.explained_variance_ratio_)
    print(pca.explained_variance_)
            
    return(principal_df)

def stacked_barplot(cell2principal_dict, output_dir = ""):
    """
    Function that plots the percentage of bad and good images as a stacked barplot for each dataset.
    
    :param cell2principal_dict: Dictionary that maps cell names to the 2D PCA dataframe
        of image quality measurements.
    :type cell2principal_dict: dict
    :param output_dir: Path to the directory for outputs, defaults to "". Plots (of the image quality PCA and
        the stacked barplot of good and bad images) are not saved unless a non-zero length string is provided.
    :type output_dir: str, optional.
    :return:
    """
    barplot_dict = dict()
    for cell in cell2principal_dict:
        principal_df = cell2principal_dict[cell]
        bad = len(principal_df[principal_df["IsolationForest"] == -1])
        ok = len(principal_df[principal_df["IsolationForest"] == 1])
        percentage_bad = (bad/(bad + ok)) * 100
        percentage_ok = (ok/(bad + ok)) * 100
        barplot_dict[cell] = [percentage_bad, percentage_ok]
        
    barplot_df = pd.DataFrame.from_dict(
        data = barplot_dict,
        orient = "index",
        columns = [
            "Percentage_BadImages",
            "Percentage_OkImages"
        ]
    )
    
    ax = barplot_df.plot(
        kind = "bar",
        stacked = True,
        color = [
            "red",
            "darkgrey"
        ]
    )
    
    ax.legend(bbox_to_anchor = (1.0, 1.0))
    
    plt.ylim(0, 20)
    
    plt.xlabel("Cell")
    plt.ylabel("% bad/ok images assigned by Isolation Forest")
    plt.title("Overall quality of images according to Isolation Forest")
    
    if output_dir != "":
        output_path = f"{output_dir}/IsolationForest_ImageQualityStackedBarPlot.png"
        plt.savefig(output_path, dpi = 200, transparent = True, bbox_inches = "tight")
        print(f"EXPORTED: Stacked bar plot showing the proportion of bad and ok images has been exported at\n{output_path}")
    plt.show()
    
    return

#######################################################
# Functions for retrieving images for visual inspection
#######################################################
class ImageManipulation():
    """
    Shows colored image(s).
    """
    
    def colorize(im, color, clip_percentile = 0.1):
        """
        Helper function to create an RGB image from a single-channel image using a specific color.
        
        :param im: Array representing an image extracted either by imageio.imread or imageio.volread.
        :type im: np.array
        :param color: Three element tuple in RGB format following emission spectra e.g. (1, 0, 0) --> red.
        :type color: tuple
        :param clip_percentile: I think this is the threshold for what pixels color, defaults to 0.1.
        :type clip_percentile: float, optional. Values accepted between 0 to 100.
        :return im_scaled * color: Colored image as an array.
        :rtype im_scaled.color: np.array
        """
        # Check that we do just have a 2D image
        if im.ndim > 2 and im.shape[2] != 1:
            raise ValueError("This function expects a single-channel image!")

        # Rescale the image according to how we want to display it
        im_scaled = im.astype(np.float32) - np.percentile(im, clip_percentile)
        im_scaled = im_scaled / np.percentile(im_scaled, 100 - clip_percentile)
        im_scaled = np.clip(im_scaled, 0, 1)

        # Need to make sure we have a channels dimension for the multiplication to work
        im_scaled = np.atleast_3d(im_scaled)

        # Reshape the color (here, we assume channels last)
        color = np.asarray(color).reshape((1, 1, -1))
        return(im_scaled * color)

    def find_ImageByNumber(
        cell,
        df,
        updatedImageNumber,
        channel = "all"
    ):
        """
        Function to retrieve images from all five channels for viewing based on
        the "UpdatedImageNumber" and dataframe containing image metadata.

        :param df: Dataframe containing the features per image.
        :type df: pd.DataFrame
        :param updatedImageNumber: UpdatedImageNumber.
        :type updatedImageNumber: int
        :param channel: Option for which channel(s) to show images of, defaults to "all".
        :type channel: str, optional. Can be: "all", "Mito", "DNA", "ER", "AGP" or "BF".
        :return:
        """
        
        # Trim down the df to just the row pertaining to the updatedImageNumber
        df = df[df["UpdatedImageNumber"] == updatedImageNumber].reset_index(drop = True)
        
        # Get the path to the directory with the raw images
        raw_image_dir = df.loc[0, "PathName_OrigDNA"]
        
        # Get the plate
        plate = df.loc[0, "Metadata_Plate"]
        
        # Get the treatment
        treatment = df.loc[0, "Metadata_Treatment"]

        # Get the well info
        well = df.loc[0, "Metadata_Well"]

        # Get the site info
        site = str(df.loc[0, "Metadata_Site"]).zfill(2)
        
        # Define the channel names
        channel_names = ["Mito", "DNA", "ER", "AGP", "BF"]
        
        # Dictionary for mapping the channel names to the color tuples
        channel2color_dict = dict()
        channel2color_dict["Mito"] = (1, 0, 0) # red
        channel2color_dict["DNA"] = (0, 0, 1) # blue
        channel2color_dict["ER"] = (1, 1, 0) # yellow
        channel2color_dict["AGP"] = (0, 1, 0) # green
        channel2color_dict["BF"] = (1, 1, 1) # white 
        
        # Show all channels
        if channel == "all":

            # Dictionary mapping the channel name to the path
            channel2im_dict = dict()
            for channel_name in channel_names:
                if channel_name == "BF":
                    file = df.loc[0, f"FileName_OrigBrightfield"]
                else:
                    file = df.loc[0, f"FileName_Orig{channel_name}"]
                channel2im_dict[channel_name] = f"{raw_image_dir}/{file}"    

            # Make a plot with all the images side by side
            fig = plt.figure(figsize = (5*len(channel_names), 5))
            for i, channel_name in enumerate(channel_names):

                plt.subplot(1, 5, i + 1)

                color = channel2color_dict[channel_name]
                im = imageio.volread(channel2im_dict[channel_name])
                im_color = ImageManipulation.colorize(im[...], color)
                plt.imshow(im_color)

                plt.title(channel_name)
                plt.axis(False)

            fig.suptitle(
                f"{cell} {plate}: {treatment} (site {site} in well {well})",
                fontsize = 16,
                y = 1.05
            )
            plt.subplots_adjust(wspace = 0.05, hspace = 0.05)

            plt.show()
            
        # Show only one channel with a histogram (to check for saturation/intensity peaking)
        elif channel in channel_names:
            channel_name = channel
            
            # Get the path to the image
            if channel == "BF":
                file = df.loc[0, f"FileName_OrigBrightfield"]
            else:
                file = df.loc[0, f"FileName_Orig{channel_name}"]
            im_path = f"{raw_image_dir}/{file}"
            print(im_path)
            
            # Initialize the figure
            fig = plt.figure(figsize = (12, 6))
            
            for subplot in [1, 2]:
                
                plt.subplot(1, 2, subplot)
                
                im = imageio.volread(im_path)
                
                # Colorize the image and plot it in the first subplot
                if subplot == 1:
                    color = channel2color_dict[channel_name]
                    im_color = ImageManipulation.colorize(im[...], color)
                    plt.imshow(im_color)
                    plt.axis(False)
                
                # Get the histogram for the image and plot it in the second subplot
                # with a red dotted line marking the 1/3 point on the channel intensity scale
                if subplot == 2:
                    color = [1, 1, 1]
                    im_color = ImageManipulation.colorize(im[...], color)
                    histogram, bin_edges = np.histogram(im_color, bins = 256, range = (0, 1), density = True)
                    plt.plot(bin_edges[0: -1], histogram, color = "grey")
                    plt.xlabel("Raw channel intensity value")
                    plt.ylabel("Percentage of pixels")
                    plt.axvline(x = 1/3, ls = "--", color = "red")
            
            fig.suptitle(f"{cell} {plate}: {treatment} (site {site} in well {well})")

            plt.show()
        
        else:
            raise TypeError(f"channel input not accepted.")
        
        return
    
######################################################
# Function for retrieving images used for segmentation
# overlaid with segmentation masks
# for visual inspection of segmentation performance
######################################################
def visualize_segmentation(
    updatedImageNumber,
    image_df,
    parent_dir,
    cell,
    batches = 10,
    figsize = (7, 7)
):
    """
    Function for showing how the segmentation went for a given well and site.
    
    :param updatedImageNumber: UpdatedImageNumber value associated to the image the user would like to retrieve.
    :type updatedImageNumber: int
    :param image_df: Dataframe containing general measurements of Image object.
    :type image_df: pd.DataFrame
    :param parent_dir: Path to the directory containing all the data pertaining to the Cell Painting Assay run.
    :type parent_dir: str
    :param cell: Name of cell line in use (e.g. "rko_wt").
    :type cell: str
    :param batches: Maximum number of batches used during the feature extraction, defaults to 10.
    :type batches: int, optional.
    :param figsize: 2-element tuple for setting the size of the images shown, defaults to (7, 7).
    :type figsize: tuple
    :return:
    """
    # Trim down the image_df just to the relevant row
    image_df = image_df[image_df["UpdatedImageNumber"] == updatedImageNumber].reset_index(drop = True)
    
    # Define the paths to the segmentation input and output directories
    module_4b_dir = f"{parent_dir}/{cell}/output_dir/module_4b_output"
    module_7b_dir = f"{parent_dir}/{cell}/output_dir/module_7b_output"
    
    # Retrieve the dual channel image used for segmenting whole cells
    plate = image_df.loc[0, "Metadata_Plate"]
    well = image_df.loc[0, "Metadata_Well"]
    site = image_df.loc[0, "Metadata_Site"]
    filename = f"{plate}_{well}_s{site}__dna_agp.tiff"
    for batch_number in range(1, batches + 1):
        search_dir = f"{module_4b_dir}/{batch_number}_dna_agp"
        if filename in os.listdir(search_dir):
            filepath = f"{search_dir}/{filename}"
            break
    else:
        raise ValueError(f"Sub-directory containing the dual channel image cannot be found.\nDual channel image filename: {filename}")
    
    # Retrieve the file name of the nuclei outlines
    nuclei_filename = f"{plate}_{well}_s{site}__nuclei_outlines.png"
    
    # Retrieve the file name of the cell outlines
    # Note: These cell outlines are those which match a nuclei outline
    cell_filename = f"{plate}_{well}_s{site}__cell_outlines.png"
    
    # Set the figure size
    fig = plt.figure(figsize = figsize)
    
    # Open the dual channel image (which will be the background image)
    background_image = cv2.imread(filepath)
    
    # List of color tuples for the outlines
    colors = [
        (255, 255, 255), # white
        (255, 255, 255)  # yellow
    ]
    
    # Overlay each outline onto the background image
    for color, filename in zip(colors, [nuclei_filename, cell_filename]):
        
        # Retrieve the full path to the outline
        outline_filepath = f"{module_7b_dir}/{batch_number}/{filename}"
        
        # Open the outline image
        outline_image = cv2.imread(outline_filepath, cv2.IMREAD_UNCHANGED)
        
        # Convert the outline_image into a boolean array
        # where 0 = background and 255 = outline
        # with no in-between values
        idx = outline_image[:, :] == 255
        
        # Reassign the pixels in the dual channel image to white or yellow
        # according to whether the pixel is part of the outline or background
        # in the idx (boolean array)
        background_image[:, :, :][idx] = color
        
    # Show the dual channel image overlaid with the outlines
    plt.imshow(background_image)
    plt.title(f"{cell}: site {site} of well {well}")
    plt.axis(False)
    plt.show()
    
    return

#################################################################
# Class for prepping the morphological profile for each cell line
#################################################################
class PrepProfile():
    """
    Class used for prepping the morphological profiles of each cell line.
    """
    
    def __init__(self, show_max_memory_usage = True, verbose = True):
        """
        Constructor method
        """
        # Inform the user that the PrepProfile class has been initialized
        print(f"INITIALIZE: {self}")
        print(f"show_max_memory_usage set as {show_max_memory_usage}")
        print(f"verbose set as {verbose}")
        print("~"*10)
        
        # Define the option for showing the maximum RAM usage
        self.show_max_memory_usage = show_max_memory_usage
        
        # Define the verbosity
        self.verbose = verbose
        
        # Define the features to add from the Image CSV file to the assembled profile
        # since these details aren't in the object CSV files
        self.featuresFromImagecsv = [
            "Metadata_Plate", "Metadata_Run", "Metadata_Site",
            "Metadata_Treatment", "Metadata_Well",
            "Count_RelatedUnfilteredCells"
        ]
        
    def _vprint(self, message, verbose = True):
        """
        Internal function that prints a message depending on verbose option set by the user.

        :param message: Message the user would like to print out.
        :type message: str
        :param verbose: Option for whether the message should be printed out or not.
        :type verbose: bool
        :return:
        """
        if verbose:
            print(message)
        return
    
    def _mergeObjectCSV(self, cells_path, nuclei_path, cytoplasm_path):
        """
        Merges the individual object CSV files while keeping
        the features measured per object unique by appending suffixes:
        cytoplasm --> _cyto
        cells --> _cells
        nuclei --> (No suffix appended)
        
        :param cells_path: Path to the Cells_merged.csv.
        :type cells_path: str
        :param nuclei_path: Path to the Nuclei_merged.csv.
        :type nuclei_path: str
        :param cytoplasm_path: Path to the Cytoplasm_merged.csv.
        :type cytoplasm_path: str
        :return df: The dataframe with the Cells, Nuclei and Cytoplasm CSV file data merged.
        :rtype df: pd.DataFrame
        """
        self._vprint("Merging object CSV files...", verbose = self.verbose)
        
        #########################################
        # Merge the cytoplasm and cells CSV files
        #########################################
        # Load the cytoplasm and cells CSV files
        cytoplasm_df = pd.read_csv(cytoplasm_path)
        cells_df = pd.read_csv(cells_path)
        
        # Drop columns that are not relevant to the merger
        cytoplasm_df = cytoplasm_df.drop(columns = ["ImageNumber"]).dropna()
        cells_df = cells_df.drop(columns = ["ImageNumber"]).dropna()
        
        # Take note of the lenth of the cytoplasm and cells CSV files
        cytoplasm_length = len(cytoplasm_df)
        cells_length = len(cells_df)
        
        # Merge the cytoplasm and cells CSV files
        df = pd.merge(
            cytoplasm_df, cells_df,
            how = "inner",
            left_on = ["UpdatedImageNumber", "Parent_Cells"],
            right_on = ["UpdatedImageNumber", "ObjectNumber"],
            suffixes = ("_cyto", "_cells")
        )

        # Free up RAM
        del cytoplasm_df
        del cells_df
        
        ############################
        # Add in the nuclei CSV file
        ############################
        # Load the nuclei CSV file
        nuclei_df = pd.read_csv(nuclei_path)
        
        # Drop columns that are not relevant to the merger
        nuclei_df = nuclei_df.drop(columns = ["ImageNumber"]).dropna()
        
        # Take note of the length of the nuclei CSV file
        nuclei_length = len(nuclei_df)

        # Merge the df and nuclei_df
        df = pd.merge(
            df, nuclei_df,
            how = "inner",
            left_on = ["UpdatedImageNumber", "Parent_Nuclei"],
            right_on = ["UpdatedImageNumber", "ObjectNumber"],
        )
        

        self._vprint(
            f"""length of cytoplasm_df: {cytoplasm_length}
length of cells_df: {cells_length}
length of nuclei_df: {nuclei_length}
length of merged object dataframe: {len(df)}
COMPLETED: Object CSV files have been merged into a dataframe\n""", verbose = self.verbose
        )
        
        # Return the dataframe from the merger of the object datasets
        return(df)
    
    def _removeFlagged(self, flagged_path, df):
        """
        Function that removes the flagged images/treatments from the object CSV file during profile assembly.
        
        :param flagged_path: Path to the <cell line>__flagged.csv prepared in 1_FlagProblems.
        :type flagged_path: str
        :param df: Object CSV as a dataframe.
        :type df: pd.DataFrame
        :return df: The object CSV data without the objects associated to the flagged images/treatments.
        :rtype df: pd.DataFrame 
        """
        self._vprint("Filtering out images/treatments flagged as problematic for downstream analysis...",
                     verbose = self.verbose)
        
        # Load the <cell line>__flagged.csv file
        flagged_df = pd.read_csv(flagged_path)
        
        # Remove the images/treatments flagged
        for reference_col in set(flagged_df["Reference_Column"].tolist()):
            temp_df = flagged_df[flagged_df["Reference_Column"] == reference_col]
            value_list = temp_df["Value"].tolist()
            df = df[~df[reference_col].isin(value_list)]
            
        self._vprint(f"length of merged object dataframe: {len(df)}", verbose = self.verbose)
        self._vprint("COMPLETED: Flagged images/objects have been removed.\n", verbose = self.verbose)
        
        # Return the dataframe after filtering out
        # the problematic images and/or treatments
        return(df)
    
    def _aggregateToSite(self, df, groupby_cols = ["UpdatedImageNumber"]):
        """
        Function that aggregates the object level data to site level data.
        
        :param df: Dataframe containing information from all object CSV files.
        :type df: pd.DataFrame
        :return df: Dataframe containing information from all object CSV file aggregated on a per image
            i.e. site basis.
        :type df: pd.DataFrame
        """
        self._vprint("Aggregating object level data to site level data...", verbose = self.verbose)
        
        # Remove irrelevant metadata columns
        prefixToRemove = [
            "ImageNumber", "ObjectNumber",
            "Children", "Parent"
        ]
        metadataToRemove = []
        for col in df.columns:
            for prefix in prefixToRemove:
                if col.startswith(prefix):
                    metadataToRemove.append(col)
        df = df.drop(columns = metadataToRemove)
        
        # Save the data from the image CSV file that are already site level measurements
        # i.e. no aggregation is required
        # These data will be added back to the profile later
        reserve_df = df[groupby_cols + self.featuresFromImagecsv].drop_duplicates()
        
        # Remove the data from the image CSV file from the df prior to aggregation
        df = df.drop(columns = self.featuresFromImagecsv)
        
        # Aggregate measurements to the median per site
        df = df.groupby(groupby_cols).agg("median").reset_index(drop = False)
        
        # Add back the details that were separated into the reserve_df
        df = pd.merge(
            df, reserve_df,
            on = "UpdatedImageNumber",
            how = "inner"
        )
        
        self._vprint(f"length of aggregated dataset: {len(df)}", verbose = self.verbose)
        self._vprint("COMPLETED: Object level data has been aggregated to site level data for all cell lines.\n",
                     verbose = self.verbose)
        
        # Return the dataframe that contains the site level data
        return(df)
    
    def _addFromImageCSV(self, image_merged_path, df, feature):
        """
        Function that maps information on a given feature from the Image_merged.csv to the
        site level profile assembled.
        
        :param image_merged_path: Path to the Image_merged.csv.
        :type image_merged_path: str
        :param df: Dataframe containing the site level information (i.e. the profile).
        :type df: pd.DataFrame
        :param feature: Name of the feature to map over to the profile from the Image_merged.csv.
        :type feature: str
        :return df: Profile containing the feature information from the Image_merged.csv as a new column.
        :rtype df: pd.DataFrame
        """
        self._vprint(f"Adding {feature} from Image_merged.csv to the site level data...",
                     verbose = self.verbose)
        
        # Load the Image_merged.csv
        image_df = pd.read_csv(image_merged_path)
        
        # Extract a dictionary for mapping the
        # UpdatedImageNumber to the feature measured in Image_merged.csv
        imageNum2value_dict = dict()
        imageNum_list = image_df["UpdatedImageNumber"].tolist()
        value_list = image_df[feature].tolist()
        for imageNum, value in zip(imageNum_list, value_list):
            imageNum2value_dict[imageNum] = value
            
        # Add the number of valid cells
        # to the df containing the rest of the site-level measurements
        df[feature] = df["UpdatedImageNumber"].map(imageNum2value_dict)
        
        self._vprint(f"COMPLETED: {feature} added to the site level data.\n", verbose = self.verbose)
        
        # Return the dataframe that contains the site level data
        # with an additional feature from the Image_merged.csv
        return(df)

    def _removeBySampleSize(self, df):
        """
        Function for removing sites/images with <10 valid cells and
        treatments with <3 valid wells.
        
        :param df: Profile.
        :type df: pd.DataFrame
        :return df: Profile without the sites/images/treatments with a low sample size.
        :rtype df: pd.DataFrame
        """
        # Define the minimum number of valid objects and wells accepted
        min_valid_objects = 10
        min_valid_wells = 3
        self._vprint(f"""Removing sites with less than {min_valid_objects} objects and treatments with less than {min_valid_wells} wells...""",
                     verbose = self.verbose)
        
        # Remove site level data if there are fewer than
        # the min_valid_objects
        df = df[df["Count_RelatedUnfilteredCells"] >= min_valid_objects]
        self._vprint(f"Number of sites with >= {min_valid_objects} objects: {len(df)}",
                     verbose = self.verbose)
        
        # Remove treatments that have fewer than
        # the min_valid_wells
        treatmentsToRemove = []
        treatment_set = set(df["Metadata_Treatment"].tolist())
        for treatment in treatment_set:
            if len(df[df["Metadata_Treatment"] == treatment]) < 3:
                treatmentsToRemove.append(treatment)
        df = df[~df["Metadata_Treatment"].isin(treatmentsToRemove)]
        valid_treatment_set = set(df["Metadata_Treatment"].tolist())
        self._vprint(f"Initial number of treatments: {len(treatment_set)}",
                     verbose = self.verbose)
        self._vprint(f"Number of treatments with >= {min_valid_wells} wells: {len(valid_treatment_set)}",
                     verbose = self.verbose)
        
        self._vprint("COMPLETED: Sites and treatments with low sample size have been removed.\n",
                     verbose = self.verbose)
        
        # Return the site level dataframe containing
        # sites and treatments with a minimum sample size
        return(df)
               
    def profile(
        self,
        cell,
        module_7b_output_dir,
        image_merged_path,
        flagProblems_output_dir,
        profile_output_dir
    ):
        """
        Function that uses the above defined internal functions to assemble the site level
        morphological profiles for a given cell line. In brief:
        For a given cell line,
            For each batch,*
            - merge the individual object CSV files into one object level dataset
            - remove any objects associated to the problematic images/treatments
              flagged in 1_FlagProblems
            - aggregate the object level dataset into site level data
            - remove any sites with fewer than 10 valid cells and/or treatments with
              less than 3 valid wells
        - Merge the profile per batch into a single profile
        - finally export the final site level dataset i.e. the profile.
        
        *A batch-wise approach is used here to reduce RAM consumption during profile assembly,
        since merger after aggregation of object -> site level reduces RAM consumption by 300x.
        
        Note:
        To keep the job within shortq usage, the maximum amount of RAM is 8Gb, so
            1. the merged object CSV per batch cannot exceed 8 Gb
            2. the profile from the merger across batches cannot exceed 8 Gb.
        Assuming that one plate produces ~8 Gb of data across objects
        (i.e. cells, nuclei and cytoplasm), up to N plates can be handled with N batches,
        where N is a number from 1 to 300. Of course, you should NOT max out on the RAM usage,
        so always leave some RAM unused.
        
        :param cell: Name of the cell line (e.g. "rko_wt").
        :type cell: str
        :param module_7b_output_dir: Path to the directory where the object CSV files
            are located for each cell line (in separate batches).
        :type module_7b_output_dir: str
        :param image_merged_path: Path to the {cell}__Image.csv files are located.
        :type image_merged_path: str
        :param flagProblems_output_dir: Path to the directory with the {cell}__flagged.csv.
        :type flagProblems_output_dir: str
        :param profile_output_dir: Path to the directory to save the assembled profile for
            the cell line to.
        :type profile_output_dir: str
        :return profile_path: Path to the profile assembled for the given cell line.
        :rtype profile_path: str
        """
        print("profile outputs:")
        
        df_list = []
        i = 0
        for batch_dir in os.listdir(module_7b_output_dir):
            batch_dir = f"{module_7b_output_dir}/{batch_dir}"
            i += 1
            
            # Merge the object CSV files
            cells_path = f"{batch_dir}/Cells.csv"
            nuclei_path = f"{batch_dir}/Nuclei.csv"
            cytoplasm_path = f"{batch_dir}/Cytoplasm.csv"
            df = self._mergeObjectCSV(cells_path, nuclei_path, cytoplasm_path)
            
            # Add the metadata columns and the cell count to the df
            for feature in self.featuresFromImagecsv:
                df = self._addFromImageCSV(image_merged_path, df, feature)
                
            # Report on the maximum RAM usage for one batch
            if self.show_max_memory_usage == True and i == 1:
                print("ESTIMATED MAXIMUM MEMORY USAGE PER BATCH")
                print(f"Number of batches: {len(os.listdir(module_7b_output_dir))}")
                print("""df.info(memory_usage = "deep")""")
                print(df.info(memory_usage = "deep"))
                print("")

            # Remove any images and/or treatments
            # that have been flagged as problematic in 1FlagProblems
            flagged_path = f"{flagProblems_output_dir}/{cell}__flagged.csv"
            if not os.path.isfile(flagged_path):
                raise ValueError(f"""The CSV file containing the images and/or treatments flagged as problematic CANNOT BE FOUND. Did you name the file/folder differently? Or have you made a typo?
If so, please change it to match the following:
{flagged_path}""")
            df = self._removeFlagged(flagged_path, df)

            # Aggregate the object data to site data
            df = self._aggregateToSite(df)
            
            # Append the df assembled for the batch to the df_list
            df_list.append(df)
        
        # Concatenate the df assembled for each batch together
        df = pd.concat(df_list, ignore_index = True)
        
        # Remove sites with too few valid cells
        # and treatments with too few valid wells
        # Note
        df = self._removeBySampleSize(df)
        
        # Inform the user that the profile has been assembled for the given cell line
        print(f"COMPLETED: Profile for {cell} assembled.")
        print(f"Number of columns: {len(df.columns)}")
        print(f"Number of rows i.e. sites retained: {len(df)}")
        number_of_treatments = len(set(df["Metadata_Treatment"].tolist()))
        print(f"Number of treatments retained: {number_of_treatments}")
        number_of_features = len(df.drop(columns = self.featuresFromImagecsv + ["UpdatedImageNumber"]).columns.tolist()) # Note: I am actually not sure if this only counts the non-morphological features or if there is still metadata in the profile
        print(f"Number of morphological features retained:  {number_of_features}")
        print(f"Note: There is 1 non-morphological feature i.e. Count_RelatedUnfilteredCells.\n")
        
        # Export the final df
        # which is essentially the site level profile
        # with sites and treatments with sufficient valid cells and wells
        # and are not associated to problematic images/treatments
        profile_path = f"{profile_output_dir}/{cell}__profile.csv"
        df.to_csv(profile_path, index = False)
        print(f"""EXPORTED: Site level profile for {cell} has been exported to
{profile_path}""")
        
        # Display the top five rows of data in two specific columns
        print("Preview of profile (the top five rows of three columns):")
        display(df.head()[["Metadata_Treatment", "AreaShape_FormFactor", "Count_RelatedUnfilteredCells"]])
        
        return(profile_path)

#####################################################
# Class and functions pertaining to feature selection
#####################################################
# Note: Most functions in this class have been tested. I haven't included failsafes to catch ALL errors though,
# so there might be some errors that go undetected. I also have not described the error messages from this class
# fully in Sphinx documentation style, which I know can be annoying so I apologize in advance. For my use cases,
# the current failsafes catch the mistakes I have a tendency of making. Hopefully, these failsafes will also be
# sufficient for other users (if any).
class SelectFeatures():
    """
    Class that carries out both global and treatment-centric feature selection using the
    concatenated profile.
    """
    def __init__(
        self,
        concatProfile_path,
        prefix = "",
        corr_threshold = 0.8,
        verbose = True,
        output_dir = "default",
        param_global = dict({
            "skip": False,
            "get_consensus": False,
            "global_cell": "c662_rko_wt",
            "vote_threshold": 0.5
        }),
        param_treatment = dict({
            "skip": False,
            "keep_byProduct": True,
            "description": "",
            "ordered_cells": ["c1141_rko_ko", "c662_rko_wt", "c1327_rko_oe"],
            "alpha": 0.05,
            "kendall_alternative": "two-sided",
            "vote_threshold": 0.5
        })
    ):
        r"""
        Constructor method that also runs the global and treatment-centric feature selection
        using a combination of methods defined within the SelectFeatures class.
        
        :param concatProfile_path: Path to the concatenated profile i.e. concatenated__profile.csv.
        :type concatProfile_path: str
        :param prefix: Prefix to add to the names of output files, defaults to "".
        :type prefix: str, optional.
        :param corr_threshold: Threshold for what counts as a true corrleation between pairs of features. If
              the correlation of two features > corr_threshold, one of the features in the pair is
              discarded by _select_deccorelateByMAD. Defaults to 0.8.
        :type corr_threshold: float, optional. Can be from 0 to 1.
        :param verbose: Option for printing out all output messages or not, defaults to True.
        :type verbose: bool, optional.
        :param output_dir: Path to the directory for all outputs or "default". If "default", the output directory
            is set as the directory the concatenated profile is in. Defaults to "default".
        :type output_dir: str, optional.
        :param param_global: Dictionary of parameters to set for the global feature selection.
            Parameters (as keys in dictionary):
             - skip (bool)
                Option for skipping the global feature selection. If True, global feature selection is NOT
                executed, defaults to False.
            - get_consensus (bool or "Both")
                Option for global feature selection. Global feature selection involves
                the de-selection of features which are redundant using _select_deocrrelateByMAD, which
                MUST be done on a per cell line basis (since the baseline morphology is different
                across cell lines). Defaults to True.
                If get_consensus = True, the global feature selection will taks a
                consensus of features to select across cell lines.
                If False, global feature selection will only use the features from one cell line
                (which must be defined the user using the keyword argument global_cell
                e.g. global_cell = "rko_wt").
                If "Both", global feature selection will be done twice --- once with all cell lines
                (like in get_consensus = True) and once with a single cell line
                (like in get_consensus = False). "global_cell" must, thus, be defined.
            - vote_threshold (float, required if get_consensus = True)
                Threshold for what counts as a consensus vote for all cell lines in the the global feature
                selection (if get_consensus = True). Can be from 0 to 1.
            - global_cell (str, required if get_consensus = False)
                The name of the cell line that will be used for retrieving the cell line-specific profile.
                This profile will be used for global feature selection.
        :type param_global: dict, optional.
        :param param_treatment: Dictionary for parameters to set for the treatment-centric feature selection.
            Parameters (as keys in dictionary):
            - skip (bool)
                Option for skipping the treatment-centric feature selection. If True, treatment-centric
                feature selection is NOT excecuted.
            - keey_byProduct (bool)
                Option for saving the tau and p values calculated using the Kendall ranked correlation test
                (done by _calculateKendall). If true, the tau and p values are saved. Else, the values are
                NOT saved.
            - description (str)
                A string describing the purpose of the treatment-centric feature selection,
                defaults to "". This parameter is more for the user or anyone reading the code to
                know what the treatment-centric feature selection was carried out for. It has no impact
                on the outputs in this function.
            - ordered_cells (list)
                List of cell line names ORDERED by genetic perturbation,
                defaults to ["rko_ko", "rko_wt", "rko_oe"] where "rko_ko" has the lowest expression
                of the protein while "rko_oe" has the highest. This parameter is used by the
                treatment-centric feature selection ONLY.
            - alpha (float)
                Threshold used in treatment-centric feature selection. Features are selected if
                it has a significant correlation with cell line context, which is calculated as
                p-value by the Kendall rank correlation test (_calculateKendall). If the
                p-value < alpha after Bonferonni correction for multiple hypothesis testing,
                the correlation between the feature and cell line context (for the given treatment) is
                considered significant. Can be from 0 to 1.
            - kendall_alternative (string)
                Test type for Kendall's Tau-b correlation coefficient (see _calculateKendall). Can be "two-sided",
                "less" or "greater".
            - vote_threshold (float)
                Serves the same purpose as in global feature selection, but done within a treatment-centric
                context i.e. only profile data across all cell lines are used for getting a consensus vote
                for what features to select across all cell lines. Defaults to 0.5. Can be from 0 to 1.
                (This implies that get_consensus is always True for treatment-centric feature selection.)
        :type param_treatment: dict, optional.
        """
        print(f"INITIALIZE: {self}")
        
        # Define the output directory
        if output_dir == "default":
            output_dir = "/".join(concatProfile_path.split("/")[:-1])
        else:
            output_dir = output_dir
        
        # Define columns NOT to use for feature selection
        # These are features which are NOT measurements of images
        # and are categorical variables
        self.doNotUseCols = [
            "UpdatedImageNumber",
            "Metadata_Plate", "Metadata_Run", "Metadata_Site",
            "Metadata_Treatment", "Metadata_Well", "Metadata_Cell"
        ]
        
        # Make a dictionary for the compulsory CSV output(s)/input(s)
        self.key2path_dict = dict()
        key_list = [
            "baseline_output",
            "featuresSelected_output",
            "treatment2absent_output"
        ]
        
        # Add non-compulsory outputs depending on the user's input
        if "keep_byProduct" in param_treatment and param_treatment["keep_byProduct"] == True:
            key_list += ["tau_output", "p_output"]
        for key in key_list:
            self.key2path_dict[key] = f"{output_dir}/{prefix}{key}.csv"
        
        # Print output path(s) for easy referencing purposes
        msg = "OUTPUT_PATHS"
        divider = "-" * len(msg)
        print(f"{divider}\n{msg}\n{divider}")
        for key, path in self.key2path_dict.items():
            print(f"{key} = {path}")
        print("")
        print("=" * 10)     
        
        # Initialize dictionary for storing the features selected
        # depending on the type of interpretation to be made
        # from the data
        # key = treatment
        # value = [aim, featuresToKeep]
        self.output_dict = dict()
        
        # Retrieve the morphological profile for all relevant cell lines
        df = pd.read_csv(concatProfile_path)
        
        # Carry out the baseline feature selection (which removes the irrelevant features)
        # and calculates the robust Z score for each feature per image
        self._vprint("1: Baseline selection & robust Z calculation", verbose)
        df = self._select_baseline(df, verbose)
        output_path = self.key2path_dict["baseline_output"]
        df.to_csv(output_path, index = False)
        self._vprint(f"EXPORTED: Profile after baseline selection and robust Z calculation exported at\n{output_path}.", verbose)
        
        ### Global feature selection ###
        if param_global["skip"] == True:
            print("User has chosen to skip global feature selection.")
            if len(prefix) == 0:
                warnings.warn("User did not provide a 'prefix' for outputs. A prefix is recommended especially if you are skipping either global or treatment-centric feature selection.")
        elif param_global["skip"] == False:
            self._vprint("\n2: Starting global feature selection...")
            print("ADJUST PARAMETERES: User has provided inputs on the parameters for global feature selection.")

            # If get_consensus = False or "Both"
            # run the global feature selection only using data from one cell line
            # which MUST be provided by the user
            if param_global["get_consensus"] in [False, "Both"]:
                if "global_cell" in param_global:
                    global_cell = param_global["global_cell"]
                    acceptable_cells = set(df["Metadata_Cell"].tolist())
                    if global_cell in acceptable_cells:
                        print(f"MODE: Global feature selection carried out with profile data from {global_cell} ONLY.")
                    else:
                        raise ValueError(f"{global_cell} cannot be found in the profile given. Did you mean one of the following cells instead?\n{acceptable_cells}")
                else:
                    raise ValueError(f"""get_consensus set to {get_consensus}. "global_cell" must be provided in "global_strategy" dictionary.""")
                cell_df = df[df["Metadata_Cell"] == global_cell]
                featuresToKeep = self._select_decorrelateByMAD(
                    df = cell_df,
                    corr_threshold = corr_threshold,
                    verbose = verbose
                )
                del cell_df
                n = len(featuresToKeep)
                featuresToKeep = " ".join(featuresToKeep)
                self.output_dict[f"all_{global_cell}"] = [f"global_{global_cell}", featuresToKeep, n]

            # if get_consensus = True (default) or "Both"
            # run global feature selection using data from all cell lines
            if param_global["get_consensus"] in [True, "Both"]:
                print("MODE: Global feature selection carried out with profile data from ALL CELL LINES.")
                self._select_consensus(
                    df = df,
                    selection_strategy = "global",
                    treatment = "all",
                    vote_threshold = param_global["vote_threshold"],
                    corr_threshold = corr_threshold,
                    verbose = verbose
            )
        
        ### treatment-centric feature selection ###
        if param_treatment["skip"] == True:
            print("User has chosen to skip treatment-centric feature selection.")
            if len(prefix) == 0:
                warnings.warn("User did not provide a 'prefix' for outputs. A prefix is recommended especially if you are skipping either global or treatment-centric feature selection.")
        elif param_treatment["skip"] == False:
            self._vprint("\n2: Starting treatment-centric feature selection...")
            print("ADJUST PARAMETERES: User has provided inputs on the parameters for treatment-centric feature selection.")
            self._select_byTreatment(
                df = df[df["Metadata_Cell"].isin(param_treatment["ordered_cells"])],
                ordered_cells = param_treatment["ordered_cells"],
                keep_byProduct = param_treatment["keep_byProduct"],
                description = param_treatment["description"],
                alpha = param_treatment["alpha"],
                kendall_alternative = param_treatment["kendall_alternative"],
                vote_threshold = param_treatment["vote_threshold"],
                corr_threshold = corr_threshold,
                verbose = verbose
             )
            
        # Raise an error if the user has chosen to skip BOTH global and treatment-centric feature selection
        if param_global["skip"] == True and param_treatment["skip"]:
            raise ValueError("You are asking me to not carry out feature selection at all.")
            
        # Export the features selected for global and/or treatment-centric feature selection
        else:
            output_path = self.key2path_dict["featuresSelected_output"]
            output_df = pd.DataFrame.from_dict(
                data = self.output_dict,
                orient = "index",
                columns = ["Selection_Strategy", "Features_Selected", "Number_Features"]
            )
            output_df["Treatment"] = output_df.index.tolist()
            output_df.to_csv(output_path, index = False)
            self._vprint(f"EXPORTED: Features selected have been exported at\n{output_path}", verbose)

            # Inform the user on the status of the overall feature selection be it global and/or
            # treatment-centric feature selection
            self._vprint("Preview of the features selected (top five rows):", verbose)
            if verbose == True:
                display(output_df.head(n = 5))
            self._vprint("COMPLETED: Feature selection completed.", verbose)
    
    ###########
    # Verbosity
    ###########
    @staticmethod
    def _vprint(message, verbose = True):
        """
        Internal function that prints a message depending on verbose option set by the user.

        :param message: Message the user would like to print out.
        :type message: str
        :param verbose: Option for whether the message should be printed out or not.
        :type verbose: bool
        :return:
        """
        if verbose:
            print(message)
        return
    
    ###########################################
    # STANDARDIZE TO DMSO: Robust Z calculation
    ###########################################
    @staticmethod
    def _calculateRobustZ(df, featuresToComputeOn):
        """
        Log transforms features and calculates the robust Z score.
        
        :param df: Profile of a specific cell line.
        :type df: pd.DataFrame
        :param featuresToComputeOn: List of features to compute on.
        :type featuresToComputeOn: list
        :return df: The profile with the morphological features log-transformed and
            converted to Robust Z score.
        :rtype df: pd.DataFrame
        """
        for feature in featuresToComputeOn:
            
            # Log transform the data so that it approximates a
            # normal distribution
            df.loc[:, feature] = np.log(df[feature] + 1 - min(df[feature]))
            
            # Standardize the data to the DMSO control
            # (essentially calculating a robust z-score)
            dmso_df = df[df["Metadata_Treatment"] == "DMSO"]
            dmso_inputList = dmso_df[feature].tolist()
            mad_dmso = stats.median_abs_deviation(dmso_inputList)
            median_dmso = np.median(dmso_inputList)
            df.loc[:, feature] = (df[feature] - median_dmso) / mad_dmso
            
        # Return the dataframe comprising the robust z scores
        return(df)
    
    #####################################################
    # PER TREATMENT: Kendall Rank Correlation calculation
    #####################################################
    def _calculateKendall(
        self,
        df,
        featuresToComputeOn,
        ordered_cells = ["rko_ko", "rko_wt", "rko_oe"],
        alternative = "two-sided",
        verbose = True
    ):
        """
        Carry out the Kendall Rank Correlation test to calculate tau and p-value,
        which quantify the correlation between the feature and the cell line contexts used.
        
        :param df: Profile of a specific cell line.
        :type df: pd.DataFrame
        :param featuresToComputeOn: List of features to compute on.
        :type featuresToComputeOn: list
        :param ordered_cells: List of cell line names ORDERED by genetic perturbation,
            defaults to ["rko_ko", "rko_wt", "rko_oe"] where "rko_ko" has the lowest expression
            of the protein while "rko_oe" has the highest.
        :type ordred_cells: list
        :param alternative: Defines the alternative hypothesis. Defaults to "two-sided".
            Can be "two-sided", "less" or "greater".
        :type alternative: str
        :param verbose: Option for printing out all output messages or not, defaults to True.
        :type verbose: bool, optional.
        :return df_list[0]: Dataframe containing the tau values calculated for each feature per treatment.
        :rtype df_list[0]: pd.DataFrame
        :return df_list[1]: Dataframe containing the p values calculated for each feature per treatment.
        :rtype df_list[1]: pd.DataFrame
        """
        # Check if the list of cell lines matches the morphological profile provided
        if set(df["Metadata_Cell"].tolist()) != set(ordered_cells):
            raise ValueError(f'''_calculateKendall cannot run.
Names of cell lines in morphological profile does NOT match {ordered_cells}''')
            
        # Add the cell line information as a ordinal variable
        name2ordinal_dict = dict()
        for ordinal, name in enumerate(ordered_cells):
            name2ordinal_dict[name] = ordinal + 1
        df["Ordinal_Cell"] = df["Metadata_Cell"].map(name2ordinal_dict)
        
        # Calculate the Kendall tau for each feature
        # per treatment
        treatment2tau_dict = dict()
        treatment2p_dict = dict()
        dictionaries = [treatment2tau_dict, treatment2p_dict]
        for treatment in set(df["Metadata_Treatment"].tolist()):
            treatment_df = df[df["Metadata_Treatment"] == treatment]
            treatment2tau_dict[treatment] = []
            treatment2p_dict[treatment] = []
            for feature in featuresToComputeOn:
                x = treatment_df[feature].tolist()
                y = treatment_df["Ordinal_Cell"].tolist()
                tau, p = stats.kendalltau(x, y, alternative = alternative)
                for dictionary, value in zip(dictionaries, [tau, p]):
                    dictionary[treatment].append(value)
        
        # Convert the dictionaries into a dataframe
        df_list = []
        for dictionary in dictionaries:
            df = pd.DataFrame.from_dict(data = dictionary, orient = "index", columns = featuresToComputeOn)
            df["Metadata_Treatment"] = df.index.tolist()
            df = df.reset_index(drop = True)
            df_list.append(df)
        
        self._vprint("COMPLETED (1/2): Kendall tau calculated.", verbose)
        self._vprint("COMPLETED (2/2): Number of features with significant Kendall p-value calculated.", verbose)
        return(df_list[0], df_list[1])
    
    ###########################
    # Feature selection methods
    ###########################  
    def _select_baseline(self, df, verbose = True):
        """
        Carry out baseline feature selection by:
        1. Discarding irrelevant features that are constant across
           all treatments.
        2. Features in DMSO controls with an MAD of 0, which compromises
           the downstream calculation of the robust Z score.
        3. Converting the feature measurements to robust Z scores
           PER cell line context.
        4. Discarding features with "NaN" or "inf" values after
           robust Z score calculation.
        
        :param df: Profile with all the features prior to any feature selection.
        :type df: pd.DataFrame
        :param verbose: Option for printing out all output messages or not, defaults to True.
        :type verbose: bool, optional.
        :returns robustZ_df: Profile with features retained after baseline feature selection
            that have been converted into robust Z score.
        :rtype df: pd.DataFrame
        """
        self._vprint("select_baseline outputs:", verbose)
        
        # Retain features with
        # non-zero median absolute deviation (MAD) across all treatments
        self._vprint("Discarding features with zero MAD across all treatments and DMSO controls...", verbose)
        dmso_df = df[df["Metadata_Treatment"] == "DMSO"]
        featuresToSelectFrom = df.drop(columns = self.doNotUseCols).columns.tolist()
        featuresToKeep = []
        for feature in featuresToSelectFrom:
            inputList = df[feature].tolist()
            mad = stats.median_abs_deviation(inputList, scale = "normal")
            if mad != 0:
                
                # which also exhibit non-zero MAD for DMSO controls
                # (kicked out to allow for subsequent standardization)
                dmso_inputList = dmso_df[feature].tolist()
                dmso_mad = stats.median_abs_deviation(dmso_inputList, scale = "normal")
                if dmso_mad != 0:
                    featuresToKeep.append(feature)
                    
        metadata_cols = [col for col in self.doNotUseCols if col in df.columns.tolist()]
        df = df[featuresToKeep + metadata_cols]
        
        # Calculate the robust Z score for each feature PER cell line context
        self._vprint("Calculating robust z score for each feature...", verbose)
        robustZ_df_list = []
        for cell in set(df["Metadata_Cell"].tolist()):
            cell_df = df[df["Metadata_Cell"] == cell]
            robustZ_df = self._calculateRobustZ(cell_df, featuresToKeep)
            robustZ_df_list.append(robustZ_df)
        robustZ_df = pd.concat(robustZ_df_list)
        
        # Free up RAM
        del robustZ_df_list
        del df
        del cell_df
        
        # Drop features with infinity or nan
        robustZ_df = robustZ_df.replace([np.inf, -np.inf], np.nan).dropna(axis = 1)    
        
        # Return the dataframe after baseline feature selection
        # and robust Z calculation
        self._vprint("COMPLETED: Basline feature selection.", verbose)
        return(robustZ_df)
    
    def _select_decorrelateByMAD(self, df, corr_threshold = 0.8, verbose = True, descendingMAD = True):
        """
        Function that discards redundant features by:
        1. Calculating the variability of each feature.
        2. Ordering the features by variability in descending/ascending order.
        3. Select the first feature and iteratively compare the subsequent feature with
           the selected feature. If the correlation between the feature and prior
           selected feature(s) > 0.8, the feature is considered redundant and not included in
           the list of selected features.
        
        :param df: Profile of all cell lines after the features have been converted to
            robust Z score.
        :type df: pd.DataFrame
        :param corr_threshold: Threshold for what counts as a true corrleation between pairs
            of features. If the correlation of two features > corr_threshold, one of the
            features in the pair is discarded by _select_deccorelateByMAD.
            Defaults to 0.8. Can be from 0 to 1.
        :type corr_threshold: float, optional.
        :param verbose: Option for printing out all output messages or not, defaults to True.
        :type verbose: bool, optional.
        :param descendingMAD: If true, features are ordered in descending (high to low) variability. If false, the
            features are ordered in ascending (low to high) variability. Defaults to True.
        :type descendingMAD: bool, optional.
        :return featuresToKeep: List of features selected.
        :rtype featuresToKeep: list
        """
        self._vprint("_select_decorrelateByMAD outputs:", verbose)
        self._vprint("This function should ONLY be used on robust Z scores!", verbose)
        
        # Retrieve features to carry out the selection on
        featuresToSelectFrom = df.drop(columns = self.doNotUseCols).columns.tolist()
        
        # Record the overall variation per feature
        featureVariation_list = []
        for feature in featuresToSelectFrom:
            featureVariation = stats.median_abs_deviation(df[feature].tolist())
            featureVariation_list.append([feature, featureVariation])
        
        # Rank the feature by variability in z-score
        featureVariation_list.sort(key = lambda row: (row[1]), reverse = descendingMAD)
        
        # Discard features with high correlation
        self._vprint("Discarding features with high correlation...", verbose)
        i = 0
        featuresToKeep = []
        for feature, variation in featureVariation_list:
            i += 1            
            if i == 1:
                featuresToKeep.append(feature)
            else:
                for ref_feature in featuresToKeep:
                    x = df[ref_feature].tolist()
                    y = df[feature].tolist()
                    corr, _ = stats.pearsonr(x, y)
                    if corr > corr_threshold:
                        break
                else:
                    featuresToKeep.append(feature)

        # Return the featuresToKeep
        self._vprint("COMPLETED: Feature selection by decorrelatedByMAD.", verbose)
        return(featuresToKeep)

    def _select_consensus(
        self,
        df,
        selection_strategy = "global",
        treatment = "all",
        descendingMAD = True,
        vote_threshold = 0.5,
        corr_threshold = 0.8,
        verbose = True
    ):
        """
        Each cell line's morphological profile is dealt with independently.
        For each profile, redundant features are pruned by
        _select_decorrelateByMAD. Each feature selected within a cell line
        context is considered a vote to select the feature. A consensus of the
        votes is taken across cell line contexts based on the vote_threshold set.
        
        The final set of features selected is then recorded.
        
        :param df: Profile of all cell lines after the features have been converted to
            robust Z score.
        :type df: pd.DataFrame
        :param selection_strategy: Description of the selection strategy currently in use. It is
            also used to fill in the featuresSelected_output to distinguish between features
            that have been selected in the "global" or "treatment-centric" manner. Defaults to "global".
        :type selection_strategy: str, optional.
        :param treatment: Indicates what treatment is currently being looked at and will be recorded
            in the featuresSelected_output. Defaults to "all".
        :type treatment: str, optional.
        :param descendingMAD: If true, features are ordered in descending (high to low) variability. If false, the
            features are ordered in ascending (low to high) variability. Defaults to True.
        :type descendingMAD: bool, optional.
        :param vote_threshold: Threshold for what counts as a consensus vote for all cell lines,
            defaults to 0.5. Can be from 0 to 1.
        :type vote_threshold: float, optional.
        :param corr_threshold: Threshold for what counts as a true corrleation between pairs
            of features. If the correlation of two features > corr_threshold, one of the
            features in the pair is discarded by _select_deccorelateByMAD.
            Defaults to 0.8. Can be from 0 to 1.
        :type corr_threshold: float, optional.
        :param verbose: Option for printing out all output messages or not, defaults to True.
        :type verbose: bool, optional.
        :return:
        """
        self._vprint("Removing redundant features...", verbose)
        
        # Initialize a dictionary for storing whether the feature
        # is retained or not after _select_decorrelateByMAD
        feature2vote = dict()
        
        # Collect votes on features to keep across all cell lines
        # in the merged profile
        self._vprint("-", verbose)
        self._vprint("Collecting votes for features to select from:", verbose)
        cell_line_list = list(set(df["Metadata_Cell"].tolist()))
        for i, cell_line in enumerate(cell_line_list):
            self._vprint(f"{i + 1} out of {len(cell_line_list)} cell lines", verbose)
            cell_df = df[df["Metadata_Cell"] == cell_line]
            featuresToKeep = self._select_decorrelateByMAD(
                df,
                corr_threshold,
                verbose = False,
                descendingMAD = descendingMAD
            )
            for feature in featuresToKeep:
                if feature not in feature2vote:
                    feature2vote[feature] = 1
                else:
                    feature2vote[feature] += 1
            self._vprint("-", verbose)
                    
        # Retain features which have votes above the vote_threshold
        featuresToKeep = []
        for feature, vote in feature2vote.items():
            if vote >= vote_threshold * len(cell_line_list):
                featuresToKeep.append(feature)
        
        # Update the dictionary of features selected
        # depending on the selection strategy and/or treatment
        n = len(featuresToKeep)
        featuresToKeep = " ".join(featuresToKeep)
        value_list = [selection_strategy, featuresToKeep, n]
        self.output_dict[treatment] = value_list
        
        # Inform the user on the completion of the feature selection
        # and the number of features selected
        number_featuresToSelectFrom = len(df.drop(columns = self.doNotUseCols).columns.tolist())
        self._vprint(f"{n} out of {number_featuresToSelectFrom} features selected.", verbose)
        self._vprint("COMPLETED: Redundant features removed.", verbose)
        return

    def _select_byTreatment(
        self,
        df,
        ordered_cells = ["rko_ko", "rko_wt", "rko_oe"],
        keep_byProduct = True,
        description = "",
        alpha = 0.05,
        kendall_alternative = "two-sided",
        vote_threshold = 0.5,
        corr_threshold = 0.8,
        verbose = True
        ):
        """
        The morphological profile across ALL cell lines are considered here
        on a PER treatment basis. Features here are selected if they exhibit a
        correlation with the cell lines, which is treated as ordinal data e.g.:
            - rko_ko --> 0
            - rko_wt --> 1
            - rko_oe --> 2
            (see _calculateKendall for more details)
            
        These selected features PER treatment are then recorded.
        
        :param df: Profile of all cell lines with features converted to robust Z score.
        :type df: pd.DataFrame
        :param ordered_cells: List of cell line names ORDERED by genetic perturbation,
            defaults to ["rko_ko", "rko_wt", "rko_oe"] where "rko_ko" has the lowest expression
            of the protein while "rko_oe" has the highest.
        :type ordered_cells: list
        :param keep_byProduct: Option for saving the outputs from the Kendall rank correlation test
            (calculated by _calculateKendall), defaults to True.
        :type keep_byProduct: bool, optional.
        :param description: Description of what the treatment-centric feature selection is for.
            This parameter is mostly for the user to see what this function is being used for in
            his/her/their code.
        :type description: str, optional.
        :param alpha: Threshold used in treatment-centric feature selection. Features are selected if
            it has a significant correlation with cell line context, which is calculated as
            p-value by the Kendall rank correlation test (_calculateKendall). If the
            p-value < alpha after Bonferonni correction for multiple hypothesis testing,
            the correlation between the feature and cell line context (for the given treatment) is
            considered significant. Defaults to 0.05. Can be from 0 to 1.
        :type alpha: float, optional.
        :param kendall_alternative: Defines the alternative hypothesis. Defaults to "two-sided".
            Can be "two-sided", "less" or "greater".
        :type alternative: str
        :param vote_threshold: Threshold for what counts as a consensus vote for all cell lines,
            defaults to 0.5. Can be from 0 to 1.
        :type vote_threshold: float, optional.
        :param corr_threshold: Threshold for what counts as a true corrleation between pairs
            of features. If the correlation of two features > corr_threshold, one of the
            features in the pair is discarded by _select_deccorelateByMAD.
            Defaults to 0.8. Can be from 0 to 1.
        :type corr_threshold: float, optional.
        :param verbose: Option for printing out all output messages or not, defaults to True.
        :type verbose: bool, optional.
        """
        self._vprint("[PER TREATMENT] Selecting features which correlate with the genetic perturbation...", verbose)
        if len(description) > 0:
            self._vprint(description)
            
        # Check if there are any treatments that aren't present across all cell line profiles
        treatment2absent_dict = dict()
        for treatment in set(df["Metadata_Treatment"].tolist()):
            temp_df = df[df["Metadata_Treatment"] == treatment]
            for cell in ordered_cells:
                if cell not in set(temp_df["Metadata_Cell"].tolist()):
                    if treatment not in treatment2absent_dict:
                        treatment2absent_dict[treatment] = cell
                    else:
                        treatment2absent_dict[treatment] += f" {cell}"
        if len(treatment2absent_dict.keys()) > 0:
            absent_df = pd.DataFrame.from_dict(
                treatment2absent_dict,
                orient = "index",
                columns = ["NotRetainedInCells"]
            )
            absent_df["Treatment"] = absent_df.index.tolist()
            output_path = self.key2path_dict["treatment2absent_output"]
            absent_df.to_csv(output_path, index = False)
            print("Some treatments were not retained across all cell lines.")
            print(f"EXPORTED: Treatments and the cell line contexts they are not retained in have been exported at\n{output_path}")
            df = df[~df["Metadata_Treatment"].isin(list(treatment2absent_dict.keys()))]
            print("Treatments not retained across cell lines discarded from profile prior to treatment-centric feature selection.")
        
        # Retrieve features to carry out the selection on
        featuresToSelectFrom = df.drop(columns = self.doNotUseCols).columns.tolist()
                
        # Calculate the correlation between
        # the robust Z score of the feature
        # and the genetic perturbation in cell line
        # per treatment
        tau_df, p_df = self._calculateKendall(
            df = df,
            featuresToComputeOn = featuresToSelectFrom,
            ordered_cells = ordered_cells,
            alternative = kendall_alternative,
            verbose = verbose
        )
        
        # Export the tau_df and p_df
        # if the user wishes to keep by products from the pipeline
        if keep_byProduct == True:
            tau_df.to_csv(self.key2path_dict["tau_output"], index = False)
            p_df.to_csv(self.key2path_dict["p_output"], index = False)
            self._vprint(f"EXPORTED: Kendall tau and p-values calculated.", verbose)
        
        # Retain features with a p-value < alpha after Bonferroni correction
        # i.e. features with a significant correlation to the genetic perturbation
        corrected_alpha = alpha/len(featuresToSelectFrom)
        treatment_list = p_df["Metadata_Treatment"].tolist()
        for i, treatment in enumerate(treatment_list):
            featuresToKeep = []
            treatment_p_df = p_df[p_df["Metadata_Treatment"] == treatment]
            for feature in featuresToSelectFrom:
                if p_df.loc[i, feature] < corrected_alpha:
                    featuresToKeep.append(feature)
            
            # Prune redundant features per treatment from the featuresToKeep
            # if there are at least 2 features that correlate with the genetic perturbation
            # for the treatment
            if len(featuresToKeep) > 1:
                treatment_df = df[df["Metadata_Treatment"] == treatment]
                treatment_df = treatment_df[self.doNotUseCols + featuresToKeep]
                self._select_consensus(
                    df = treatment_df,
                    selection_strategy = "treatment-centric",
                    treatment = treatment,
                    vote_threshold = vote_threshold,
                    corr_threshold = corr_threshold,
                    descendingMAD = False,
                    verbose = False
                )
                
            # Update the dictionary for features selected directly
            # if there are less than 2 features that correlate with the genetic perturbation
            else:
                n = len(featuresToKeep)
                featuresToKeep = " ".join(featuresToKeep)
                value_list = ["treatment-centric", featuresToKeep, n]
                self.output_dict[treatment] = value_list
                
            # Update the user on the progress every 20 treatments
            if verbose == True:
                if (i + 1) % 20 == 0:
                    print(f"{i + 1} out of {len(treatment_list)}")
                    
        self._vprint("COMPLETED: [PER TREATMENT] Features that correlate with the genetic perturbation have been selected.", verbose)
        return
    
###############################################################
# Set up the interfacing of function(s) for use in command line
###############################################################
if __name__ == "__main__":

    # If the user types --help
    parser = argparse.ArgumentParser(description = "Functions in post_feature_extraction_modules.py: modules for pre-data processing in between feature extraction and profile interpretation. Note that only one function from this module has been set up for use in command line.")
    
    # If the user types --help for a specific function
    parser.add_argument(
        "--modify_featureExtractionCSV",
        nargs = 3,
        help = """
Usage: python3 post_feature_extraction_modules.py --modify_featureExtractionCSV ${cell} ${module_7b_output_dir} ${output_dir}
Returns: Path to the merged image CSV file. 
"""
    )
    
    # Retrieve the arguments from the user's input in command line
    args = parser.parse_args()
    
    # Map the user's inputs from command line to the appropriate function
    if args.modify_featureExtractionCSV:
        cell, module_7b_output_dir, output_dir = args.modify_featureExtractionCSV
        modify_featureExtractionCSV(cell, module_7b_output_dir, output_dir)