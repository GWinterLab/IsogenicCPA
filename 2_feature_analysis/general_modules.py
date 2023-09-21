#############################################
# !! WORKING VERSION !!
# Author(s): Amanda Ng R.H.
# Created on: 31 Mar 2023
# Last updated on: 28 Jul 2023 (added a function for adjusting colors)
# Documentation status:
# In progress
# See https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html
# for more information on the Sphinx documentation style.
#
# General functions
#############################################
# Import package(s)
import os
from matplotlib import pylab as plt
import matplotlib.colors as mc
import colorsys

# Function for making a directory if it does not exist already
def makeDirectory(directory):
    """
    Function that checks if the user-provided path to a directory exists. If it does not exist, the
    directory will be made.
    
    :parem directory: Path to a directory (existing/non-existing).
    :type directory: str
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"The following directory has been made:\n{directory}")
    else:
        print(f"The following directory already exists:\n{directory}")
    return

# Function for saving plots if the user has provided an output_path
def savePlot(output_path = ""):
    """
    Function that exports plots to a user-provided output path if provided.
    
    :param output_path: Path to output the plot to.
    :type output_path: str
    """
    if len(output_path) > 0:
        plt.savefig(output_path, dpi = 200, transparent = True, bbox_inches = "tight")
        print(f"EXPORTED: The plot has been exported at\n{output_path}")
    else:
        print("WARNING: The plot has NOT been exported.")
    return

# Function that darkens or lightens a given color
def adjust_colorLuminosity(color, adjustBy):
    """
    Adjusts the given color by multiplying (1-luminosity).
    Input can be matplotlib color, hex string or RGB tuple.
    From: https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return(colorsys.hls_to_rgb(c[0], 1 - adjustBy * (1 - c[1]), c[2]))