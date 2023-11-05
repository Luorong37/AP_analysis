# AP Analysis Toolkit

This toolkit is designed to process and analyze action potential (AP) data from neuronal recordings. It offers a suite of MATLAB functions to load raw data, select regions of interest (ROIs), apply photobleaching correction, find peaks, analyze sensitivity and signal-to-noise ratio (SNR), and calculate firing rates.

## Getting Started

Clone this repository to your local machine using:

git clone https://github.com/Luorong37/AP_analysis.git


### Prerequisites

Ensure you have MATLAB installed with the following toolboxes:

- Image Processing Toolbox
- Curve Fitting Toolbox
- Signal Processing Toolbox

### Installation

No additional installation is required. Ensure you set your MATLAB path to include the directory where you cloned this repository.

### Usage

Run the main script `AP_Analysis.m`. The script is structured as follows:

1. Loading raw data: Specify the path to your .tif, .tiff, or .bin files containing neuronal recording data.
2. Saving the loaded movie: Optionally save the loaded data to a .mat file for quicker loading in future sessions.
3. Creating and saving a sensitivity map: This is an optional step to visualize the sensitivity of your recorded data.
4. Selecting ROIs: Manually select ROIs for further analysis.
5. Correcting for photobleaching: Apply a highpass filter to correct for photobleaching in the fluorescence signal.
6. Finding peaks: Use defined thresholds and window widths to find peaks representing action potentials.
7. Sensitivity and SNR analysis: Calculate and plot the sensitivity and SNR of the detected action potentials.
8. Statistic analysis of APs: Generate statistics for
