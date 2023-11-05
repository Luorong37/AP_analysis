# AP Analysis Workflow Documentation

## Overview

This workflow is designed for analyzing action potential (AP) recordings using MATLAB. It includes scripts for loading data, selecting regions of interest (ROIs), correcting for photobleaching, finding peaks, and conducting sensitivity and signal-to-noise ratio (SNR) analysis. It concludes with statistical analysis and visualization of AP characteristics and firing rates.

Developed by: Liu-Yang Luorong

## Getting Started

### Prerequisites

To run this workflow, ensure that you have the following MATLAB toolboxes installed:

- Image Processing Toolbox
- Curve Fitting Toolbox
- Signal Processing Toolbox

### Installation

Clone the repository from [GitHub](https://github.com/Luorong37/AP_analysis) or download the ZIP file and extract it to your local machine.

### Usage

1. Define the path to your data folder and the frame rate in the script.
2. Run the script section-by-section by pressing `Ctrl+Enter` in MATLAB, observing prompts, and interacting with figures as necessary.

## File Structure

- `load_movie`: Function to load raw imaging data.
- `select_ROI`: Function to select regions of interest (ROIs) on the data.
- `calculate_firing_rate`: Function to calculate the firing rate based on detected spikes.
- `calculate_FWHM`: Function to calculate the full width at half maximum for APs.
- `calculate_SNR`: Function to calculate signal-to-noise ratio.

## Workflow Steps

1. **Loading raw data:** Load imaging files or raw data using provided dialogues.
2. **Creating a map:** Generate sensitivity maps to visualize ROI correlation coefficients.
3. **Selecting ROIs:** Manually select ROIs for further analysis.
4. **Correcting Photobleaching:** Apply single exponential function to correct for photobleaching artifacts.
5. **Finding Peaks:** Determine the peaks in the fluorescence traces that correspond to APs.
6. **Analyzing Sensitivity and SNR:** Calculate sensitivity and SNR for each trace.
7. **Statistical Analysis:** Generate statistics on AP characteristics across all ROIs.
8. **Plotting:** Visualize average sensitivity, SNR, and firing rate data for each ROI.

## Output

The script will generate figures and save them in a specified analysis folder. Data points such as firing rates, peak amplitudes, and SNR are also saved in .mat and .xlsx formats.

## Support

For issues or questions about the workflow, open an issue on the GitHub repository or contact the maintainer directly.

## Acknowledgments

Thanks to the contributors and users of the AP_analysis toolset for feedback and contributions to the development of this workflow.

