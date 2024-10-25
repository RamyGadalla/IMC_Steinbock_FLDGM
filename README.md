# IMC Steinbock FLDGM

## Overview

This repository contains scripts and resources for working with Imaging Mass Cytometry (IMC) data using Steinbock and related tools. The project involves preprocessing, segmentation, and analysis of IMC images, leveraging both Python and R scripts, as well as HistoCAT and CellProfiler pipelines.

## Repository Structure

- **Python Scripts**: Scripts for data preprocessing and analysis.
- **R\_spatial**: Contains R scripts related to spatial analysis of IMC data.
- **CellProfiler Pipelines (`.cppipe`): Workflows for preparing and segmenting images.
- **CSV Files**:
  - `images.csv`: List of images and metadata.
  - `panel.csv`: Details about markers used in the IMC.
- **Shell Script (`steinbock.sh`): Script to automate Steinbock processing steps.
- **RDS File (`spe_steinbock.rds`): R data object for spatial experiment data.
- **Ilastik Classifier (`pixel_classifier.ilp`): Pixel classification model for Ilastik.

## Usage

1. **Preprocessing**: Use `1_prepare_ilastik.cppipe` to prepare images for Ilastik.
2. **Segmentation**: Run `2_segment_ilastik.cppipe` to segment images after pixel classification.
3. **Analysis**: Use Python and R scripts for further analysis, including spatial and statistical evaluations.
4. **Automation**: Execute `steinbock.sh` to streamline the processing workflow.

## Requirements

- **Python** (version >= 3.8)
- **R** (version >= 4.0)
- **CellProfiler**
- **Ilastik**


![](https://github.com/RamyGadalla/IMC_Steinbock_FLDGM/blob/main/Tissue_archit_graph.png)

##
