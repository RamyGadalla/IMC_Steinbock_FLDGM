#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:38:44 2022

@author: 
"""

#remember to set the env to 'imcsegpipe'
#$ conda activate imcsegpipe - in terminal. this environment is in 'anaconda3/envs'
#should later move to the project folder.


# Load necessary modules
from pathlib import Path
from urllib import request

#create path object using pathlib module. ".." means parent directory
raw_folder = Path("../raw")
#make a directory in the specified path
raw_folder.mkdir(exist_ok=True, parents=True)

for example_file_name, example_file_url in [
    (
        "Patient1.zip",
        "https://zenodo.org/record/5949116/files/Patient1.zip",
    ),
    (
        "Patient2.zip",
        "https://zenodo.org/record/5949116/files/Patient2.zip",
    ),
    (
        "Patient3.zip",
        "https://zenodo.org/record/5949116/files/Patient3.zip",
    ),
    (
        "Patient4.zip",
        "https://zenodo.org/record/5949116/files/Patient4.zip",
    ),
    (
        "panel.csv",
        "https://zenodo.org/record/5949116/files/panel.csv",
    )
]:
    example_file = raw_folder / example_file_name
    if not example_file.exists():
        request.urlretrieve(example_file_url, example_file)


# Ilastik project
ilastik_project = Path("..") / "IMCWorkflow.ilp"
if not ilastik_project.exists():
    request.urlretrieve("https://zenodo.org/record/6449127/files/IMCWorkflow.ilp", ilastik_project)

# Sample metadata
sample_metadata = Path("..") / "sample_metadata.xlsx"
if not sample_metadata.exists():
    request.urlretrieve("https://zenodo.org/record/5949116/files/sample_metadata.xlsx", sample_metadata)
