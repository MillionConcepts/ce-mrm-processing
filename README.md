# ce-mrm-processing

## introduction

This repository contains software and ancillary materials used to create the PDS4 bundle
`chang_e_microwave_processed`. This bundle includes Planetary Data System 4 (PDS4) versions 
of calibrated (L2C) data products from the Chinese Lunar Exploration Program (CLEP)'s 
Chang'e-1 (CE-1) and Chang'e-2 (CE-2) Microwave Radiometer (MRM) instruments, along 
with new products derived from those data. It is hosted by the 
[Geosciences Node of the Planetary Data System](https://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-chang_e_microwave_processed/).

This repository is not a mirror of that bundle, which is too large to distribute on 
GitHub. This repository is intended principally as an extension of the documentation 
provided in that bundle.

Please refer to the in-code documentation for complete descriptions of the software's
functionality and intended usage.

For a more detailed discussion of the software's context and background, please refer
to `primary_bundle_documentation.md` in the root directory of this repository.

## repository contents

### /label_templates

Label templates (marked-up XML, not syntactically valid until processed) used to 
programmatically generate PDS4 labels.

### /label_templates/stubs

Text files used to help populate some templates.

### /mrm

A Python module containing all executable code. The root directory of the module contains
handler scripts used to execute various steps of the processing and labeling pipeline. 
Some are essentially independent of one another, and others must be run in order. Please 
refer to their in-code documentation for usage instructions and configuration options.

### /mrm/converter

Code for programmatically constructing PDS4 labels from templates.

### /mrm/modeling

Code for modeling physical and brightness temperatures.

### /mrm/modeling/heat1d

An implementation of a planetary heatflow model based on [K.-Michael Aye's Python version 
of Paul Hayne's `heat1d` package](https://github.com/phayne/heat1d/). 

### /mrm/processing

Analysis and data manipulation code. Includes modules for deconvolving and map-projecting 
point samples, ingesting source data products, constructing tables, writing browse products,
and more.

### /mrm/shared

Assorted utility code used by other submodules of `mrm`.

### /static/fonts

Font files used in browse product generation.

## code usage

All executable code in this repository is standard Python. The provided environment.yml
file may be used to construct a conda environment suitable for running all included code.

Some steps of the processing and labeling pipeline require source data files that are 
too large to distribute via GitHub. Please email mstclair@millionconcepts.com if you 
would like copies of these files.

NOTE: The parallel processing features of some portions of the code may not be fully 
compatible with MacOS or Windows.

## licensing

All templates and Python files in this repository are subject to the BSD 3-Clause License,
copyright Million Concepts.

The code in /mrm/processing/heat1d/ is a rewrite of [K.-Michael Aye's Python version of 
Paul Hayne's `heat1d` package](https://github.com/phayne/heat1d/), and, as such, is 
additionally subject to the MIT License, copyright Paul Hayne (included in this repository's
LICENSE file).

The files in /static/fonts are subject to the Open Font License under multiple copyrights.
Copies of their licenses are included in /static/fonts.

`primary_bundle_documentation.md` is in the public domain.
