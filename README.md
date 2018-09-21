# Wintertime_melt

This folder contains scripts used for processing automatic weather station (AWS) data and Met Office Unified Model data using python.

Some of this work was published in Kuipers Munneke et al. (2018) Geophysical Research Letters 45, 7615–7623. doi: http://doi.org/10.1029/2018GL077899. 
The remainder comprises work completed for the first part of my PhD thesis.

Observational data: 

AWS data are taken from the Cabinet Inlet AWS (AWS18), operated and managed by IMAU. The dataset is publicly available online as "Dataset from iWS 18 in Cabinet Inlet, Larsen C Ice Shelf, Antarctica, 2014–2017". It measures surface meteorological variables and surface energy balance parameters. 
* Link to metadata and information: https://www.projects.science.uu.nl/iceclimate/aws/files_oper/oper_29121.php

* Link to data: http://doi.org/10.5285/05c9124b-7119-4d99-8e17-ab754eb3f51c

+++++++++++

Model data:
Code is written to process model data from the Met Office Unified Model (vn10.4), in their proprietary format, .pp, which results in a smaller file size. However, this could very easily be adapted to work with netCDF files, and model output from any version of the model from vn7 onwards (provided the right parameters are available). The scripts require the installation of the python package Iris (https://scitools.org.uk/iris/docs/latest/), which is a tool for processing CF-compliant climate and model data. 