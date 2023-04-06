## How to use the toolbox

### Life on the edge – a new toolbox for climate change vulnerability assessments at the population level: Vignette
Christopher D. Barratt, Renske E. Onstein, Malin Pinsky, Sebastian Steinfartz, Hjalmar Kuehl, Brenna R. Forester

## Overview
Life on edge (hereafter LotE) is a new climate change vulnerability assessment toolbox, facilitating the integration of environmental, molecular and ecological data. With the increasing availability of high-quality georeferenced genome-wide datasets published in open access online repositories, as well as constantly improving climate model simulations, the LotE framework offers a range of tools that can be used to investigate intraspecific responses to global change, thus providing empirical results from large genomic and spatial datasets to inform and assist biodiversity conservation in our rapidly changing world. The toolbox uses the concepts defined in Razgour et al. (2018) based on the IPCC AR5 report (2014) to assess the Exposure, Sensitivity and Range shift potential 
of intraspecific populations across species, creating a vulnerability index per population that can be compared within and across species to identify the early warning signals of potential population declines due to global change.

Steps 1-5 below provide details on the initial setup of the toolbox and guidelines for formatting the underlying datasets to analyse. Steps 6-12 walk the user through a typical complete LotE analysis. Step 13 demonstrates how the toolbox can be used modularly so that only requested parts of the analyses will be undertaken, for example if no genomic data is available for your species and you wish to analyse exposure (the magnitude of predicted change) and range shift potential (ability of populations to move) together.

## Setup
### 1.	Installation
To use LotE you’ll first need to do the following:
1)	Download the github directory for LotE, and move it to where you want to run the toolbox from (in your HPC environment)
2)	Install Circuitscape, Stacks 2, Singularity, and Julia (or ask your system administrators very nicely)
3)	Install R version 4.1.3 or download the Singularity container containing that version of R. Move the singularity container to your LotE working directory
4)	Run the 00_setup_life_on_the_edge.sh script in your HPC environment. The .sh scripts throughout LotE are designed for HPC job queue systems using SLURM, if you use SGE/UGE systems or otherwise then please talk to your HPC cluster administrator to translate them. You’ll need to define your working directories, emails, job logs in the submit script, as well as providing links to your R libraries. The script will install all necessary R packages and dependencies in your containerized version of R that runs in Singularity. 
5)	Open R and install package dependencies (00_setup.R). Many dependencies are required in R and it’s important that none have error messages (which is why the Singularity container may be useful here, this is tried and tested!)

We highly recommend using LotE in a HPC cluster environment and via the suppled Singularity contained due to the heavy computational load of many of the toolbox modules, particularly if pre-processing raw genomic data using Stacks 2, or even if dealing with processed genomic datasets with larger numbers of individuals/SNPs. Species Distribution Modelling and Circuitscape analyses can also be particularly computationally demanding if environmental data is at high spatial resolution and/or spanning large geographic areas (i.e. most of the time).

### 2.	Preparing input data 
Genomic input data should be a PLINK formatted .map and .ped file, with each row representing an individual, and each column representing a SNP. Individual names of each sample should be listed in the second column, containing no whitespaces. If processing genomic data yourself, Stacks 2 (and other data processing programs) have a tendency to output a header line in the .map and .ped files containing some info on the program used to generate the files, this will need to be deleted before using the files for LotE. The .map and .ped files should be placed in a folder named your species (‘Genus_species’) separated by an underscore (e.g. Homo_sapiens) in /-data-/genomic_data.

Spatial data should be a .csv file with three columns; first column named ‘Sample’ should contain the individual names. These names should exactly match the individuals in the genomic .map and .ped files. The second and third columns (‘LONG’, ‘LAT’) should be populated with the georeferenced coordinates of where the individual is from, in decimal degrees. The .csv file should be placed in a folder named your species (Genus_species) separated by an underscore (e.g. Homo_sapiens) in /-data-/spatial_data.

### 3.	Preparing environmental data 
Environmental data should be downloaded for the present and future time periods in order to run several parts of the LotE toolbox and make future predictions (i.e. SDMs, Exposure, investigating local adaptation with GEAs). We recommend freely available high-resolution data (ideally at 30 arc seconds, ~1km2 resolution) from Worldclim2 or CHELSA. The environmental data for the time period of future projections is your choice (e.g. 2070), and you also may select a specific global change scenario of your choice (e.g. Shared Socioeconomic Pathway SSP5 – worst case scenario) for which you can either download the full dataset and store it, or crop it to a smaller region (see script 00_process_environmental_data.R). 

The environmental data for present and future should be stored in their own directories (‘current’ and ‘future’), with separate files representing the data for each predictor (e.g .tif or .asc files). It’s important that the same predictors are available for both time periods, and are at the same spatial resolution and geographic extent. The folder locations of these data are defined in the params file (‘current_environmental_data_path’ and ‘future_environmental_data_path’). 

### 4.	Populating the params file
The params file (params.tsv) will control parameters for all analyses. The params file contains a row per species, and up to 41 parameters that can be used, each line of the params file is thus independent for each species analysis. Most parameters are essential to specify, so the toolbox will fail without them (e.g. the species name, ‘species_binomial’), but some are optional (e.g. mapping extent, ‘geographic_extent’), see Table 3 in the manuscript for an overview of which params are optional. We recommend best practice of populating almost all parameters so that you can be sure your analysis will not fail, or at least copying the example params file and modifying it to suit your own species.

### 5.	Ensuring LotE knows where to find your scripts and programs
Assuming that you did all the above steps correctly, you are almost ready to begin analysing data. First there are a couple of final things to do:

i) Edit your 01_run_life_on_the_edge.sh submit script to match your own HPC details and setup:
*	In the header (first 9 lines) you’ll need to modify the email address and output directory
*	In the main body of the script (lines 11 onward) you’ll need to call your version of PLINK (line 13), your version of Julia (lines 16-18) and your version of Java (lines 21-23). Check with your HPC cluster administrator if you are unsure how to do this correctly, as this will be different for every user
*	Also in the main body of the script, change the working directory (line 26) and the export link to local user R libraries – this will ensure that your R libraries are correctly exported

ii) Edit your run_LOE_exposure.R, -LFMM-.R, run_LOE_sensitivity.R, run_LOE_range_shift_potential.R, run_LOE_population_vulnerability.R scripts (in -scripts-/) so that the first line points towards the directory where your toolbox is located. You shouldn’t need to modify anything else in these scripts unless you want to selectively choose which parts of the pipeline are run (see step 13).




