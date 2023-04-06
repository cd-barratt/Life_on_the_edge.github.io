## How to use the toolbox: Life on the edge vignette

## Overview
Life on edge (hereafter LotE) is a new climate change vulnerability assessment toolbox, facilitating the integration of environmental, molecular and ecological data. With the increasing availability of high-quality georeferenced genome-wide datasets published in open access online repositories, as well as constantly improving climate model simulations, the LotE framework offers a range of tools that can be used to investigate intraspecific responses to global change, thus providing empirical results from large genomic and spatial datasets to inform and assist biodiversity conservation in our rapidly changing world. The toolbox uses the concepts defined in [Razgour et al. (2018)](https://doi.org/10.1111/1755-0998.12694) and [Razgour et al. (2019)](https://doi.org/10.1073/pnas.1820663116) based on the [IPCC AR5 report (2014)](https://www.ipcc.ch/report/ar5/syr/) to assess the Exposure, Sensitivity and Range shift potential of intraspecific populations across species, creating a vulnerability index per population that can be compared within and across species to identify the early warning signals of potential population declines due to global change

**Steps 1-5** below provide details on the initial setup of the toolbox and guidelines for formatting the underlying datasets to analyse. **Steps 6-12** walk the user through a typical complete LotE analysis. **Step 13** demonstrates how the toolbox can be used modularly so that only requested parts of the analyses will be undertaken, for example if no genomic data is available for your species and you wish to analyse exposure (the magnitude of predicted change) and range shift potential (ability of populations to move) together

## Setup
### 1.	Installation
To use LotE you’ll first need to do the following:
* Download the [github directory](https://github.com/cd-barratt/Life_on_the_edge) for LotE, and move it to where you want to run the toolbox from (in your HPC environment)
* Install [Circuitscape](https://circuitscape.org/downloads/), [Stacks 2](https://catchenlab.life.illinois.edu/stacks/), [Singularity](https://sylabs.io/singularity/), and [Julia](https://julialang.org/) (...or ask your system administrators very nicely)
* Install R version 4.1.3 or download the [Singularity container](https://cloud.sylabs.io/library/sinwood/bioconductor/bioconductor_3.14) containing that version of R. Move the singularity container to your LotE working directory
* Run the `00_setup_life_on_the_edge.sh` script in your HPC environment. The .sh scripts throughout LotE are designed for HPC job queue systems using SLURM, if you use SGE/UGE systems or otherwise then please talk to your HPC cluster administrator to translate them. You’ll need to define your working directories, emails, job logs in the submit script, as well as providing links to your R libraries. The script will install all necessary R packages and dependencies in your containerized version of R that runs in Singularity
* Open R and install package dependencies (`00_setup.R`). Many dependencies are required in R and it’s important that none have error messages (which is why the Singularity container may be useful here, this is tried and tested!)

We highly recommend using LotE in a HPC cluster environment and via the suppled Singularity contained due to the heavy computational load of many of the toolbox modules, particularly if pre-processing raw genomic data using Stacks 2, or even if dealing with processed genomic datasets with larger numbers of individuals/SNPs. Species Distribution Modelling and Circuitscape analyses can also be particularly computationally demanding if environmental data is at high spatial resolution and/or spanning large geographic areas (i.e. most of the time)

### 2.	Preparing input data 
Genomic input data should be a PLINK formatted *.map* and *.ped* file, with each row representing an individual, and each column representing a SNP. Individual names of each sample should be listed in the second column, containing no whitespaces. If processing genomic data yourself, Stacks 2 (and other data processing programs) have a tendency to output a header line in the *.map* and *.ped* files containing some info on the program used to generate the files, this will need to be deleted before using the files for LotE. The *.map* and *.ped* files should be placed in a folder named your species (‘Genus_species’) separated by an underscore (e.g. Homo_sapiens) in `/-data-/genomic_data/`

Spatial data should be a *.csv* file with three columns; first column named ‘Sample’ should contain the individual names. These names should exactly match the individuals in the genomic *.map* and *.ped* files. The second and third columns (‘LONG’, ‘LAT’) should be populated with the georeferenced coordinates of where the individual is from, in decimal degrees. The *.csv* file should be placed in a folder named your species (Genus_species) separated by an underscore (e.g. Homo_sapiens) in `/-data-/spatial_data/`

### 3.	Preparing environmental data 
Environmental data should be downloaded for the present and future time periods in order to run several parts of the LotE toolbox and make future predictions (i.e. SDMs, Exposure, investigating local adaptation with GEAs). We recommend freely available high-resolution data (ideally at 30 arc seconds, ~1km<sup>2</sup> resolution) from [Worldclim2](http://www.worldclim.com/version2) or [CHELSA](https://chelsa-climate.org/). The environmental data for the time period of future projections is your choice (e.g. 2070), and you also may select a specific [global change scenario](https://www.carbonbrief.org/explainer-how-shared-socioeconomic-pathways-explore-future-climate-change/) of your choice (e.g. Shared Socioeconomic Pathway SSP5 – worst case scenario) for which you can either download the full dataset and store it, or crop it to a smaller region (see script `00_process_environmental_data.R`)

The environmental data for present and future should be stored in their own directories (`current` and `future`), with separate files representing the data for each predictor (e.g *.tif* or *.asc* files). It’s important that the same predictors are available for both time periods, and are at the same spatial resolution and geographic extent. The folder locations of these data are defined in the params file (‘**current_environmental_data_path**’ and ‘**future_environmental_data_path**’)

### 4.	Populating the params file
The params file, `params.tsv`, stored in the root of your main directory, controls relevant parameters for all analyses. The params file contains a row per species, and up to 41 parameters that can be used, each line of the params file is thus independent for each species analysis. Most parameters are essential to specify, so the toolbox will fail without them (e.g. the species name, ‘**species_binomial**’), but some are optional (e.g. mapping extent, ‘**geographic_extent**’), see Table 3 in the manuscript for an overview of which params are optional. We recommend best practice of populating almost all parameters so that you can be sure your analysis will not fail, or at least copying the example params file and modifying it to suit your own species

### 5.	Ensuring LotE knows where to find your scripts and programs
Assuming that you did all the above steps correctly, you are almost ready to begin analysing data. First there are a couple of final things to do:

i) Edit your `01_run_life_on_the_edge.sh` submit script to match your own HPC details and setup:
*	In the header (first 9 lines) you’ll need to modify the email address and output directory
*	In the main body of the script (lines 11 onward) you’ll need to call your version of PLINK (line 13), your version of Julia (lines 16-18) and your version of Java (lines 21-23). Check with your HPC cluster administrator if you are unsure how to do this correctly, as this will be different for every user
*	Also in the main body of the script, change the working directory (line 26) and the export link to your own user R libraries – this will ensure that your R libraries are correctly exported

ii) Edit your `run_LOE_exposure.R`, `-LFMM-.R`, `-RDA-.R`, `run_LOE_sensitivity.R`, `run_LOE_range_shift_potential.R`, `run_LOE_population_vulnerability.R` scripts (in `-scripts-/`) so that the first line points towards the directory where your toolbox is located. You shouldn’t need to modify anything else in these scripts unless you want to selectively choose which parts of the pipeline are run (Step 13)

## Analysing data
### 6.	Overview of analyses
Finally, you are ready for performing analyses! In the first part of this vignette we will perform a complete analysis of some example data for an East African spiny reed frog (*Afrixalus fornasini*) from start to finish. We provide example *.map*, *.ped*, *.csv* and *params.tsv* files as part of the LotE package. We will run through all steps of the toolbox one at a time and pause at certain parts to check outputs when decisions need to be made that affect subsequent steps. This example uses already newly processsed genomic data from [Barratt et al. (2018)](https://doi.org/10.1111/mec.14862), using Stacks 2,  with RAD-seq data for 7309 SNPs genotyped across 43 individuals

Best practice for using LotE on an unknown (novel) dataset are to run in parts and check outputs before running through the subsequent parts of the toolbox. For example, knowledge of population structure is needed in order to select a reasonable estimate of the number of populations (k) for Genotype-Environment Association (GEA) analyses, and if this is not correctly accounted for then GEA results may be confounded by structure in the data that is unaccounted for. Furthermore, some exploration of the GEA analyses themselves is necessary, where you will likely want to rerun them to explore the data and set a tolerance for true and false positives for your species

We follow this best practice here, running the following components of the LotE toolbox as separate shell scripts that call the relevant software, R scripts and functions. Essentially, all analyses will be passing the relevant parameters from the *params.csv* file to control each analytical step:

*	Running Exposure analyses, including Species Distribution Modelling (SDM) - [biomod2](https://github.com/biomodhub/biomod2)
*	Running Sensitivity analyses, incuding GEAs - [Latent Factor Mixed Models (LFMM)](https://doi.org/10.1093/molbev/msz008) and [RDA (Redundancy Analysis)](https://cran.r-project.org/web/packages/vegan/vegan.pdf)
*	Running Range shift potential analyses, including [Circuitscape](https://circuitscape.org/) modelling 
*	Running Population vulnerability analyses and summarising all results

All relevant details will be reported in the log file: `./-outputs-/log_files/Afrixalus_fornasini.log`

### 7.	Exposure
To run the exposure analyses, we will simply run the following code embedded in a shell script:

``` singularity exec ~/barratt_software/Singularity_container/bioconductor_3.14.sif Rscript ./-scripts-/run_LOE_exposure.R ‘Afrixalus_fornasini’ ``` 

To give you an idea of what this is doing - this will read the contents of the run_LOE_exposure.R script, running through each line in sequence. The script itself sources all the internal LotE functions on lines 3-4 and then calls them on each new line. Below, you can see that each function will be passed the species name (‘**species_binomial**’) to run the analyses, and each internal function will extract the parameters from the relevant line of the params file that matches the species name

``` 
prepare_spatial_data(species_binomial)
prepare_environmental_data(species_binomial)
spatially_rarefy_presences_and_create_background_data(species_binomial)
sdms_biomod2(species_binomial)
exposure(species_binomial)
impute_missing_data(species_binomial) 
```

This code will take a little while to run everything depending on your parameters, the longest step being the SDMs. Below is an overview of what is happening at each stage and which outputs are generated by each function/script: 

**prepare_spatial_data()** prepares the spatial data prior to building SDMs. The script will download existing [GBIF](https://www.gbif.org/) (Global Biodiversity Information Facility) data for your defined species, and clean it using the [CoordinateCleaner R package](https://doi.org/10.1111/2041-210X.13152) (Zizka et al. 2019). Data will be cleaned to remove any records without coordinates, those that are representing country centroids, biodiversity institutions (e.g. museums, university collections), and those that are geographic outliers from the rest of the species range (>1000km, this can be modified in the script). The cleaned GBIF data will be combined with the existing georeferenced genomic data and written as a new file ending in *_presence_data.csv* which will be used for spatially rarefying the presence data prior to building SDMs. The output file from this will look something like below (image truncated):
 
**prepare_environmental_data()** reads the environmental data stored in the current and future directories and crops it to the geographic extent for the target species. This is important to ensure that SDMs in particular are modelling an area that is ecologically relevant for the study species. If the variable '**geographic_extent**' is defined in the params file (xmin, xmax, ymin, ymax) then this will be taken as the geographic modelling extent, and if this is not populated, the function will take the extent covered by all presence samples written above by prepare_spatial_data() to define the modelling extent. Secondly, the function will then extract the relevant environmental data for all predictors at each sampling location of all individuals (= populations) and store it in a new file ending *_full_env_data.csv* which will be used for GEAs later. The output file from this will look something like below (image truncated):
 
**spatially_rarefy_presences_and_create_background_data()** reads in the *_presence_data.csv* written by **prepare_spatial_data()** and spatially rarefies the data so that spatial autocorrelation (e.g. sampling bias) will not confound predictions made by SDMs. The '**sp_rare_dist_km**' parameter in the params file represents the minimum distance that two presences are allowed to be, if two presences are less than this distance apart (i.e. highly spatially clustered) then one of the presences will be removed. The spatially rarefied presence data will be written as a new file (*_thinned.csv*) for building SDMs and plotted with the original data for visualization. The function will then generate a number of background (pseudoabsence) points from a given buffer around presence points (defined in the params file – '**n_background_points**' and '**buffer_distance_degrees**'), writing these points as a new file (*_background_points.csv*) and plotting them for visualization

The output files from this will look something like below (raw data vs. spatially rarefied data and a summary of background points), and the relevant *.csv* files will be stored for later use by the toolbox:
        
**sdms_biomod2()** will use the spatially rarefied presence data to build species distribution models (SDMs) for each species, for present conditions, but also for forecasted future conditions.  Here, the framework follows SDM best practices, see [Araujo et al. (2019)](DOI: 10.1126/sciadv.aat4858), using reduced spatial autocorrelation in presence data, as well as accounting for multi-collinearity in predictor variables using Variance Inflation Factors (VIF) if defined in the params file. Alternatively, users may specify a subset of predictors that are ecologically relevant for the species in question (‘**subset_predictors**’). The biomod2 R package is used to evaluate models built using available modelling algorithms which is also subsettable via the params file, ‘**biomod_algorithms**', and retaining only ‘good’ models (i.e. TSS>0.5, modifiable in the params file '**TSS_min**' or '**ROC_min**') for the final ensemble species distribution model prediction. If using MAXENT in your list of SDM modelling algorithms, a copy of the *maxent.jar* file from the [MAXENT website](https://biodiversityinformatics.amnh.org/open_source/maxent/) will need to be placed in a suitable location, and this location set as the ‘**maxent_path**’ parameter in the params file. The model will then be projected onto the future environmental data layers to forecast the species distribution in the future. Variable importances will also be tracked across all retained models so that the user has a sense of which predictor variables have the most influence on the species distribution. Several SDM parameters are modifiable in the params file, including if VIF should be used ('**perform_vif**'), or a subset of variables should be selected ('**subset_predictors**'), which SDM algorithms to use ('**biomod_algorithms**'), number of replicates per algorithm ('**sdm_reps_per_algorithm**', data split percentage between model training and testing ('**data_split_percentage**'). The function will create a self-contained analysis subfolder within the SDMs folder containing all the models and data, summarizing the models themselves as well as various metrics of their performance

The output files from the SDMs will look something like below (current and future SDM predictions, response curves per predictor, variable importances of model predictors and overall performance of different modelling algorithms that passed defined thresholds to generate final ensemble models)
   
**exposure()** firstly compares how selected climate variables change between the present and predicted future conditions. It takes this information as well as the predicted change in the SDM between the same two time periods and uses it to calculate ‘Exposure’ (ranging from 0-10), measuring the magnitude of change each population/locality will experience between the two time periods. The output files from the exposure analysis will look something like below (predicted change between current and future SDM predictions – range expansions in green and range contractions in orange, predicted change in environmental predictor 1, predicted change in environmental predictor 2)

From these data, a final *Exposure.csv* file will be generated that summarises all the changes in the SDM and your selected environmental predictors for each population (i.e. unique LONG/LAT)

**impute_missing_data()** A final step of the ‘Exposure’ submit script is to set up files for the next step (‘Sensitivity’). Firstly, because GEA analyses (LFMM and RDA) typically require complete data matrices of called SNP genotypes, something which is not achievable using reduced representation library approaches (e.g. RAD-seq/ddRAD-seq), it is necessary to impute missing data, which will be done in two ways, firstly using [LEA]((https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12382)) (Frichot & François 2015), and secondly based on the mean frequencies of genotypes per population cluster (see Razgour et al. 2018). This will generate additional files (*_imp.csv*, *_imp.geno* and *_imp.lfmm*) in the `-data-/genomic_data/` directory for use with the GEA methods in the next part of the toolbox. Lastly, some basic population structure analyses will be performed in order to understand the spatial population structure of the data so that the GEAs are accounting for this. This will be performed as part of this function and will be output in the `-outputs-/genomic_analyses/` folder – with a *_pop_structure_pca.png* file (including a PCA plot of the data and the inertia of each PC axis) as well as an *_snmf_barplot.png* and *_snmf_cross_entropy.png* file which together can be used to determine likely structure. The output files will look like below, with a PCA plot and inertia of each PC axis, a summary of the cross-entropy values for varying values of k specified in the params file, and a barplot of the population structure selected by the best k
   
Finally, a population assignment file is made (*_pop_assignment_LEA.csv*) which has the population assignment for each individual:

At this point it is important to update the params file for your species with the most appropriate number of population clusters represented in the data ('**k**').

### 8.	GEA- LFMM
To run the first part of the sensitivity analyses (LFMM), run the following code embedded in a shell script:

``` singularity exec ~/barratt_software/Singularity_container/bioconductor_3.14.sif Rscript ./-scripts-/-LFMM-.R ‘Afrixalus_fornasini’ ```

LFMM is a univariate method that will account for population structure from the underlying data, and SNP genotypes will be statistically evaluated against your defined environmental predictors to select candidate SNPs that are potentially under selection. LFMM will read in the environmental data you prepared using **prepare_environmental_data()**, and subset the variables of choice (defined as '**env_predictor_1**' and '**env_predictor_2**' in the params file). As highly colinear variables are problematic for GEA it will check the Variance Inflation Factor and report the correlation between the variables in a pair plot (*_env_correlations.png*). If your variables are highly correlated (e.g. >0.8) it may be worth considering alternative variables from your predictor set that are less strongly correlated

After this step, the data will be analysed using LFMM (reading your newly populated value of k from your params file). LFMM will read the imputed LFMM formatted file (the one generated by **impute_missing_data()**), and perform a GEA to search for candidate SNPs that show strong statistical correlations with your environmental predictors. For each SNP, a p-value will be assigned, which will be used to calculate a q-value. This q-value will ultimately be used to decide which candidate SNPs pass the pre-defined FDR (False Discovery Rate) thresholds of 0.1, 0.05 and 0.01. Firstly the p-values and calibrated p-values will be calculated, and then the GIF (Genomic Inflation Factor) will be used to rescale the p-values into an acceptable distribution. We recommend that the '**scale_gif_lfmm**' parameter in params be set to 1 (i.e. no scaling) on first running LFMM, and then some investigation is required by the user. The LFMM script will summarise the p-values, the calibrated p-values and the scaled p-values so that you can look at their distributions. Ideally you want a histogram distribution more leaning towards the left (i.e. low p-values), but not too liberal (i.e over-representative high frequencies of low p-values). The plot should be investigated to ensure an acceptable p value distribution (i.e. not to conservative, not too liberal). See [this great tutorial](https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html) on GIF modifcation in LFMM for more details

Output plot of p-values, calibrated p-values and re-adjusted p-values. Here the re-adjusted p-values after modifying the GIF using the scale_gif_lfmm parameter in the params file look acceptable, so we will go with this for now

Ideally, you should experiment with the GIF, running LFMM a few times until you are satisfied with the readjusted p-value distributions (lower panel). In the log file, the GIF is reported each time, so you can use this to adjust the '**scale_gif_lfmm**' parameter which modifies the GIF. A well calibrated set of p-values should show a GIF of around 1, too liberal <1 and too conservative >1, so if your GIF is around 1.8 for example, and your p-values are very skewed towards high values, you could set your '**scale_gif_lfmm**' parameter to 0.7 which would bring the newly calculated GIF to 1.26 (=1.8 x 0.7) and increase the frequency of lower p-values. Each time, your list of candidate SNPs that are selected to be below the defined FDR (False Discovery Rate) thresholds (0.1, 0.05, 0.01) will change, and it is worth keeping in mind that only a fraction of your total loci (in this case the total is 7309) should realistically show signals of local adaptation. Thus, histogram distributions that are too liberal will detect high numbers of false positives, and too conservative approaches will result in zero detections. P-values are often not well behaved in empirical datasets so you should modify the GIF to an extent that the numbers of candidate SNPs and their p-value distribution is tolerable for you and believable for your study species – there is no right or wrong way to do this, it is subjective, and as long as you report your criteria exactly it is perfectly acceptable

The output files containing your candidate SNPs will be written per predictor and also summarized in *_LFMM_candidate_SNPs.csv* (truncated file below). The numbers of SNPs below each FDR threshold will be reported in the log file

After LFMM has completed and you are satisfied with your candidate SNPs, we will proceed with running RDA before performing the rest of the ‘Sensitivity’ analysis, which is somewhat less convoluted

### 9.	GEA- RDA
To run the second part of the sensitivity analyses (RDA), run the following code embedded in a shell script:

``` singularity exec ~/barratt_software/Singularity_container/bioconductor_3.14.sif Rscript ./-scripts-/-RDA-.R ‘Afrixalus_fornasini’ ```

Similar to LFMM, RDA (Redundancy analysis), a multivariate method, will be implemented using the vegan package in R. The process is similar to LFMM whereby population structure will be accounted for in the underlying data, and SNP genotypes will be statistically evaluated against your defined environmental predictors to select candidate SNPs that are potentially under selection. The RDA method implemented here applies the same framework as LFMM, whereby the GIF can be adjusted (using the '**scale_gif_rda**' parameter) to select thresholds for candidate SNPs, but it will also select outlier candidate SNPs using a function that measures the standard deviation of each SNP from the mean loading value across all SNPs (Razgour et al. 2019). We recommend setting this standard deviation ('**rda_sd**') in the params file to 3, but you can make this threshold less conservative by reducing it (e.g. 2.5)

Once the RDA script has run you can investigate the outputs; of particular interest here is the first plot (*_RDA_plot.png*) which will show in the left panel the clustering of individuals in ordination space and their position relative to the biplot arrows of each predictor, and in the right panel will show the eigenvalues of the PC axes

Secondly, just like with LFMM we can see the histogram of the unadjusted and adjusted p-values (after modifying the GIF) of the SNPs, but this time for each predictor (*_RDA_p_value_distribution_env_1.png*, *_RDA_p_value_distribution_env_2.png*). Like with LFMM, the distributions of these p-values can be modified by updating the GIF ('**scale_gif_rda**') in the params file and rerunning the script 

The outlier approach to selecting candidate SNPs is also adjustable as mentioned, outputting a two panel (i.e. a panel for each environmental predictor) histogram of the loadings of each SNP in the dataset. The outer head and tails of the distribution are classed as outlier candidate SNPs (defined using the '**rda_sd**' parameter in params)

An output plot of the SNPs in the RDA ordination space coloured by their environmental predictor association will be made (outlier approach using standard deviation)

Similarly, an output plot of the SNPs in the RDA ordination space with FDR<0.05 will be made, as well as a list of all candidate SNPs across the different FDR thresholds (<0.1, <0.05, <0.01) for each predictor
 
### 10.	Sensitivity
If you are satisfied with the LFMM and RDA analyses (having explored how changing the GIF and/or standard deviation parameters affects the output candidate SNPs), we can continue with the rest of the analysis. To continue with the rest of the ‘Sensitivity’ analysis, we submit the following code embedded in a shell script:

``` singularity exec ~/barratt_software/Singularity_container/bioconductor_3.14.sif Rscript ./-scripts-/run_LOE_sensitivity.R ‘Afrixalus_fornasini’ ```

Again, to elaborate what this is doing - this will read the contents of the `run_LOE_sensitivity.R` script, running through each line in sequence. Again it will pass the species name to each of the functions and run them line by line:

```
gea_rda_individual_categorisation(species_binomial)
adaptive_diversity(species_binomial)
neutral_diversity(species_binomial)
sensitivity(species_binomial)
create_circuitscape_inputs(species_binomial) 
```

**gea_rda_individual_categorisation(species_binomial)** categorises each individual based on their relative position in the constrained RDA ordination space relative to the environmental predictors (see Razgour et al. 2019). The script automatically reads in putative candidate SNPs that have been identified by either/both LFMM and RDA (definable in params file, see below) and then performs a new RDA based on only these candidate SNPs and the environmental data

To select which SNPs you wish to retain for the individual categorization analysis (we recommend using option 1 as the midpoint between being conservative and liberal),  set the '**which_loci**' parameter in the params file to one of the following:

*  1: present in either RDA (SD<3) or LFMM (FDR<0.05)
*  2: present in either RDA (FDR<0.05) or LFMM (FDR<0.05)
*  3: only LFMM loci (FDR<0.05)
*  4: only RDA loci (SD<3)
*  5: only RDA loci (FDR<0.05)
*  6: present in both RDA (SD<3) and LFMM (FDR<0.05)
*  7: present in both RDA (FDR<0.05) and LFMM (FDR<0.05) - not recommended, often 0 loci

It will write these loci to a file named *_adaptive_loci.txt* (below) and use this as a basis to categorise individuals using RDA

The built-in functions automatically select individuals that are falling within the range certain conditions following Razgour et al. (2019). For example, in a two predictor model (e.g. rainfall and temperature), certain individuals may be classed as 'hot_dry' or 'cold_wet' adapted. Individuals that fall in between the ordination categorisations of either of these categories may be categorised as ‘intermediate’. The distribution of these categorised individuals is plotted in the ordination space, and each individual, which category it belongs to, and its geographic location is saved in a separate *.csv* file. Categories (e.g. ‘hot_dry’/’cold_wet’) may be defined in the params file ('**category_1**', '**category_2**')

**adaptive_diversity()** quantifies the adaptive diversity across all populations by plotting the proportions of individuals that are adapted to each category defined in the params file (e.g. 'hot_dry'/'cold_wet'/'intermediate'). It uses the outputs generated by the individual categorisation (above), and plots a summary map of all samples (individuals and populations) and which conditions they are adapted to

Individual categorization summary
 
Population categorization summary (only sites with > individuals)
 
Output categorization maps (individual left panel, population right panel)
 
**neutral_diversity()** calculates neutral (i.e. non-adaptive) genetic diversity by masking out the putatively adaptive loci. It requires [PLINK](https://www.cog-genomics.org/plink/) (Purcel et al. 2007) to be installed (read the binary location from the params file, ‘**plink_executable**’), then calls PLINK via R. PLINK will generate the output files and the script here will automatically count the populations, number of individuals and neutral heterozygosity to calculate ‘Neutral sensitivity’ (ranging from 0-10)

Example output heterozygosity file *_neutral_sensitivity.csv*
 
**sensitivity()** integrates the calculated neutral diversity and exposure with the proportions of adaptively categorised individuals across each population. Based on these it will generate an ‘Adaptive sensitivity’ metric based on the % predicted change in each predictor and the proportion of individuals locally adapted to that predictor. The adaptive sensitivity metric will range between 0-10; being lower if a population has many individuals adapted to a condition that is forecast to change minimally, and higher if the conditions are forecast to change substantially

**create_circuitscape_inputs()** will setup the necessary files and structure for the range shift potential analyses. Because the subsequent analyses will use Circuitscape to model pairwise connectivity, memory and time requirements can become substantial with large datasets (i.e many populations/localities). For this reason it is highly recommended that these are performed in a HPC environment. This script will write the required files (*.ini*, *.jl*) for Circuitscape to run in Julia, along with the necessary files (*.sh*) to submit these jobs via HPC. The script copies relevant files to the Circuitscape directory for analysis, and generates the specific required file for the spatial points input (derived from the input spatial genomic data), which is used to specify the ‘nodes’ to model connectivity between populations. In the params file, an option to transform all 0’s to values of 0.001 is provided (‘**circuitscape_transform_zeros**’) (recommended if running range shift potential analyses on SDM outputs) so that Circuitscape does not interpret unsuitable areas as completely impermeable barriers

### 11.	Range shift potential
To run the the Range shift potential analyses (Circuitscape) run the following code embedded in a shell script:

``` julia --startup-file=no '/work/barratt/Life_on_the_edge_pipeline/-outputs-/'$1'/Range_shift_potential/circuitscape/'$1'_ensemble_SDM_current.jl'
julia --startup-file=no '/work/barratt/Life_on_the_edge_pipeline/-outputs-/'$1'/Range_shift_potential/circuitscape/'$1'_ensemble_SDM_future.jl'
singularity exec ~/barratt_software/Singularity_container/bioconductor_3.14.sif Rscript ./-scripts-/run_LOE_range_shift_potential.R 'Afrixalus_fornasini'  ```

The circuitscape analyses will take a long time to run, especially if there are many samples over a large geographic area (it will make pairwise comparisons across all sampling localities). The circuitscape analyses here are for analysing connectivity through the SDM, though you can (and probably should) explore other possible drivers of gene flow such as environmental predictors, forest cover etc.). A widely used R package for optimizing resistance surfaces is [ResistanceGA](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12984) (Peterman, 2018).

Once the Circuitscape analyses are complete, the R function **range_shift_potential()** will automatically summarise and plot the data, as well as quantifying range shift potential. **range_shift_potential()** reads in the output data from the Circuitscape analyses and plots maps of modelled present and future connectivity. It extracts the mean connectivity of each population to all other populations within a maximum dispersal distance (defined for each species in the params file) for both time periods, calculates the % change in connectivity and then uses this to calculate ‘Range shift potential’ (0-10), high reduction in range shift potential from present-future = 10, minimal reduction, stable or increase in range shift potential = 0. When calculating mean connectivity, a maximum dispersal distance in kilometres ('**max_dispersal_distance_km**') may be set in the params file for defining which populations are within geographic reach of one another (i.e. this avoids unrealistically distant populations being considered when calculating mean connectivity for each population). 
 
### 12.	Population_vulnerability
To run the final part of the toolbox, run the following code embedded in a shell script:
``` singularity exec ~/barratt_software/Singularity_container/bioconductor_3.14.sif Rscript ./-scripts-/run_LOE_population_vulnerability.R ‘Afrixalus_fornasini’ ```

To elaborate the last part of this code, it will do the same as previously, running the **population_vulnerability()** function to create the final files, and then using **summary_pdfs()** will generate the final *.pdf* output which captures all the relevant outputs and information in a single pdf

```
population_vulnerability(species_binomial)
summary_pdfs(species_binomial) 
``` 

**population_vulnerability()** integrates the exposure, adaptive and neutral sensitivity and range shift potential results to create a summary ‘Population vulnerability’ metric per population ranging from low (0) to high (10). The function will create a *.csv* output file of the observed data (i.e. where genomic samples are located, and use an inverse distance weighted interpolation to predict Exposure, neutral sensitivity, adaptive sensitivity and Range shift potential results across geographical space, with the underlying assumption that locations closer together share similar properties than those further apart. The function produces summary maps of non-interpolated (observed) and interpolated (predicted) data for each of the metrics (separately and also a composite 4-panel map). The params file allows the user to decide how population vulnerability is calculated by weighting the exposure, adaptive and neutral sensitivity and range shift potential metrics

**summary_pdfs()** uses all outputs generated and information in the log file to paste results together into a final summary PDF sheet using the [grobblR](https://cran.r-project.org/web/packages/grobblR/vignettes/grobblR.html) package (Floyd, 2020). Results can be identified and probed by rerunning the individual functions with modified parameter settings, and we recommend thorough reporting and transparency in all publications that use this toolbox.

13.	Partial analyses (running modularly)

