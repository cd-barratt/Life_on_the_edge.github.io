We have recently developed a novel informatic toolbox applicable to any species or geographic area to predict vulnerability to global change, ultimately aimed at improving conservation prioritisation efforts. 

The toolbox facilitates the integration of environmental (e.g. climate, land use), ecological (e.g. spatial occurrences), and evolutionary (e.g. genome-wide SNP) data via a series of modular scripts. The toolbox can be run from start to finish (i.e. raw spatial, environmental and genomic data) through to final population vulnerability maps), or specific modules can be used separately (e.g. if you just want to build Species Distribution Models, look at population striucture or perform Genotype Environment Association analyses).

We follow the frameworks of two main papers, [Razgour et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12694) and [Razgour et al. 2019](https://www.pnas.org/doi/10.1073/pnas.1820663116), via a series of scripts and functions that have been generalised with flexible code to accomodate any species input data from any geographic region. The toolbox estimates three main metrics for each population/sampling locality:

* EXPOSURE - the magnitude of predicted change at future climatic conditions
* SENSITIVITY - the neutral and adaptive capacity of each individual/population based on genomic and environmental data
* RANGE SHIFT POTENTIAL - the current and predicted change in landscape connecivity for future climatic conditions

These three metrics are then assessed to calculate final POPULATION VULNERABILITY (see gure below for workflow and final output)

<img src="https://cd-barratt.github.io/Life_on_the_edge.github.io/workflow.png"  align="left" width="600">
<img src="https://cd-barratt.github.io/Life_on_the_edge.github.io/pop_vulnerability.png"  align="right" width="600">

The toolbox runs from a params.tsv file (up to 40 parameters which may be defined/modified), and all you need to provide are the spatial and genomic input data (though you can filter the input spatial, environmental and molecular data that are included based on your own requirements).

Below the main functionailty of the toolbox is briefly listed:
* Download and process of genome-wide data (e.g. from SRA or ENI, or your own raw data)
* Perform QC on these data to help select the most appropriate processing parameters in Stacks 2 for each dataset
* Optimise and finalise parameters following [best practices](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) and generate PLINK format output files ready for the Life on the edge toolbox
* Download and process spatial data (e.g. GBIF and georeferenced genomic data), clean it and prepare environmental data (e.g. Workdclim) for Species Distribution Modelling (including clipping to relevant extents per species and preparing background pseudoabsence data)
* Run ensemble species distribution models and predict future changes based on user specified parameters
* Calculate the magnitude of environmental change predicted to occur for each population (based on environmental dissimilarity and SDMs)
* Perform Genotype-Environment Association Analyses with user specified environmental predictors (RDA, LFMM)
* Individually categorise individuals in GEA ordination space to assess the extent of local adaptation within each sampled population
* Calculate neutral and adaptive genetic diversity per population
* Build adaptive SDMs using these locally adapted individuals to gain a more thorough understanding of differential responses amongst locally adapted populations
* Evaluate the current and future landscape connectivity of populations (using Circuitscape) to assess the potential for evolutionary rescue of isolated populations
* Combine multiple analyses to integrate EXPOSURE, adaptive and neutral SENSITIVITY and RANGE SHIFT POTENTIAL to make POPULATION VULNERABILITY predictions
* Create summary PDFs with completely transparent log files recording all steps within the toolbox

