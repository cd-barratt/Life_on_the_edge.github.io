We have recently developed a novel informatic toolbox applicable to any species or geographic area to predict vulnerability to global change, ultimately aimed at improving conservation prioritisation efforts. 

The toolbox facilitates the integration of environmental (e.g. climate, land use), ecological (e.g. spatial occurrences), and evolutionary (e.g. genome-wide SNP) data via a series of modular scripts and functions. The toolbox can be run from start to finish (i.e. raw spatial, environmental and genomic data) through to final population vulnerability maps), or specific modules can be used separately (e.g. if you just want to build Species Distribution Models, look at population striucture or perform Genotype Environment Association analyses).

We follow the frameworks of two main papers, [Razgour et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12694) and [Razgour et al. 2019](https://www.pnas.org/doi/10.1073/pnas.1820663116), via a series of scripts and functions that have been generalised with flexible code to accomodate any species input data from any geographic region. The toolbox estimates four main metrics for each population/sampling locality:

* EXPOSURE - the magnitude of predicted change at future climatic conditions (i.e. environmental dissimilarity + changes in habitat suitability using species distribution models)
* NEUTRAL SENSITIVITY - the neutral sensitivity of each individual/population to global change based on genomic data (i.e. lower heterozygosity = higher sensitivity)
* ADAPTIVE SENSITIVITY - the adaptive capacity of each individual/population to global change based on genomic and environmental data (i.e. how locally adapted populations are)
* RANGE SHIFT POTENTIAL - the potential for populations to be able to move (i.e. based on current landscape connectivity)

These four metrics are then assessed to calculate final POPULATION VULNERABILITY (see figures below for workflow and final output)

![image](https://cd-barratt.github.io/Life_on_the_edge.github.io/workflow_general.png)

![image](https://cd-barratt.github.io/Life_on_the_edge.github.io/workflow.png)

![image](https://cd-barratt.github.io/Life_on_the_edge.github.io/pop_vulnerability.png)

The toolbox runs from a Params.tsv file (up to 53 parameters which may be defined/modified), and all you need to provide are the spatial, enviromental and genomic input data per species

Below a list of some of the things you can do using the Life on the edge toolbox:
* Download and process of genome-wide data (e.g. from SRA or ENI, or your own raw data)
* Perform QC on these data to help select the most appropriate processing parameters in Stacks for each dataset
* Optimise and finalise parameters following [best practices](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) and generate PLINK format output files ready for the Life on the edge toolbox
* Download and process spatial data (e.g. GBIF and georeferenced genomic data), clean it and prepare environmental data (e.g. Worldclim) for Species Distribution Modelling (including clipping to relevant extents per species, spatially rarefying input data to control for spatial autocorrelation, and preparing background pseudoabsence data)
* Run ensemble species distribution models and predict future changes based on user specified parameters
* Calculate the magnitude of environmental change predicted to occur for each population (based on environmental dissimilarity and SDMs)
* Perform Genotype-Environment Association Analyses with user specified environmental predictors (RDA, LFMM)
* Individually categorise individuals in GEA ordination space to assess the extent of local adaptation within each sampled population
* Calculate neutral and adaptive genetic diversity per population based on standing genetic diversity (heterozygosity) and local adaptation of individuals and populations
* Build adaptive SDMs using these locally adapted individuals to gain a more thorough understanding of differential responses amongst locally adapted populations
* Evaluate the range shift potential of populations (using Circuitscape) to assess the potential for evolutionary rescue of isolated populations and spread of neutral/adaptive genetic diversity
* Combine multiple analyses to integrate EXPOSURE, adaptive and neutral SENSITIVITY and RANGE SHIFT POTENTIAL to make POPULATION VULNERABILITY predictions
* Create summary PDFs with completely transparent log files recording all steps within the toolbox

