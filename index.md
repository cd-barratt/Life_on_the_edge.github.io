## Life on the edge

### A new informatic tolbox to quantify intraspecific vulnerability to global change
We have developed a novel informatic toolbox applicable to any species or geographic area to predict vulnerability to global change, aimed at improving conservation prioritisation efforts.

The toolbox facilitates the integration of environmental (e.g. climate, land use), ecological (e.g. spatial occurrences), and evolutionary (e.g. genome-wide SNP) data via a series of modular scripts. The toolbox can be run from start (i.e. raw spatial, environmental and genomic data) to finish (i.e. population vulnerability maps), or specific modules can be used separately (e.g. if you just want to build Species Distribution Models, look at population striucture or perform Genotype Environment Association analyses.

We follow the frameworks of two main papers, [Razgour et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12694) and [Razgour et al. 2019](https://www.pnas.org/doi/10.1073/pnas.1820663116), and a series of 22 scripts have been generalised with flexible code to accomodate any species input data from any geographic region.

Here's an overview of the main parts of the pipeline:
<img src="https://cd-barratt.github.io/Life_on_the_edge.github.io/workflow.png"  align="center" width="1700">

To briefly summarise, the toolbox can do the following things, running from a params.tsv file (up to 36 parameters which may be defined)
* Download annd process genome-wide data (e.g. from SRA or ENI, or your own raw data)
* Perform QC on these data to help select the most appropriate processing parameters in Stacks 2 for each dataset
* Finalise the parameters following [best practices](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) and generate PLINK format output files for the Life on the edge toolbox
* Download and process spatial data (e.g. GBIF and georeferenced genomic data), clean it and prepare environmental data (e.g. Workdclim) for Species Distribution Modelling (including clipping to relevant extents per species and preparing background pseudoabsence data)
* Run ensemble species distribution models and predict future changes based on user specified parameters
* Calculate the magnitude of environmental change predicted to occur for each population (based on environmental dissimilarity and SDMs)
* Perform Genotype-Environment Association Analyses with user specified environmental predictors (RDA, LFMM)
* Individually categorise individuals to assess the extent of local adaptation within each given sampled population
* Calculate neutral and adaptive genetic diversity per population
* Build refined SDMs using these locally adapted individuals to gain a more thorough understanding of differential responses amongst locally adapted populations
* Evaluate the current and future landscape connectivity of populations (using Circuitscape) to assess the potential for evolutionary rescue of isolated populations
* Combine multiple analyses to integrate Exposure, Adaptive and Neutral sensitivity and Range shift potential to make final population vulnerability predictions
* Create completely transparent logs and summary PDFs of all steps within the toolbox

Example output:
<img src="https://cd-barratt.github.io/Life_on_the_edge.github.io/pop_vulnerability.png"  align="center" width="1700">
