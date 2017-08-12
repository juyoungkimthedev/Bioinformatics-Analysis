# Bioinformatics-Analysis
Typical Bioinformatics Gene Expression Analysis using R

Analysis of publicly available Microarray Data - Bioinformatics Workflow
(this is focus on a re-analysis of public data)

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Introduction

Although microarrays have been superseded by high-throughput sequencing technologies for gene expression
profiling such as RNASeq, years of experience gained from analysing microarray data has led to a 
variety of analysis techniques and datasets that can be exploited in other contexts. In this section
we will focus on retrieving and exploring microarray data from public repositories such as GEO

Aims: 
1. Importing Data from public repositories
2. Exploratory data analysis techniques for high-throughput data
3. Constructing pipeline workflows for the analysis of Affymetrix gene expression data
4. Preprocessing and normalization of gene expression data
5. Visualization of gene expression data and identification of clusters of similar genes/drugs
6. Principal component analysis and hierachical clustering of gene expression data
7. Differential Gene Expression Analysis using linear-modeling techniques
8. Biological pathway identification/interpretation of the results w.r.t gene expression data

Objectives:
What will you do to achieve the aims?
a. import gene expression datasets from 'The Connectivity Map Project - Broad Institute'
b. Set up experimental conditions, biological replicates necessary for the analysis
c. Assess the quality of dataset in a repository and undergo necessary preprocessing steps
d. Identify, and correct for batch effects to remove any technical bias in the data
e. Construct 'workable' standard gene expression matrix (list of genes x therapeutic drugs). This also includes biological gene names annotation using necessary annotation libraries
f. Create a heat map to visualize the regulations of genes against drugs applied and identify clusters of similar genes / acting drugs
g. Perform standard Differential Gene Expression Analysis to get a ranked list of 'differentially' expressed genes
h. Use un-supervised / clustering methods to explore the dataset
i. Interrogate particular genes of interest and identify its biological pathways and interpret in terms of Breast Cancer related carbohydrate gene perspectives

Prerequisites
1. Must be able to code in R programming language and utilize number of computational biology related
2. packages from the Bioconductor project in Cran. Available: source("http://www.bioconductor.org/biocLite.R")
number of packages pre-installed must include: Affy, Biobase, limma, longevityTools, devtools (for GitHub installation) etc
3. Must be able to utilize the codes and work flow procedures suggested by the "LOngevityTools" project by Thomas Girke from University of Riverside California USA
"Longevitytools" project offer number of useful packages for microarray analysis and also some of their code snippets are very helpful
4. Must install annotation library packages for the Affymetrix gene chips such as HGU133A/HTHGU133A annotation libraries

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

A. Data and Script

Please download publicly available gene expression data in response to drug treatment from the Connectivity Map Project
this is available at: http://portals.broadinstitute.org/cmap/ or you can use T. Girke longevitytools project code to automate the download
the code is: getCmapCEL(rerun = FALSE). Note that this has to run in RStudio which is IDE for R programming development. User must define working directory in advance to have the data in the working directory specified.
Download of files will take some time

B. Overview of Microarray Technologies:
for theoretical background, only important concepts are covered as well as actual procedures followed in the M.Sc project 2016-2017
single channel microarrays can be produced to measure the absolute expression level of every genes of interest in a given sample
the fluorescence of each feature is a measure of the expression level of a particular gene. The most popular microarray technology was Affymetrix
Affymetrix use 25 base-pair probes that are synthesized on the array surface. Each gene of interest is interrogated by a collection of 11-20 probe pairs known as probe set
The expression level for a gene is then derived by combining all measurements from a particular probe set

Typical workflow:

Despite differences in array construction, there are a few commonalities in the way that raw data from a microarray experiment are processed
1. Image processing: microarray surface is scanned by a laser to produce an image representation of the fluorescence emitted by it.
These images are usually processed by the manufacturer's software. However, the pixel intensities measured on the image
may be influenced by factors other than hybridization etc. Therefore, a background intensity is estimated for each feature to account for such factors
The background and foreground estimates generally act as a starting point for statistical analysis

2. Data Processing: the intensities of the features on a microarray are influenced by many sources of noise and repeated measurements made on different microarrays may also appear to disagree
Therefore, a number of data-cleaning, or pre-processing steps, must take place before being able to draw a valid biological conclusions from the experiment
Background correction, quality assessments, transformation, normalization, annotation steps are necessary
refer to: https://rawgit.com/bioinformatics-core-shared-training/microarray-analysis/master/intro.html for more information

Annotation: microarray manufacturers use their own identifiers schemes that don't relate to biology such as '1000_at', '1001_at' ...
We need to map these IDs to gene names, genome position etc.
Sometimes the mappings can be wrong, this is the main reason why sequencing is better

Side note on Microarrays vs Sequencing

The 'death' of microarrays was predicted as early as 2008. We have recently reached the tipping point where RNA-Seq has taken over from gene expression arrays
Regardless, many of the same issues and techniques apply to NGS data
Experimental design, quality assessments, normalization, Statistics (testing for RNA-Seq is build upon the knowledge from microarrays)

The Bioconductor Project:
Packages analyse all kinds of Genomic Data, documentation, course materials, example data and workflow, re-usable frameworks.
i.e. citations for Bioconductor software packages: limma, Smyth (2005), 2417 (N)
Downloading a package: each package has its own landing page which includes all installation script, vignettes and manuals, details of package maintainer
reading data using Bioconductor: All data can be read into R using read.csv/read.table/read.delim but for this project we have read the .CEL files using affy package
Dataset can be split into different components such as: matrix of expression values, sample information, annotation for the probes
