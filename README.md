# barseq-Leishmania-2025
A repository containing all scripts and data for the 2025 paper "Bi-directional highways and super-seeder tissues underpin parasite dissemination in experimental visceral leishmaniasis", by Loughrey et al. All R scripts were developed and run using R version 4.4.3. 

### Scripts

Scripts for the analysis included within the paper above.

Processing folder contains the initial barcode extraction Rust script, used on the raw FASTQ files, and the subseqnent R script for counting the reads and creating the data table for downstream analysis. These scripts were developed with Alastair Droop and are also included in his University of York Data Science repository at: https://github.com/uoy-research/fqbarcode alongside further relevant information. The output data from the Rust script is required in order to run the R script. Raw FASTQ data is available at NCBI BioProject with Project ID PRJNA1392034 (http://www.ncbi.nlm.nih.gov/bioproject/1392034). Dependency details for the Rust script are as follows:
* Rust version = "0.2.1"; edition = "2021"
* Dependencies: regex = "1.7.1"; rand = "0.8.5"; simple-eyre = "0.3.1"; clap = { version = "4.1.10", features = ["derive", "usage"] }; stderrlog = "0.5.4"; log = "0.4.17"; flate2 = "1.0.25"; levenshtein = "1.0.5"

Barseq folder contains all analysis scripts for the analysis of the barcode counts dataset, covering Figures 2 to 6 and Extended Data Figures 3-7. The scripts titled "GD_calculation_script", "TSFP_FFP_calculations_script" and "Shannon_index_analysis_script" need to be run before the other scripts, as the others rely on outputs from these.

In_vitro folder contains the R scripts for analysis of the in vitro infection data (Poisson fitting, histogram plots and Poisson plots) in Figure 1. It also contains the in vitro parasite culture data analysis R script, also in Figure 1. The "Infection_Poisson_fitting_analysis" script needs to be run before the "Plotting_Poisson_parameters" script. 

IVIS folder contains the R script needed to analysis the ex vivo IVIS data and create the figures in Figure 6 relating to this. 

Simulation folder contains the MATLAB scripts for running the simulation described in Figure 1. There are five scripts, one for the initial simulation, one for the simulation with growth rate variation and three for the three clonal expansion scenarios. Also included is the R script for the correlation analysis performed on the simulatied datasets and for producing the figure in the paper. If you do not wish to run the MATLAB code, the generated CSV files are available in the Data_files folder, so you can run just the "Analysis_script" on these files instead. 

Validation folder contains the R scripts for the analysis discussed in the Validation section of Methods and shown in Supplementary Figures 4-6. For data re-analysis, we used published data from Hotinger et al., Hullahalli and Waldor and Lebrun-Corbin et al. In the case of Hotinger et al, this data was obtained directly from the paper authors via email correspondence. In the case of Hullahalli and Waldor, the data was the 1_to_54_OrderedFrequencies csv file available with the eLife manuscript. In the case of Lebrun-Corbin et al., the data was the Frequencies_M1toM13 csv file available with the mBio manuscript. The other two scripts "Migration_validation_TSFP_FFP_script" and "Randomised_simulations_validation_script" use the 2025-JCM-002-barcodes.csv dataset in the Barseq folder inside Data_files.

### Data files

Data files for the analysis performed in the paper are found in the Data_files subfolder.

Barseq contains the read counts for the barcode sequencing data comprising most of the paper data (Figures 2-6 and Extended Data Figures 3-7).

In_vitro contains the in vitro parasite library culture data and the in vitro macrophage infection datasets, covered in Figure 1. 

IVIS contains the IVIS ex vivo imaging for the secondary infection tissue burdens, covered in Figure 6 and Extended Data Figure 8.

Simulation contains the MATLAB simulation output datasets, covered in Figure 1.

### Notes

The scripts should download as R files and run in any application which can run R-based code. 

Please ensure you add your working directory details to the scripts before running. 

You will need to save the datasets into a folder called Data (or alter the scripts accordingly) to run the scripts without errors. These should all save as CSV files when downloaded. 

Time for downloading scripts and data is minimal on a standard computer. Time to run all of the code in full should amount to less than four hours on a standard computer. 

Please ensure you run scripts in the correct order, as some scripts depend on results from prior scripts in order to run. See comments above and at start of individual scripts for details.

R may ask for permission to create relevant directories to save results into whilst running these scripts. 

We used the following packages for our analysis:

* tidyverse 2.0.0
* ggraph 2.2.1
* ggpubr 0.6.0
* igraph 2.1.4
* fitdistrplus 1.2-2
* Seurat 5.3.99.9000
* rstatix 0.7.2
* cplm 0.7-12.1
* tweedie 2.3.5
* vegan 2.6-10

We cannot validate that the code will run correctly on other versions without testing.

### References

Hotinger, J. A., Campbell, I. W., Hullahalli, K., Osaki, A. & Waldor, M. K. Quantification of Salmonella enterica serovar Typhimurium population dynamics in murine infection using a highly diverse barcoded library. eLife 13, RP101388 (2025).

Hullahalli, K. & Waldor, M. K. Pathogen clonal expansion underlies multiorgan dissemination and organ-specific outcomes during murine systemic infection. eLife 10, e70910 (2021).

Lebrun-Corbin, M. et al. Pseudomonas aeruginosa population dynamics in a vancomycin-induced murine model of gastrointestinal carriage. mBio 16, e0313624 (2025).
