# barseq-Leishmania-2025
A repository containing all scripts and data for the 2025 paper

### Scripts

Scripts for the analysis included within the paper above.

Processing folder contains the initial barcode extraction Rust script, used on the raw FASTQ files, and the subseqnent R script for counting the reads and creating the data table for downstream analysis. These scripts were developed with Alastair Droop and are also included in his University of York Data Science repository at: https://github.com/uoy-research/fqbarcode alongside further relevant information. 

Barseq folder contains all analysis scripts for the analysis of the barcode counts dataset, covering Figures 2 to 6 and Extended Data Figures 3-7.

In_vitro folder contains the R scripts for analysis of the in vitro infection data (Poisson fitting, histogram plots and Poisson plots) in Figure 1. It also contains the in vitro parasite culture data analysis R script, also in Figure 1. 

IVIS folder contains the R script needed to analysis the ex vivo IVIS data and create the figures in Figure 6 relating to this. 

Simulation folder contains the MATLAB scripts for running the simulation described in Figure 1. There are five scripts, one for the initial simulation, one for the simulation with growth rate variation and three for the three clonal expansion scenarios. Also included is the R script for the correlation analysis performed on the simulatied datasets and for producing the figure in the paper. 

Validation folder contains the R scripts for the analysis discussed in the Validation section of Methods and shown in Supplementary Figures 4-6. For data re-analysis, we used pubished data from Hotinger et al., Hullahalli and Waldor and Lebrun-Corbin et al. In the case of Hotinger et al, this data was obtained directly from the paper authors via email correspondence. In the case of Hullahalli and Waldor, the data was the 1_to_54_OrderedFrequencies csv file available with the eLife manuscript. In the case of Lebrun-Corbin et al., the data was the Frequencies_M1toM13 csv file available with the mBio manuscript. 

### Data files

Data files for the analysis performed in the paper are found in the Data_files subfolder.

Barseq contains the read counts for the barcode sequencing data comprising most of the paper data (Figures 2-6 and Extended Data Figures 3-7).

In_vitro contains the in vitro parasite library culture data and the in vitro macrophage infection datasets, covered in Figure 1. 

IVIS contains the IVIS ex vivo imaging for the secondary infection tissue burdens, covered in Figure 6 and Extended Data Figure 8.

Simulation contains the MATLAB simulation output datasets, covered in Figure 1.


### References

Hotinger, J. A., Campbell, I. W., Hullahalli, K., Osaki, A. & Waldor, M. K. Quantification of Salmonella enterica serovar Typhimurium population dynamics in murine infection using a highly diverse barcoded library. eLife 13, RP101388 (2025).

Hullahalli, K. & Waldor, M. K. Pathogen clonal expansion underlies multiorgan dissemination and organ-specific outcomes during murine systemic infection. eLife 10, e70910 (2021).

Lebrun-Corbin, M. et al. Pseudomonas aeruginosa population dynamics in a vancomycin-induced murine model of gastrointestinal carriage. mBio 16, e0313624 (2025).
