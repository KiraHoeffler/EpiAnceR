# EpiAnceR
A method to adjust for genetic ancestry in DNA methylation studies (run on Illumina arrays) when genotyping data is not available.

Ancestry_PCs:
- Ancestry_functions: ancestry_info() extracts the ancestry information from the DNA methylation dataset and ancestry_PCA() uses this information to calculate ancestry PCs
- Resources: SNP_cgs files contain CpGs that overlap with SNPs (0bp distance) at MAF < 0.05 for the different arrays

Manuscript:
- Ancestry_pipeline: pipeline used for the analyses and tests in the manuscript
- Resources - contains all files needed to run the ancestry PC calculation & follow-up analyses
        - SNP_cgs: CpGs that overlap with SNPs (0bp distance) at MAF < 0.05 for the different arrays -> needed for the ancestry calculations
        - General_functions: contains functions that are used for the ancestry PC calculations and analyses in the manuscript
        - theme: ggplot2 themes used to generate figures for the manuscript


Contact:
kira.hoeffler@uib.no
