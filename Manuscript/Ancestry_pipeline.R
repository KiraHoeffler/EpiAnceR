
# Author: Kira Hoeffler
# Aim: script to test the performance of different methods to adjust for ancestry
# please let me know if you struggle with any part of the script

################################################################################
# LOAD LIBERIES
################################################################################

# Please make sure you have the packages below installed before you start running the script. 

suppressMessages({
  # LOAD LIBRARIES
  
  # These packages need to be loaded manually on our server, probably not a problem on your server:
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  #library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  #library(shiny, lib.loc = "C:/Users/kihof7027/AppData/Local/R/win-library/4.4")
  
  library(ggplot2)
  library(readxl)
  library(writexl)
  library(plotly)
  library(htmlwidgets)
  library(tidyverse)
  library(stats)
  library(cowplot)
  library(car)
  library(rstatix)
  library(effectsize)
  library(cluster)
  library(factoextra)
  library(minfi)
  library(wateRmelon)
  library(ChAMP)
  library(limma)
  library(patchwork)
})



################################################################################
# SET-UP (NEEDS TO BE ADAPTED!!)
################################################################################

# SET WORKING DIRECTORY
# IMPORTANT: Make sure the Resources folder we provide is in the working directory!
dir_gen <- "S:/Project/WP-epigenetics/08_Ancestry/Ancestry_Pipeline/" #working directory
setwd(dir_gen)

# SAMPLESHEET PATH AFTER QC (samples that did not pass the QC are excluded)
# IMPORTANT: your samplesheet should be an Excel file
# it should have the following columns:
#     <anc_classif>: with genetically predicted ancestry (or reported ancestry, preferably short labels)
#     <Basename> with the sample name (overlapping with the sample names in the RGset)
#     <Indiv_ID> individual ID if you have multiple samples per individual
#     Sex
#     Age
#     columns containing estimated cell type proportions, for example:
#             blood: "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"
#             saliva: "Epi", "Fib", "comb_ICs"
samplesheet_path <- "S:/Project/WP-epigenetics/08_Ancestry/samplesheet_after_QC.xlsx"

# PATH TO BACKGROUND-CORRECTED AND EXTENDED RGSET
#     -> background-corrected using bgcorrect.illumina(RGset))
#     -> extended using read.metharray(basenames, extended = TRUE) or similar function
RGset_path <- "S:/Project/WP-epigenetics/04_Pipeline_old/output/RData/RGset_bgcorrected.Rdata"

# DO YOU HAVE MULTIPLE SAMPLES PER INDIVIDUAL? ("yes" or "no")
multiple_samples_per_indiv <- "yes"

################################################################################
# SET UP FOLDER STRUCTURE FOR THE OUTPUT FOLDER
################################################################################

if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/RData")) dir.create("output/RData")
if (!dir.exists("output/Figures")) dir.create("output/Figures")
if (!dir.exists("output/Tables")) dir.create("output/Tables")

################################################################################
# IMPORT
################################################################################

# FUNCTIONS
source("Resources/General_functions.R")

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel(samplesheet_path))
samplesheet_test(samplesheet)

# IMPORT BACKGROUND-CORRECTED RGSET
# IMPORTANT: please call the object RGset!
# if you have not background-corrected your RGset yet, please use: RGset <- bgcorrect.illumina(RGset))
load(RGset_path)

# IF ALREADY AVAILABLE: DETECTION P VALUES (goes faster later)
load("S:/Project/WP-epigenetics/04_Pipeline_old/output/RData/detectionPvalue.RData")

################################################################################
# ADAPT SAMPLESHEET
################################################################################

# REMOVE ROWS WITH INCOMPLETE ANCESTRY INFORMATION
    # this can be adapted based on your dataset, please also exclude individuals with the label "mixed"
samplesheet_compl <- samplesheet[which(!is.na(samplesheet$anc_classif) & samplesheet$anc_classif != "missing" & samplesheet$anc_classif != "unknown"), ]

print("Please check that there is no missing ancestry information:")
table(samplesheet_compl$anc_classif)


# ORDER LEVELS OF ANCESTRIES ALPHABETICALLY
samplesheet_compl$anc_classif <- factor(samplesheet_compl$anc_classif, levels = sort(unique(samplesheet_compl$anc_classif)))


# OVERLAP FILEs
overlap_samples <- intersect(colnames(RGset), samplesheet_compl$Basename)
print(length(overlap_samples)) #check how many samples overlap
samplesheet_compl <- samplesheet_compl[which(samplesheet_compl$Basename %in% overlap_samples), ]
RGset <- RGset[, colnames(RGset) %in% overlap_samples]
dp <- dp[, colnames(dp) %in% overlap_samples] #detection p values (if you have them already)


# ORDER IMPORT FILES

# samplesheet:
samplesheet_compl <- samplesheet_compl[order(samplesheet_compl$Basename), ]

# RGset:
sample_names <- order(sampleNames(RGset))
RGset <- RGset[, sample_names]

# detection p values (if you have them already)
dp <- dp[, sort(colnames(dp))]



# SAVE FINAL SAMPLESHEET
write_xlsx(samplesheet_compl, path = "output/Tables/Samplesheet_final.xlsx")


################################################################################
# RUN PIPELINE TO GET ANCESTRY INFORMATION
################################################################################

# PLEASE ADAPT THE NAME OF YOUR RGset OBJECT, THE ARRAY TYPE, AND THE COLUMNS OF THE CELL TYPE PROPORTIONS IN THE FUNCTION

# array_type: "450K" OR "EPICv1" OR "EPICv2"
# cell_cols: column names of cell type proportions in the samplesheet, e.g. for saliva c("Epi", "Fib", "comb_ICs")
# dp = detection p value matrix (if not available, just leave it out. Then it just takes a lot longer to run the pipeline.)
ancestry_info(background_corr_RGset = RGset, samplesheet = samplesheet_compl, array_type = "EPICv2", path_SNP_folder = "Resources/", 
              cell_cols = c("Epi", "Fib", "comb_ICs"), output_folder = "output/RData/", dp = dp)


rm(RGset)

# LOAD ANCESTRY DATA
load("output/RData/comb_SNPs_NewMethod.RData")
load("output/RData/impSNPs_OldMethod.RData")
load("output/RData/resid_beta_NewMethod.RData")
  

################################################################################
# PERFORMANCE PIPELINE
################################################################################

# MAKE LIST OF ALL THE SNP ANCESTRY OBJECTS YOU WANT TO TEST AND NAME THEM
# please write it in the order you want to plot it later
# please do not test more than 5 methods in parallel
SNP_list <- list(beta_SNP0bp_imp, resid_beta, comb_SNPs)
SNP_list_names <- c("EpiAncOrig", "EpiAnceR", "EpiAnceR+")

# RUN PIPELINE
performance_pipeline(SNP_list = SNP_list, SNP_list_names = SNP_list_names, samplesheet_compl = samplesheet_compl, output_folder = "output/")
