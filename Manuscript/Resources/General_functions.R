
# Author: Kira Hoeffler
# Functions to be sourced by the Performance_test script to test the performance of different ancestry-adjustment methods
# IMPORTANT: many of the functions or building on each other, load them together


################################################################################
# DEFINE THEME AND COLOURS

th <- theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(size = 14, colour = "black", face = "bold"),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
            axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
            axis.title.x = element_text(vjust = -0.5, size = 12, colour = "black"),
            axis.title.y = element_text(vjust = -0.5, size = 12, colour = "black", margin = margin(r = 10)),
            axis.title = element_text(face = "bold"),
            panel.border = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size =12),
            legend.key = element_blank())

th_transparent <- theme(panel.background = element_rect(fill = "transparent"),
                        plot.background = element_rect(fill = "transparent", colour = NA),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.background = element_rect(fill = "transparent", colour = NA),
                        legend.box.background = element_rect(fill = "transparent"))


# COLOUR PALETTE FOR PLOTS
my_colours <- c("#084081", "#126605", "#d19600", "#542788", "darkred")



#####################################################################################
                        # FUNCTIONS
####################################################################################

#' TEST IF SAMPLESHEET HAS THE NEEDED COLUMNS
#' @param samplesheet needs to have the columns Basename, Indiv_ID, and anc_classif

samplesheet_test <- function(samplesheet){
  if (!"Basename" %in% colnames(samplesheet)){
    stop("The samplesheet does not have the needed Basename column. Please check the comments")
  }
  if (!"anc_classif" %in% colnames(samplesheet)){
    stop("The samplesheet does not have the needed anc_classif column. Please check the comments")
  }
  if (multiple_samples_per_indiv == "yes"){
    if (!"Indiv_ID" %in% colnames(samplesheet)){
      stop("The samplesheet does not have the needed Indiv_ID column. Please check the comments")
    }
  }
  print("The samplesheet contains the needed columns.")
}

##################

#' SCALE TO 0-1
#' @param x vector
#' @return scaled vector (0-1)

scale_min_max <- function(x) {
  (x- min(x)) / (max(x) - min(x))
}

#################

#' CALCULATE EUCLIDEAN DISTANCE
#' @param a a vector of coordinates in 2D or 3D (e.g. X1, Y1, Z1)
#' @param b a vector of coordinates in 2D or 3D (e.g. X2, Y2, Z2)
#' @return euclidean distance

euclidean_distance <- function(a,b){
  sqrt(sum((a-b)^2))
}



##################


#' ADD COLOURS FROM COLOUR PALETTE TO PLOTS ACCORDING TO THE NUMBER OF methodS THAT WERE TESTED
#' @param methods vector with the names of the methods that were tested
#' @param use_case "fill" for scale_fill_manual in plot, "group" for scale_colour_manual in plot
#' @return line of code that can be added when setting up a ggplot 
#' 
add_costum_colours <- function(methods, use_case){
  if (length(methods) > 5) {
    stop("The functions to generate the figures work only for up to 5 tested methods in the same plot.")
  }
  
  selected_colours <- my_colours[1:length(methods)]
  
  if (use_case == "fill"){
    return(scale_fill_manual(values = selected_colours))
  }
  
  if (use_case == "group"){
    return(scale_colour_manual(values = selected_colours))
  }
  
}


####################

#' EXTRACT ANCESTRY INFORMATION
#' 
#' from DNA methylation of CpGs that overlap with SNPs
#' residualised to correct for technical and biological factors
#' combined with called genotypes from rs probes
#'
#' @param background_corr_RGset background-corrected RGset (using bgcorrect.illumina(RGset))
#' @param samplesheet containing samples that survived the QC; df; samples as rows, variables as columns; needs to have columns: 
#' Basename (=sample name also used in the RGset), Sex, Age, and columns with the cell type proportions
#' @param array_type "450K" OR "EPICv1" OR "EPICv2"
#' @param path_SNP_folder path to directory that contains the provided files: "SNP_cgs_EPICv2.RData", "SNP_cgs_EPICv1.RData", "SNP_cgs_450K.RData"
#' @param cell_cols column names of cell type proportions in the samplesheet
#' e.g. for saliva: c("Epi", "Fib", "comb_ICs")
#'
#' @export the ancestry df gets exported automatically
#' 
ancestry_info <- function(background_corr_RGset, samplesheet, array_type, path_SNP_folder, cell_cols, output_folder, dp = NA){
  
  ################################################################################
  # 0) TESTS
  ################################################################################
  
  print("Note: use the samplesheet that contains only samples that survived the overall QC.")
  print("Tests before starting the pipeline")
  
  
  # PACKAGES
  if (!requireNamespace("minfi", quietly = TRUE)){
    stop("The package minfi is not installed. Please install it before proceeding.")
  }
  if (!requireNamespace("wateRmelon", quietly = TRUE)){
    stop("The package wateRmelon is not installed. Please install it before proceeding.")
  }
  if (!requireNamespace("ChAMP", quietly = TRUE)){
    stop("The package ChAMP is not installed. Please install it before proceeding.")
  }
  
  # SAMPLESHEET
  if (class(samplesheet) != "data.frame"){
    stop("Check that the samplesheet object is a data.frame.")
  }
  if (!"Basename" %in% colnames(samplesheet)){
    stop("The <Basename> column is missing in the samplesheet. It refers to the sample name (and the colnames of the RGset)")
  }
  if (!"Sex" %in% colnames(samplesheet)){
    stop("The <Sex> column is missing in the samplesheet.")
  }
  if (!"Age" %in% colnames(samplesheet)){
    stop("The <Age> column is missing in the samplesheet.")
  }
  
  cell_cols_in_samplesheet <- cell_cols %in% colnames(samplesheet)
  if (FALSE %in% cell_cols_in_samplesheet){
    stop("The cell columns are not all in the samplesheet.")
  }
  rm(cell_cols_in_samplesheet)
  
  # ARRAY TYPE
  if (!array_type %in% c("450K", "EPICv1", "EPICv2")){
    stop("The array_type has to be <450K>, <EPICv1>, or <EPICv2>")
  }
  
  files_in_dir <- list.files(path_SNP_folder, recursive = TRUE, full.names = FALSE)
  
  if (array_type == "EPICv2" & !"SNP_cgs_EPICv2.RData" %in% files_in_dir){
    stop("The provided SNP files are not in the path_SNP_folder.")
  }
  if (array_type == "EPICv1" & !"SNP_cgs_EPICv1.RData" %in% files_in_dir){
    stop("The provided SNP files are not in the path_SNP_folder.")
  }
  if (array_type == "450K" & !"SNP_cgs_450K.RData" %in% files_in_dir){
    stop("The provided SNP files are not in the path_SNP_folder.")
  }
  
  # SAMPLE NAMES
  if ((length(intersect(colnames(background_corr_RGset), samplesheet$Basename)) < 5)){
    stop("The colnames of the RGset and the names in samplesheet$Basename do not overlap.")
  }
  
  
  print("All tests passed.")
  
  
  ### COMPLETE SAMPLESHEET ###
  
  samplesheet_compl <- samplesheet[which(!is.na(samplesheet$Age) & !is.na(samplesheet$Sex) & !is.na(samplesheet$Basename)), ]
  
  for (cell_col in c(1:length(cell_cols))){
    samplesheet_compl <- samplesheet_compl[which(!is.na(samplesheet_compl[, cell_col])), ]
  }
  
  if(nrow(samplesheet) != nrow(samplesheet_compl)){
    print("WARNING!! The Baseline and/or Age and/or Sex and/or cell proportion columns in the samplesheet have missing values.")
    print("Samples with missing values will not be included in the ancestry prediction")
  }
  
  
  ################################################################################
  # 1) FILTER AND ADAPT DATA FRAMES
  ################################################################################
  
  print("Step 1 - filter and adapt dataframes")
  
  # ORDER samplesheet_compl
  samplesheet_compl <- samplesheet_compl[order(samplesheet_compl$Basename), ]
  
  
  # LOAD SNP PROBES
  if (array_type == "EPICv2"){
    load(paste0(path_SNP_folder, "SNP_cgs_EPICv2.RData"))
  }
  if (array_type == "EPICv1"){
    load(paste0(path_SNP_folder, "SNP_cgs_EPICv1.RData"))
  }
  if (array_type == "450K"){
    load(paste0(path_SNP_folder, "SNP_cgs_450K.RData"))
  }
  
  
  # ANNOTATE RGset
  if (array_type == "EPICv2"){
    array <- "IlluminaHumanMethylationEPICv2" #array
    anno_type <- "20a1.hg38"
  }
  if (array_type == "EPICv1"){
    array <- "IlluminaHumanMethylationEPIC" #array
    anno_type <- "ilm10b4.hg19" #select type of annotation
  }
  
  if (array_type == "450K"){
    array <- "IlluminaHumanMethylation450k" #array
    anno_type <- "ilmn12.hg19" #select type of annotation
  }
  
  background_corr_RGset@annotation <- c(array = array, annotation = anno_type)
  
  
  
  
  # OVERLAP FILEs
  overlap_samples <- intersect(colnames(background_corr_RGset), samplesheet_compl$Basename)
  samplesheet_compl <- samplesheet_compl[which(samplesheet_compl$Basename %in% overlap_samples), ]
  background_corr_RGset <- background_corr_RGset[, colnames(background_corr_RGset) %in% overlap_samples]
  
  # ORDER IMPORT FILES
  
  # samplesheet:
  samplesheet_compl <- samplesheet_compl[order(samplesheet_compl$Basename), ]
  
  # background_corr_RGset:
  sample_names <- order(sampleNames(background_corr_RGset))
  background_corr_RGset <- background_corr_RGset[, sample_names]
  
  # detection p values 
  if (is.matrix(dp)){
    dp <- dp[, colnames(dp) %in% overlap_samples]
    dp <- dp[, sort(colnames(dp))]
  }
  
  ################################################################################
  # 2) EXTRACT FROM RGset
  ################################################################################
  
  print("Step 2 - Extract from RGset: snp probes, ctrl probes, bead counts, detection p values (this step can take quite long, depending on the size of the data set)")
  
  # GET BETA VALUES OF SNP PROBES AND ACCORDING SAMPLE NAME
  snp_betas <- getSnpBeta(background_corr_RGset)
  
  # CONTROL PROBES
  controls=getProbeInfo(background_corr_RGset, type = "Control")
  types=unique(controls$Type)
  types=types[types!='NEGATIVE']
  ctrl.names=controls[controls$Type %in% types,'ExtendedType']
  ctrl.address=controls[controls$Type %in% types,'Address']
  ctrl.Green <- matrix(NA_real_, ncol = ncol(background_corr_RGset), nrow = length(ctrl.names), dimnames = list(ctrl.names, sampleNames(background_corr_RGset)))
  ctrl.Green[ctrl.names,] <- getGreen(background_corr_RGset)[ctrl.address,]
  ctrl.Red <- matrix(NA_real_, ncol = ncol(background_corr_RGset), nrow = length(ctrl.names), dimnames = list(ctrl.names, sampleNames(background_corr_RGset)))
  ctrl.Red[ctrl.names,] <- getRed(background_corr_RGset)[ctrl.address,]
  rownames(ctrl.Green)=paste(rownames(ctrl.Green), '.grn', sep='')
  rownames(ctrl.Red)=paste(rownames(ctrl.Red), '.red', sep='')
  ctrl = rbind(ctrl.Green, ctrl.Red)
  rm(ctrl.names, ctrl.address, ctrl.Green, ctrl.Red, types, controls)
  
  # GET BEAD COUNTS
  bc <- beadcount(background_corr_RGset)
  bc_probe <- rowSums(is.na(bc))/ncol(bc) #NAs represent probes with beadcount < 3
  CpGs_lowbeadcount = names(bc_probe[bc_probe > 0.05])
  rm(bc, bc_probe)
  
  # DETECTION P VALUES
  if (is.matrix(dp)){
    dp <- dp[, colnames(dp) %in% samplesheet_compl$Basename]
  } else {
    dp <- minfi::detectionP(background_corr_RGset)
    dp <- dp[, colnames(dp) %in% samplesheet_compl$Basename]
  }
  
  if (length(setdiff(colnames(dp), samplesheet_compl$Basename)) > 0){
    stop("colnames(dp) and samplesheet$Basename do not overlap.")
  }
  if (length(setdiff(samplesheet_compl$Basename, colnames(dp))) > 0){
    stop("colnames(dp) and samplesheet$Basename do not overlap.")
  }
  
  
  
  ################################################################################
  # 3) EXTRACT INTENSITIES FROM RGset
  ################################################################################
  
  print("Step 3 - Extract intensities from RGset")
  
  # TYPE II PROBES
  TypeII.Name <- getProbeInfo(background_corr_RGset, type = "II")$Name
  TypeII.Green <- getGreen(background_corr_RGset)[getProbeInfo(background_corr_RGset, type = "II")$AddressA,]
  TypeII.Red <- getRed(background_corr_RGset)[getProbeInfo(background_corr_RGset, type = "II")$AddressA,]
  rownames(TypeII.Red) <- TypeII.Name
  colnames(TypeII.Red) <- sampleNames(background_corr_RGset)
  rownames(TypeII.Green) <- TypeII.Name
  colnames(TypeII.Green) <- sampleNames(background_corr_RGset)
  
  # TYPE I PROBES, SPLIT INTO GREEN AND RED CHANNELS
  # green
  TypeI.Green.Name <- getProbeInfo(background_corr_RGset, type = "I-Green")$Name
  TypeI.Green.M <- getGreen(background_corr_RGset)[getProbeInfo(background_corr_RGset, type = "I-Green")$AddressB,]
  rownames(TypeI.Green.M) <- TypeI.Green.Name
  colnames(TypeI.Green.M) <- sampleNames(background_corr_RGset)
  TypeI.Green.U <- getGreen(background_corr_RGset)[getProbeInfo(background_corr_RGset, type = "I-Green")$AddressA,]
  rownames(TypeI.Green.U) <- TypeI.Green.Name
  colnames(TypeI.Green.U) <- sampleNames(background_corr_RGset)
  
  # red
  TypeI.Red.Name <- getProbeInfo(background_corr_RGset, type = "I-Red")$Name
  TypeI.Red.M <- getRed(background_corr_RGset)[getProbeInfo(background_corr_RGset, type = "I-Red")$AddressB,]
  rownames(TypeI.Red.M) <- TypeI.Red.Name
  colnames(TypeI.Red.M) <- sampleNames(background_corr_RGset)
  TypeI.Red.U <- getRed(background_corr_RGset)[getProbeInfo(background_corr_RGset, type = "I-Red")$AddressA,]
  rownames(TypeI.Red.U) <- TypeI.Red.Name
  colnames(TypeI.Red.U) <- sampleNames(background_corr_RGset)
  
  samples <- overlap_samples
  
  rm(background_corr_RGset)
  
  ################################################################################
  # 4) APPLY DETECTION P VALUE THRESHOLD
  ################################################################################
  
  print("Step 4 - Apply detection p value threshold")
  
  thres <- 1E-16
  
  # SET VALUES ABOVE DETECTION P VALUE THRESHOLD TO NA IN INTENSITIES DFS
  d=dp[rownames(TypeII.Green),colnames(TypeII.Green)]
  TypeII.Green.d = ifelse(d<thres,TypeII.Green,NA)
  TypeII.Red.d = ifelse(d<thres,TypeII.Red,NA)
  d=dp[rownames(TypeI.Green.M),colnames(TypeI.Green.M)]
  TypeI.Green.M.d = ifelse(d<thres,TypeI.Green.M,NA)
  TypeI.Green.U.d = ifelse(d<thres,TypeI.Green.U,NA)
  d=dp[rownames(TypeI.Red.M),colnames(TypeI.Red.M)]
  TypeI.Red.M.d = ifelse(d<thres,TypeI.Red.M,NA)
  TypeI.Red.U.d = ifelse(d<thres,TypeI.Red.U,NA)
  rm(dp,d)
  
  
  ################################################################################
  # 5) FILTER INTENSITIES ON SNP CPGS, EXCLUDE PROBES WITH LOW BEAD COUNT
  ################################################################################
  
  print("Step 5 - Filter intensities on SNP CPGs, exclude probes with low bead count")
  
  SNP_cgs <- setdiff(SNP_cgs, CpGs_lowbeadcount)
  
  markers=as.matrix(intersect(rownames(TypeII.Green.d), SNP_cgs))
  TypeII.Green = TypeII.Green.d[markers,samples]
  TypeII.Red = TypeII.Red.d[markers,samples]
  markers=intersect(rownames(TypeI.Green.M.d), SNP_cgs)
  TypeI.Green.M = TypeI.Green.M.d[markers,samples]
  TypeI.Green.U = TypeI.Green.U.d[markers,samples]
  markers=intersect(rownames(TypeI.Red.M.d), SNP_cgs)
  TypeI.Red.M = TypeI.Red.M.d[markers,samples]
  TypeI.Red.U = TypeI.Red.U.d[markers,samples]
  
  
  ################################################################################
  # 6) EXCLUDE PROBES WITH LOW CALL RATE 
  ################################################################################
  
  print("Step 6 - Exclude probes with low call rate")
  
  # CALCULATE BETAS FOR DIFFERENT PROBE TYPES & COMBINE
  TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
  TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
  TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
  beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))
  
  marker.call=rowSums(!is.na(beta))/ncol(beta)
  
  # IDENTIFY CPGs WITH A LOW MARKER CALL RATE & EXPORT
  CpGs_low_call_rate <- names(marker.call[marker.call < 0.90])
  
  rm(beta)
  
  # FILTERED PROBES
  SNP_cgs_red <- setdiff(SNP_cgs, CpGs_low_call_rate) #to keep
  
  # FILTER ON PROBES
  markers=intersect(rownames(TypeII.Green.d), SNP_cgs_red)
  TypeII.Green = TypeII.Green.d[markers,samples]
  TypeII.Red = TypeII.Red.d[markers,samples]
  markers=intersect(rownames(TypeI.Green.M.d), SNP_cgs_red)
  TypeI.Green.M = TypeI.Green.M.d[markers,samples]
  TypeI.Green.U = TypeI.Green.U.d[markers,samples]
  markers=intersect(rownames(TypeI.Red.M.d), SNP_cgs_red)
  TypeI.Red.M = TypeI.Red.M.d[markers,samples]
  TypeI.Red.U = TypeI.Red.U.d[markers,samples]
  
  rm(TypeII.Green.d, TypeII.Red.d, TypeI.Green.M.d, TypeI.Green.U.d, TypeI.Red.M.d, TypeI.Red.U.d)
  
  
  ################################################################################
  # 7) QUANTILE NORMALIZATION
  ################################################################################
  
  print("Step 7 - Quantile normalization of intensities")
  
  # NORMALISE
  TypeII.Green = normalizeQuantiles(TypeII.Green)
  TypeII.Red = normalizeQuantiles(TypeII.Red)
  TypeI.Green.M = normalizeQuantiles(TypeI.Green.M)
  TypeI.Green.U = normalizeQuantiles(TypeI.Green.U)
  TypeI.Red.M = normalizeQuantiles(TypeI.Red.M)
  TypeI.Red.U = normalizeQuantiles(TypeI.Red.U)
  
  # CALCULATE BETA VALUES
  TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
  TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
  TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
  beta_SNP0bp = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))
  
  rm(TypeII.Green,TypeII.Red,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U,TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
  
  
  
  ################################################################################
  # 8) IMPUTE MISSING VALUES
  ################################################################################
  
  print("Step 8 - Impute missing values")
  
  imputed_beta <- champ.impute(beta = beta_SNP0bp, pd = samplesheet_compl, method = "KNN")
  beta_SNP0bp_imp <- imputed_beta$beta
  rm(beta_SNP0bp, imputed_beta)
  
  
  
  ################################################################################
  # 9) CALL GENOTYPES OF SNP PROBES
  ################################################################################
  
  print("Step 9 - call genotypes of SNP probes")
  
  # FUNCTION: calculate beta from genotypes (from meffil package)
  calculate.beta.genotypes <- function(snp.betas, centers=c(0.2,0.5,0.8)) {
    x <- t(apply(snp.betas,1,function(x) {
      tryCatch(kmeans(x, centers=centers)$cluster - 1,
               error=function(e) {
                 cluster <- rep(1,ncol(snp.betas))
                 cluster[which(x < min(centers))] <- 0
                 cluster[which(x > max(centers))] <- 2
                 cluster
               })
    }))
    dimnames(x) <- dimnames(snp.betas)
    x
  }
  
  # exclude samples from rs SNPs that were excluded in filtering step
  snp_betas <- snp_betas[, colnames(snp_betas) %in% samplesheet_compl$Basename]
  
  rs_genotypes <- calculate.beta.genotypes(snp_betas)
  rs_genotypes[rs_genotypes == 1] <- 0.5
  rs_genotypes[rs_genotypes == 2] <- 1
  
  rm(snp_betas, calculate.beta.genotypes)
  
  
  ################################################################################
  # 10) CELL TYPE PCs AND CONTROL PROBE PCs
  ################################################################################
  
  print("Step 10 - calculate cell type PCs and control probe PCs")
  
  
  ### EXTRACT FROM samplesheet_compl ###
  Age <- samplesheet_compl$Age
  Sex <- samplesheet_compl$Sex
  
  
  ### CONTROL PROBE PCA ###
  ctrl <- ctrl[, colnames(ctrl) %in% overlap_samples] 
  ctrl_pca <- prcomp(na.omit(t(ctrl))) # run pca
  ctrlprobe_PCAscores <- as.data.frame(ctrl_pca$x) #extract PCA scores
  
  Ctrl_PC1 <- scale(ctrlprobe_PCAscores$PC1)[, 1]
  Ctrl_PC2 <- scale(ctrlprobe_PCAscores$PC2)[, 1]
  Ctrl_PC3 <- scale(ctrlprobe_PCAscores$PC3)[, 1]
  Ctrl_PC4 <- scale(ctrlprobe_PCAscores$PC4)[, 1]
  Ctrl_PC5 <- scale(ctrlprobe_PCAscores$PC5)[, 1]
  Ctrl_PC6 <- scale(ctrlprobe_PCAscores$PC6)[, 1]
  Ctrl_PC7 <- scale(ctrlprobe_PCAscores$PC7)[, 1]
  Ctrl_PC8 <- scale(ctrlprobe_PCAscores$PC8)[, 1]
  Ctrl_PC9 <- scale(ctrlprobe_PCAscores$PC9)[, 1]
  Ctrl_PC10 <- scale(ctrlprobe_PCAscores$PC10)[, 1]
  rm(ctrl, ctrl_pca, ctrlprobe_PCAscores) 
  
  ### CELL TYPE PCA ###
  cell_type_df <- samplesheet_compl[, cell_cols]
  cellcounts_pca <- prcomp(cell_type_df) # run pca
  cellcounts_PCAscores <- as.data.frame(cellcounts_pca$x) #extract PCA scores
  
  cellcounts_PC1 <- scale(cellcounts_PCAscores$PC1)[, 1]
  cellcounts_PC2 <- scale(cellcounts_PCAscores$PC2)[, 1]
  
  nr_cell_type_prop_PCs <- length(cell_cols) - 1
  if (nr_cell_type_prop_PCs >= 3){
    cellcounts_PC3 <- scale(cellcounts_PCAscores$PC3)[, 1]
  }
  if (nr_cell_type_prop_PCs >= 4){
    cellcounts_PC4 <- scale(cellcounts_PCAscores$PC4)[, 1]
  }
  if (nr_cell_type_prop_PCs >= 5){
    cellcounts_PC5 <- scale(cellcounts_PCAscores$PC5)[, 1]
  }
  
  rm(cellcounts_PCAscores, cellcounts_pca, cell_type_df)
  
  
  ################################################################################
  # 11) RESIDUALIZE BETA SNP0BP AND COMBINE WITH RS BETA
  ################################################################################
  
  print("Step 11 - residualise beta SNP0bp and combine with rs beta")
  
  ### RESIDUALISED BETA SNP0bp ###
  beta_SNP0bp_imp <- beta_SNP0bp_imp[, colnames(beta_SNP0bp_imp) %in% samplesheet_compl$Basename]
  resid_beta <- as.data.frame(matrix(data = NA, nrow = nrow(beta_SNP0bp_imp), ncol = ncol(beta_SNP0bp_imp)))
  
  for (i in c(1:nrow(beta_SNP0bp_imp))){
    
    if (nr_cell_type_prop_PCs == 2){
      model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + 
                    Sex + Age + 
                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
    }
    
    if (nr_cell_type_prop_PCs == 3){
      model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 +
                    Sex + Age + 
                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
    }
    
    if (nr_cell_type_prop_PCs == 4){
      model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 + cellcounts_PC4 + 
                    Sex + Age + 
                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
    }
    
    if (nr_cell_type_prop_PCs >= 5){
      model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 + cellcounts_PC4 + cellcounts_PC5 + 
                    Sex + Age + 
                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
    }
    
    
    
    
    residuals <- resid(model)
    resid_beta[i, ] <- residuals
  }
  rownames(resid_beta) <- rownames(beta_SNP0bp_imp)
  colnames(resid_beta) <- colnames(beta_SNP0bp_imp)
  
  
  save(resid_beta, file = paste0(output_folder, "resid_beta_Newmethod.RData"))
  save(beta_SNP0bp_imp, file = paste0(output_folder, "impSNPs_Oldmethod.RData"))
  
  
  rm(model)
  
  
  ### COMBINE RESIDUALISED SNP0bp with rs beta ###
  
  rs_genotypes <- rs_genotypes[, colnames(rs_genotypes) %in% colnames(resid_beta)]
  comb_SNPs <- rbind(resid_beta, rs_genotypes)
  
  save(comb_SNPs, file = paste0(output_folder, "comb_SNPs_Newmethod.RData"))
  
  rm(resid_beta, beta_SNP0bp_imp)
  
}



##################

#' Pipeline to test the performance of the ancestry methods
#'
#' @param SNP_list list of data frames SNPs/CpGs as rows and sample names as columns
#' @param SNP_list_names names of the data frames in SNP_list (names of ancestry information to test)
#' @param samplesheet_compl samplesheet with columns Basename, anc_classif (and Indiv_ID if several samples per individual). 
#' anc_classif should not have missing values
#'
#' @export 
#'
performance_pipeline <- function(SNP_list, SNP_list_names, samplesheet_compl, output_folder){
  
  # 1) DATA FRAMES TO ADD RESULTS
  print ("Step 1) Set up data frames to save results.")
  
    centroids_df_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 9))
    colnames(centroids_df_comb) <- c("anc_classif", "centroid_x", "centroid_y", "centroid_z", "mean_within_dist_3D", "mean_betw_dist_3D", "mean_within_dist_2D", "mean_betw_dist_2D", "method")
    
    dist_matrix_2D_long_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
    colnames(dist_matrix_2D_long_comb) <- c("Ancestry1", "Ancestry2", "Distance", "method")
    
    dist_matrix_3D_long_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
    colnames(dist_matrix_3D_long_comb) <- c("Ancestry1", "Ancestry2", "Distance", "method")
    
    
    sample_dist_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
    colnames(sample_dist_comb) <- c("Basename", "Indiv_ID", "anc_classif", "dist_3D", "dist_2D", "method")
    
    statistics_results_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
    colnames(statistics_results_comb) <- c("PC", "test", "p_value", "effect_size", "method")
  
    sil_df_3D_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
    colnames(sil_df_3D_comb) <- c("cluster", "neighbor", "sil_width", "label", "method")
    
    sil_df_2D_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
    colnames(sil_df_2D_comb) <- c("cluster", "neighbor", "sil_width", "label", "method")
    
    mean_sil_df_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
    colnames(mean_sil_df_comb) <- c("dimension", "method", "mean_sil")
    
    mean_silhouette_3D_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
    colnames(mean_silhouette_3D_comb) <- c("label", "mean_silhouette", "method")
    
    mean_silhouette_2D_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
    colnames(mean_silhouette_3D_comb) <- c("label", "mean_silhouette", "method")
  
    
    
    if (multiple_samples_per_indiv == "yes"){
      mean_within_dist_df_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
      colnames(mean_within_dist_df_comb) <- c("dimension", "mean_value", "method")
      
      mean_sil_df_within_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
      colnames(mean_sil_df_within_comb) <- c("dimension", "method", "mean_sil")
      
      indiv_centroids_df_comb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
      colnames(indiv_centroids_df_comb) <- c("Indiv_ID", "centroids_x", "centroid_y", "centroid_z", "mean_within_dist_3D", "mean_within_dist_2D",  "method")
    }
    
    
  # 2) GET PERFORMANCE RESULTS 
  print("Step 2) Get performance results for each method.")
  
  for (i in c(1:length(SNP_list_names))){
    
    print(paste0("running ", SNP_list_names[i], " method"))
    
    
    
    
    # A) RUN PCA
    print("A) Run PCA.")
    
      # check overlap between samples in samplesheet_compl and SNP_list object
      overlapping_samples <- intersect(colnames(SNP_list[[i]]), samplesheet_compl$Basename)
      
      if (length(overlapping_samples) < 5){
        stop("The samplenames do not overlap: colnames(SNP_list[[i]]) and samplesheet_compl$Basename.")
      }
      
      if(length(overlapping_samples) != length(colnames(SNP_list[[i]]))){
        print("No exact overlap between the samples of the samplesheet and the SNP_info object. Only overlapping samples are used.")
      }
      
      # select overlapping samples
      SNP_list[[i]] <- SNP_list[[i]][, overlapping_samples]
      samplesheet_compl <- samplesheet_compl[which(samplesheet_compl$Basename %in% overlapping_samples), ]
      
      # make sure the sample names are in the same order
      SNP_list[[i]] <- SNP_list[[i]][, order(colnames(SNP_list[[i]]))]
      samplesheet_compl <- samplesheet_compl[order(samplesheet_compl$Basename), ]
      
      # run PCA and extract the first 5 PCs
      pc <- prcomp(t(SNP_list[[i]]))
      top5pc <- as.data.frame(pc$x[,1:5])
      
      #scale the PCs
      top5pc <- as.data.frame(lapply(top5pc, scale_min_max))
      
      # add the Basename
      top5pc$Basename <- samplesheet_compl$Basename
      
      # remove objects that are not needed
      rm(pc)
    

    
      
    # B) SCATTER PLOTS (OPEN IN FIREFOX)  
    
    print("B) Generate scatterplots.")
      
      # generate 3D plot
      scatterplot_3D <-  plot_ly(top5pc, x = ~PC1, y = ~PC2, z = ~PC3, color = samplesheet_compl$anc_classif,
                                 type = "scatter3d", mode = "markers", marker = list(size = 4, opacity = 1))
      scatterplot_3D <- scatterplot_3D %>% layout(scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")), title = paste0("PCA - ", SNP_list_names[i]))
      
      # save 3D plot
      saveWidget(scatterplot_3D, file = paste0(output_folder, "Figures/Scatterplot3D_", SNP_list_names[i], ".html"))
      
      # merge PCA with samplesheet
      comb_df <- merge(top5pc, samplesheet_compl, by = "Basename")
      
      # PC1 vs PC2
      PC12 <- ggplot(comb_df, aes(x= PC1, y = PC2, colour = anc_classif)) + 
        geom_point() + 
        labs(colour = "Genetic Ancestry", x = "ancestry PC1", y = "ancestry PC2") + 
        th + th_transparent + 
        theme(legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_text(face = "bold")) + 
        ggtitle(SNP_list_names[i])
      ggsave(paste0(output_folder, "Figures/PC12_", SNP_list_names[i], ".pdf"), plot = PC12, device = "pdf", width =6, height = 4.5, dpi = 400)
      
      # PC2 vs PC3
      PC23 <- ggplot(comb_df, aes(x= PC2, y = PC3, colour = anc_classif)) + 
        geom_point() + 
        labs(colour = "Genetic Ancestry", x = "ancestry PC2", y = "ancestry PC3") + 
        th + th_transparent + 
        theme(legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_text(face = "bold")) + 
        ggtitle(SNP_list_names[i])
      ggsave(paste0(output_folder, "Figures/PC23_", SNP_list_names[i], ".pdf"), plot = PC23, device = "pdf", width =6, height = 4.5, dpi = 400)
      
      # combine the two plots (PC1/2 and PC2/3)
      legend <- get_legend(PC23)
      legend_plot <- ggdraw(legend)
      PC12 <- PC12 + theme(legend.position = "none", plot.title = element_blank())
      PC23 <- PC23 + theme(legend.position = "none", plot.title = element_blank())
      design <- "AAABBBC"
      layout <- (PC12 + PC23 + legend_plot) + patchwork::plot_layout(design = design) + 
        patchwork::plot_annotation(title = paste0("Ancestry PCs - ", SNP_list_names[i]),
                        theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      ggsave(paste0(output_folder, "Figures/PC123_", SNP_list_names[i], ".pdf"), layout, units = "cm", width = 27, height = 13, dpi = 400)
      
      
      
      
      
    # C) TEST ASSOCIATION BETWEEN PCs AND ANCESTRY
    print("C) Test association between PCs and ancestry.")
      
      # create df for statistics results
      statistics_results <- data.frame(PC = c(1,2,3,4,5), test = NA, p_value = NA, effect_size = NA, method = SNP_list_names[i])
      
      # test for every PC
      for (k in c(1:5)){
        
        # extract PC values and ancestry info
        PC <- paste0("PC", k)
        PC_values <- comb_df[, PC]
        ancestry_values <- factor(comb_df$anc_classif)
        
        # Shapiro test
        shapiro_result <- shapiro.test(residuals(aov(PC_values ~ ancestry_values)))
        shapiro_p <- shapiro_result$p.value
        
        # Levene test
        levene_result <- leveneTest(PC_values ~ ancestry_values)
        levene_p <- levene_result$`Pr(>F)`[1]
        
        # Kruskal-Wallis or ANOVA depending on results of Shapiro and Levene test
        if (shapiro_p < 0.05 | levene_p < 0.05){
          
          statistics_results$test[k] <- "Kruskal-Wallis"
          
          kruskal_test <- kruskal.test(PC_values ~ ancestry_values)
          kruskal_p <- format(kruskal_test$p.value, scientific = TRUE)
          
          statistics_results$p_value[k] <- kruskal_p
          
          kruskal_effect_size <- round(rstatix::kruskal_effsize(data.frame(PC = PC_values, ancestry = ancestry_values), formula = (PC ~ ancestry))$effsize, 5)
          statistics_results$effect_size[k] <- kruskal_effect_size
          
        } else {
          
          statistics_results$test[k] <- "ANOVA"
          
          anova_test <- aov(PC_values ~ ancestry_values)
          anova_sum <- summary(anova_test)
          anova_p <- format(anova_sum[[1]]["Pr(>F)"][1, 1], scientific = TRUE)
          
          statistics_results$p_value[k] <- anova_p
          
          eta_squar <- round(effectsize::eta_squared(anova_test)$Eta2, 5)
          statistics_results$effect_size[k] <- eta_squar
        }
      }
    
      
      
    
    # D) CALCULATE CENTROIDS
    
      print("D) Calculate centroids.")
      
        centroids_df <- as.data.frame(comb_df %>%
                                     group_by(anc_classif) %>%
                                     summarise(
                                       centroid_x = mean(PC1),
                                       centroid_y = mean(PC2),
                                       centroid_z = mean(PC3)
                                     ))
          
        centroids_df$mean_within_dist_3D <- NA
        centroids_df$mean_betw_dist_3D <- NA
        centroids_df$mean_within_dist_2D <- NA
        centroids_df$mean_betw_dist_2D <- NA
        centroids_df$method <-  SNP_list_names[i]
        
    
    
        
      # E) CALCULATE THE DISTANCE FROM EACH SAMPLE TO THE CENTROID
        
        print("E) Calculate the distance from each sample to the centroid.")
        
          # Calculate distance from centroid to sample
          sample_dist <- as.data.frame(matrix(data=NA, nrow = 0, ncol = ncol(comb_df)))
          
          for (k in c(1:nrow(centroids_df))){
            anc <- centroids_df$anc_classif[k]
            anc_df <- comb_df[which(comb_df$anc_classif == anc), ]
            anc_df$dist_3D <- NA
            anc_df$dist_2D <- NA
            for (j in c(1:nrow(anc_df))){
              anc_df$dist_3D[j] <- euclidean_distance(
                c(anc_df$PC1[j], anc_df$PC2[j], anc_df$PC3[j]),
                c(centroids_df$centroid_x[k], centroids_df$centroid_y[k], centroids_df$centroid_z[k])
              )
            }
            
            for (j in c(1:nrow(anc_df))){
              anc_df$dist_2D[j] <- euclidean_distance(
                c(anc_df$PC1[j], anc_df$PC2[j]),
                c(centroids_df$centroid_x[k], centroids_df$centroid_y[k])
              )
            }
            sample_dist <- rbind(sample_dist, anc_df)
          }
          
          sample_dist$method <- SNP_list_names[i]
          
          # sample distance df with fewer columns
          sample_dist_red <- sample_dist[, c("Basename", "Indiv_ID", "anc_classif", "dist_3D", "dist_2D", "method")]
    
          # add mean sample distance to centroids_df
          unique_anc <- unique(centroids_df$anc_classif)
          for (anc in unique_anc){
            centroids_df$mean_within_dist_3D[which(centroids_df$anc_classif == anc)] <- mean(sample_dist$dist_3D[which(sample_dist$anc_classif == anc)])
            centroids_df$mean_within_dist_2D[which(centroids_df$anc_classif == anc)] <- mean(sample_dist$dist_2D[which(sample_dist$anc_classif == anc)])
          }
          
    
          
          
    # F) CREATE DISTANCE MATRIX BETWEEN ANCESTRIES
    
    print("F) Create distance matrix between ancestries.")
        
      # Distance matrix 3D
      dist_matrix_3D <- as.data.frame(matrix(data = NA, nrow = nrow(centroids_df), ncol = nrow(centroids_df)))
      rownames(dist_matrix_3D) <- centroids_df$anc_classif
      colnames(dist_matrix_3D) <- centroids_df$anc_classif
      
      for (k in c(1:nrow(centroids_df))){
        for (j in c(1:nrow(centroids_df))){
          dist_matrix_3D[k, j] <- euclidean_distance(
            c(centroids_df$centroid_x[k], centroids_df$centroid_y[k], centroids_df$centroid_z[k]),
            c(centroids_df$centroid_x[j], centroids_df$centroid_y[j], centroids_df$centroid_z[j])
          )
        }
      }
      
      # Distance matrix 2D
      dist_matrix_2D <- as.data.frame(matrix(data = NA, nrow = nrow(centroids_df), ncol = nrow(centroids_df)))
      rownames(dist_matrix_2D) <- centroids_df$anc_classif
      colnames(dist_matrix_2D) <- centroids_df$anc_classif
      
      for (k in c(1:nrow(centroids_df))){
        for (j in c(1:nrow(centroids_df))){
          dist_matrix_2D[k, j] <- euclidean_distance(
            c(centroids_df$centroid_x[k], centroids_df$centroid_y[k]),
            c(centroids_df$centroid_x[j], centroids_df$centroid_y[j])
          )
        }
      }
    
    
      # Add mean distance between clusters to centroids_df
      for (k in c(1:nrow(dist_matrix_3D))){
        centroids_df$mean_betw_dist_3D[k] <- mean(as.numeric(dist_matrix_3D[k, ]))
      }
      for (k in c(1:nrow(dist_matrix_2D))){
        centroids_df$mean_betw_dist_2D[k] <- mean(as.numeric(dist_matrix_2D[k, ]))
      }
      
      
      # long format of distance matrix - 3D
      dist_matrix_3D_long <- as.data.frame(dist_matrix_3D %>%
                                          rownames_to_column(var = "Ancestry1") %>%
                                          pivot_longer(cols = -Ancestry1, names_to = "Ancestry2", values_to = "Distance"))
      dist_matrix_3D_long$method <- SNP_list_names[i]
      dist_matrix_3D_long <- dist_matrix_3D_long[which(dist_matrix_3D_long$Distance != 0), ]
        
      # long format of distance matrix - 2D
      dist_matrix_2D_long <- as.data.frame(dist_matrix_2D %>%
                                             rownames_to_column(var = "Ancestry1") %>%
                                             pivot_longer(cols = -Ancestry1, names_to = "Ancestry2", values_to = "Distance"))
      dist_matrix_2D_long$method <- SNP_list_names[i]
      dist_matrix_2D_long <- dist_matrix_2D_long[which(dist_matrix_2D_long$Distance != 0), ]
      
    
    # G) WITHIN INDIVIDUALS
    
    print("G) distance within individuals (if multiple samples per individual)")
    
    if (multiple_samples_per_indiv == "yes"){
      
      # select individuals that have more than 1 sample
      duplic_entries <- data.frame(table(comb_df$Indiv_ID))
      duplic_entries <- duplic_entries$Var1[which(duplic_entries$Freq > 1)]
      comb_df2 <- comb_df[which(comb_df$Indiv_ID %in% duplic_entries), ]
      
      # calculate centroids of individuals
      indiv_centroids_df <- as.data.frame(comb_df2 %>%
                                            group_by(Indiv_ID) %>%
                                            summarise(
                                              centroid_x = mean(PC1),
                                              centroid_y = mean(PC2),
                                              centroid_z = mean(PC3)
                                            ))
      indiv_centroids_df$mean_within_dist_3D <- NA
      indiv_centroids_df$mean_within_dist_2D <- NA
      indiv_centroids_df$method <- SNP_list_names[i]
      
      # calculate distance of samples to centroid per individual
      for (k in c(1:nrow(indiv_centroids_df))){
        individual <- indiv_centroids_df$Indiv_ID[k]
        indiv_df <- comb_df2[which(comb_df2$Indiv_ID == individual), ]
        indiv_df$dist_3D <- NA
        indiv_df$dist_2D <- NA
        
        for (j in c(1:nrow(indiv_df))){
          indiv_df$dist_3D[j] <- euclidean_distance(
            c(indiv_df$PC1[j], indiv_df$PC2[j], indiv_df$PC3[j]),
            c(indiv_centroids_df$centroid_x[k], indiv_centroids_df$centroid_y[k], indiv_centroids_df$centroid_z[k])
          )
          indiv_centroids_df$mean_within_dist_3D[k] <- mean(indiv_df$dist_3D)
        }
        for (j in c(1:nrow(indiv_df))){
          indiv_df$dist_2D[j] <- euclidean_distance(
            c(indiv_df$PC1[j], indiv_df$PC2[j]),
            c(indiv_centroids_df$centroid_x[k], indiv_centroids_df$centroid_y[k])
          )
          indiv_centroids_df$mean_within_dist_2D[k] <- mean(indiv_df$dist_2D)
        }
      }
      
      # combine results
      mean_within_dist_df <- data.frame(dimension = c("2D", "3D"), 
                                        mean_value = c(mean(indiv_centroids_df$mean_within_dist_2D), 
                                                       mean(indiv_centroids_df$mean_within_dist_3D)), 
                                        method = SNP_list_names[i])
      
      
      
    }
    
    
    # H) SILHOUETTE SCORES
    
    print("H) Silhouette scores")
    
    
      # silhouette scores 3D
      pca_data <- data.frame(PC1 = comb_df$PC1, PC2 = comb_df$PC2, PC3 = comb_df$PC3)
      labels <- as.numeric(factor(comb_df$anc_classif))
      sil <- silhouette(labels, dist(pca_data))
      sil_df_3D <- as.data.frame(sil)
      sil_df_3D$label <- factor(comb_df$anc_classif)
      sil_df_3D$method <- SNP_list_names[i]
      
      # silhouette scores 2D
      pca_data <- data.frame(PC1 = comb_df$PC1, PC2 = comb_df$PC2)
      labels <- as.numeric(factor(comb_df$anc_classif))
      sil <- silhouette(labels, dist(pca_data))
      sil_df_2D <- as.data.frame(sil)
      sil_df_2D$label <- factor(comb_df$anc_classif)
      sil_df_2D$method <- SNP_list_names[i]
    
      # mean silhouette scores
      mean_sil_df <- data.frame(dimension = c("3D", "2D"), method = SNP_list_names[i], 
                                mean_sil = c(mean(sil_df_3D$sil_width, na.rm = TRUE), mean(sil_df_2D$sil_width, na.rm = TRUE)))
    
      # mean silhouette scores per ancestry group in 3D
      mean_silhouette_3D <- sil_df_3D %>%
        group_by(label) %>%
        summarise(mean_silhouette = mean(sil_width))
      mean_silhouette_3D <- as.data.frame(mean_silhouette_3D)
      mean_silhouette_3D$method <- SNP_list_names[i]
   
      # mean silhouette scores per ancestry group in 2D
      mean_silhouette_2D <- sil_df_2D %>%
        group_by(label) %>%
        summarise(mean_silhouette = mean(sil_width))
      mean_silhouette_2D <- as.data.frame(mean_silhouette_2D)
      mean_silhouette_2D$method <- SNP_list_names[i]
   
      
    # I) SILHOUETTE SCORES
      
      print("I) Silhouette scores within individuals (if multiple samples per individual)")
      
      if (multiple_samples_per_indiv == "yes"){
        # silhouette scores 3D
        pca_data <- data.frame(PC1 = comb_df$PC1, PC2 = comb_df$PC2, PC3 = comb_df$PC3)
        labels <- as.numeric(factor(comb_df$Indiv_ID))
        sil <- silhouette(labels, dist(pca_data))
        sil_df_3D_within <- as.data.frame(sil)
        sil_df_3D_within$label <- factor(comb_df$Indiv_ID)
        sil_df_3D_within$method <- SNP_list_names[i]
        
        # silhouette scores 2D
        pca_data <- data.frame(PC1 = comb_df$PC1, PC2 = comb_df$PC2)
        labels <- as.numeric(factor(comb_df$Indiv_ID))
        sil <- silhouette(labels, dist(pca_data))
        sil_df_2D_within <- as.data.frame(sil)
        sil_df_2D_within$label <- factor(comb_df$Indiv_ID)
        sil_df_2D_within$method <- SNP_list_names[i]
        
        # mean silhouette scores
        mean_sil_df_within <- data.frame(dimension = c("3D", "2D"), method = SNP_list_names[i], 
                                         mean_sil = c(mean(sil_df_3D_within$sil_width, na.rm = TRUE), mean(sil_df_2D_within$sil_width, na.rm = TRUE)))
      }
     
   
      
    # J) COMBINE AND SAVE RESULTS
      
      print("J) Combine results of the different methods.")
    
      centroids_df_comb <- rbind(centroids_df_comb, centroids_df)
      dist_matrix_2D_long_comb <- rbind(dist_matrix_2D_long_comb, dist_matrix_2D_long)
      dist_matrix_3D_long_comb <- rbind(dist_matrix_3D_long_comb, dist_matrix_3D_long)
      sample_dist_comb <- rbind(sample_dist_comb, sample_dist_red)
      statistics_results_comb <- rbind(statistics_results_comb, statistics_results)
      sil_df_3D_comb <- rbind(sil_df_3D_comb, sil_df_3D)
      sil_df_2D_comb <- rbind(sil_df_2D_comb, sil_df_2D)
      mean_sil_df_comb <- rbind(mean_sil_df_comb, mean_sil_df)
      mean_silhouette_3D_comb <- rbind(mean_silhouette_3D_comb, mean_silhouette_3D)
      mean_silhouette_2D_comb <- rbind(mean_silhouette_2D_comb, mean_silhouette_2D)
      if (multiple_samples_per_indiv == "yes"){
        mean_within_dist_df_comb <- rbind(mean_within_dist_df_comb, mean_within_dist_df)
        mean_sil_df_within_comb <- rbind(mean_sil_df_within_comb, mean_sil_df_within)
        indiv_centroids_df_comb <- rbind(indiv_centroids_df_comb, indiv_centroids_df)
      }
      
  }
      # DEFINE method LEVELS FOR RESULTS FROM QUANTIFICATION
      centroids_df_comb$method <- factor(centroids_df_comb$method, levels = SNP_list_names)
      dist_matrix_2D_long_comb$method <- factor(dist_matrix_2D_long_comb$method, levels = SNP_list_names)
      dist_matrix_3D_long_comb$method <- factor(dist_matrix_3D_long_comb$method, levels = SNP_list_names)
      sample_dist_comb$method <- factor(sample_dist_comb$method, levels = SNP_list_names)
      statistics_results_comb$method <- factor(statistics_results_comb$method, levels = SNP_list_names)
      sil_df_3D_comb$method <- factor(sil_df_3D_comb$method, levels = SNP_list_names)
      sil_df_2D_comb$method <- factor(sil_df_2D_comb$method, levels = SNP_list_names)
      mean_sil_df_comb$method <- factor(mean_sil_df_comb$method, levels = SNP_list_names)
      mean_silhouette_3D_comb$method <- factor(mean_silhouette_3D_comb$method, levels = SNP_list_names)
      mean_silhouette_2D_comb$method <- factor(mean_silhouette_2D_comb$method, levels = SNP_list_names)
      if (multiple_samples_per_indiv == "yes"){
        mean_within_dist_df_comb$method <- factor(mean_within_dist_df_comb$method, levels = SNP_list_names)
        mean_sil_df_within_comb$method <- factor(mean_sil_df_within_comb$method, levels = SNP_list_names)
        indiv_centroids_df_comb$method <- factor(indiv_centroids_df_comb$method, levels = SNP_list_names)
      }
      
      
      
      # save results
      write_xlsx(centroids_df_comb, path = paste0(output_folder, "Tables/Centroids_df.xlsx"))
      save(centroids_df_comb, file = paste0(output_folder, "RData/Centroids_df.RData"))
      
      write_xlsx(dist_matrix_2D_long_comb, path = paste0(output_folder, "Tables/Distance_matrix_2D.xlsx"))
      save(dist_matrix_2D_long_comb, file = paste0(output_folder, "RData/Distance_matrix_2D.RData"))
      
      write_xlsx(dist_matrix_3D_long_comb, path = paste0(output_folder, "Tables/Distance_matrix_3D.xlsx"))
      save(dist_matrix_3D_long_comb, file = paste0(output_folder, "RData/Distance_matrix_3D.RData"))
      
      save(sample_dist_comb, file = paste0(output_folder, "RData/Sample_distance.RData"))
      
      write_xlsx(statistics_results_comb, path = paste0(output_folder, "Tables/Statistics.xlsx"))
      save(statistics_results_comb, file = paste0(output_folder, "RData/Statistics.RData"))

      save(sil_df_3D_comb, file = paste0(output_folder, "RData/Silhouette_scores_3D.RData"))
      
      save(sil_df_2D_comb, file = paste0(output_folder, "RData/Silhouette_scores_2D.RData"))
      
      write_xlsx(mean_sil_df_comb, path = paste0(output_folder, "Tables/Mean_silhouette_scores.xlsx"))
      save(mean_sil_df_comb, file = paste0(output_folder, "RData/Mean_silhouette_scores.RData"))
      
      write_xlsx(mean_silhouette_3D_comb, path = paste0(output_folder, "Tables/Mean_silhouette_scores_ancestry_3D.xlsx"))
      save(mean_silhouette_3D_comb, file = paste0(output_folder, "RData/Mean_silhouette_scores_ancestry_3D.RData"))
      
      write_xlsx(mean_silhouette_2D_comb, path = paste0(output_folder, "Tables/Mean_silhouette_scores_ancestry_2D.xlsx"))
      save(mean_silhouette_2D_comb, file = paste0(output_folder, "RData/Mean_silhouette_scores_ancestry_2D.RData"))
      
      if (multiple_samples_per_indiv == "yes"){
        
        write_xlsx(mean_sil_df_within_comb, path = paste0(output_folder, "Tables/Mean_silhouette_scores_within.xlsx"))
        save(mean_sil_df_within_comb, file = paste0(output_folder, "RData/Mean_silhouette_scores_within.RData"))
        
        write_xlsx(mean_within_dist_df_comb, path = paste0(output_folder, "Tables/Mean_within_dist.xlsx"))
        save(mean_within_dist_df_comb, file = paste0(output_folder, "RData/Mean_within_dist.RData"))
        
        write_xlsx(indiv_centroids_df_comb, path = paste0(output_folder, "Tables/Indiv_centroids_df.xlsx"))
        save(indiv_centroids_df_comb, file = paste0(output_folder, "RData/Indiv_centroids_df.RData"))
      }
      
      
      # K) PLOTS
      # One thing that might need to be adapted here is the width if you add more SNP datasets you want to test
      
      print("K) Plots.")
      
      # WITHIN INDIVIDUALS
      if (multiple_samples_per_indiv == "yes"){
        
        #3D
        within_plot_3D <- ggplot(indiv_centroids_df_comb, aes(x = method, y = mean_within_dist_3D, fill = method)) +
          geom_violin(scale = "width", alpha = 0.3) + 
          geom_boxplot(width = 0.1, fill = "white") +
          th + th_transparent +
          ggtitle("Within Individual Distance - 3D") +
          labs(x = "method", y = "average distance to centroid per indivual") + 
          add_costum_colours(levels(indiv_centroids_df_comb$method), "fill") +
          theme(legend.position = "none") +
          theme(plot.title = element_text(hjust = 0.5),
                legend.title = element_blank())
        ggsave(paste0(output_folder, "Figures/Within_indiv_dist_3D.pdf"), within_plot_3D, units = "cm", width = 11, height = 12, dpi = 400)
        
        #2D
        within_plot_2D <- ggplot(indiv_centroids_df_comb, aes(x = method, y = mean_within_dist_2D, fill = method)) +
          geom_violin(scale = "width", alpha = 0.3) + 
          geom_boxplot(width = 0.1, fill = "white") +
          th + th_transparent +
          ggtitle("Within Individual Distance - 2D") +
          labs(x = "method", y = "average distance to centroid per indivual") + 
          add_costum_colours(levels(indiv_centroids_df_comb$method), "fill") +
          theme(legend.position = "none") +
          theme(plot.title = element_text(hjust = 0.5),
                legend.title = element_blank())
        ggsave(paste0(output_folder, "Figures/Within_indiv_dist_2D.pdf"), within_plot_2D, units = "cm", width = 11, height = 12, dpi = 400)
        
      }
      
    # BETWEEN CLUSTER
      
      #3D 
      df_means_between_3D <- as.data.frame(dist_matrix_3D_long_comb %>%
                                             group_by(Ancestry1, method) %>%
                                             summarize(Mean = mean(Distance, na.rm = TRUE),
                                                       se = sd(Distance, na.rm = TRUE) / sqrt(n()),
                                                       sd = sd(Distance)))
      
      between_cluster_3D <- ggplot(df_means_between_3D, aes(x=Ancestry1, y = Mean, colour = method, group = method)) +
        geom_point(size = 3, position = position_dodge(width = 0.5)) + 
        geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width = 0.2, position = position_dodge(width = 0.5)) + 
        geom_point(data = dist_matrix_3D_long_comb, aes(y = Distance, x = Ancestry1, colour = method), position = position_dodge(width = 0.5), alpha = 0.5) +
        th + th_transparent + 
        ggtitle("Between Cluster Distance - 3D") +
        labs(x = "genetic ancestry", y = "distance to centroids of other clusters", fill = "method") +
        add_costum_colours(levels(df_means_between_3D$method), "group") +
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Between_cluster_dist_3D.pdf"), between_cluster_3D, units = "cm", width = 12, height = 12, dpi = 400)
      
      #2D 
      df_means_between_2D <- as.data.frame(dist_matrix_2D_long_comb %>%
                                             group_by(Ancestry1, method) %>%
                                             summarize(Mean = mean(Distance, na.rm = TRUE),
                                                       se = sd(Distance, na.rm = TRUE) / sqrt(n()),
                                                       sd = sd(Distance)))
      
      between_cluster_2D <- ggplot(df_means_between_2D, aes(x=Ancestry1, y = Mean, colour = method, group = method)) +
        geom_point(size = 3, position = position_dodge(width = 0.5)) + 
        geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width = 0.2, position = position_dodge(width = 0.5)) + 
        geom_point(data = dist_matrix_2D_long_comb, aes(y = Distance, x = Ancestry1, colour = method), position = position_dodge(width = 0.5), alpha = 0.5) +
        th + th_transparent + 
        ggtitle("Between Cluster Distance - 2D") +
        labs(x = "genetic ancestry", y = "distance to centroids of other clusters", fill = "method") +
        add_costum_colours(levels(df_means_between_2D$method), "group") +
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Between_cluster_dist_2D.pdf"), between_cluster_2D, units = "cm", width = 12, height = 12, dpi = 400)
      
      
    # WITHIN CLUSTER
      
      # 3D
      df_means_within_3D <- as.data.frame(sample_dist_comb %>%
                                            group_by(anc_classif, method) %>%
                                            summarize(Mean = mean(dist_3D, na.rm = TRUE),
                                                      se = sd(dist_3D, na.rm = TRUE) / sqrt(n()),
                                                      sd = sd(dist_3D)))
      
      
      within_cluster_3D <- ggplot(sample_dist_comb, aes(x=anc_classif, y = dist_3D, fill = method)) +
        geom_violin(scale = "width", width = 0.4, position = position_dodge(width = 0.5), alpha = 0.25, trim = TRUE) +
        geom_point(data = df_means_within_3D, aes(y = Mean, x = anc_classif, fill = method), size = 2, position = position_dodge(width = 0.5)) +
        geom_errorbar(data = df_means_within_3D, aes(y = Mean, ymin = Mean - sd, ymax = Mean + sd), position = position_dodge(width = 0.5), width = 0.2) + 
        th + th_transparent + 
        ggtitle("Within Cluster Distance - 3D") +
        labs(x = "genetic ancestry", y = "distance to centroid", fill = "method") +
        scale_colour_manual(values = c("old" = "black", "improved" = "black")) +
        add_costum_colours(levels(sample_dist_comb$method), "fill") + 
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Within_cluster_dist_3D.pdf"), within_cluster_3D, units = "cm", width = 14, height = 12, dpi = 400)
      
      
      # 2D
      df_means_within_2D <- as.data.frame(sample_dist_comb %>%
                                            group_by(anc_classif, method) %>%
                                            summarize(Mean = mean(dist_2D, na.rm = TRUE),
                                                      se = sd(dist_2D, na.rm = TRUE) / sqrt(n()),
                                                      sd = sd(dist_2D)))
      
      
      within_cluster_2D <- ggplot(sample_dist_comb, aes(x=anc_classif, y = dist_2D, fill = method)) +
        geom_violin(scale = "width", width = 0.4, position = position_dodge(width = 0.5), alpha = 0.25, trim = TRUE) +
        geom_point(data = df_means_within_2D, aes(y = Mean, x = anc_classif, fill = method), size = 2, position = position_dodge(width = 0.5)) +
        geom_errorbar(data = df_means_within_2D, aes(y = Mean, ymin = Mean - sd, ymax = Mean + sd), position = position_dodge(width = 0.5), width = 0.2) + 
        th + th_transparent + 
        ggtitle("Within Cluster Distance - 2D") +
        labs(x = "genetic ancestry", y = "distance to centroid", fill = "method") +
        scale_colour_manual(values = c("old" = "black", "improved" = "black")) +
        add_costum_colours(levels(sample_dist_comb$method), "fill") + 
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Within_cluster_dist_2D.pdf"), within_cluster_2D, units = "cm", width = 14, height = 12, dpi = 400)
      

    # ASSOCIATION
      statistics_results_comb$p_value_log <- -log10(as.numeric(statistics_results_comb$p_value))
      statistics_results_comb_red <- statistics_results_comb[which(statistics_results_comb$PC == 1 |
                                                                     statistics_results_comb$PC == 2 |
                                                                     statistics_results_comb$PC == 3), ]
      
      statistics_plot <- ggplot(statistics_results_comb_red, aes(x = PC, y = p_value_log, fill = method)) + 
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.5, width = 0.8, color = "black") +
        labs(title = "Ancestry PCs ~ Ancestry Categories", 
             x = "ancestry PC",
             y = "-log10(p)") + 
        th + th_transparent + 
        scale_y_continuous(expand = c(0,0)) + 
        add_costum_colours(levels(statistics_results_comb_red$method), "fill") +
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Statistics_plot.pdf"), statistics_plot, units = "cm", width = 10.2, height = 10, dpi = 400)
      
      
      # SILHOUETTE SCORE
      
      # 3D
      df_means_silhouette_3D <- as.data.frame(sil_df_3D_comb %>%
                                            group_by(label, method) %>%
                                            summarize(Mean = mean(sil_width, na.rm = TRUE),
                                                      se = sd(sil_width, na.rm = TRUE) / sqrt(n()),
                                                      sd = sd(sil_width)))
      
      
      silhouette_3D <- ggplot(sil_df_3D_comb, aes(x=label, y = sil_width, fill = method)) +
        geom_violin(scale = "width", width = 0.4, position = position_dodge(width = 0.5), alpha = 0.25, trim = TRUE) +
        geom_point(data = df_means_silhouette_3D, aes(y = Mean, x = label, fill = method), size = 2, position = position_dodge(width = 0.5)) +
        geom_errorbar(data = df_means_silhouette_3D, aes(y = Mean, ymin = Mean - sd, ymax = Mean + sd), position = position_dodge(width = 0.5), width = 0.2) + 
        th + th_transparent + 
        ggtitle("Silhouette Scores - 3D") +
        labs(x = "genetic ancestry", y = "silhouette scores", fill = "method") +
        scale_colour_manual(values = c("old" = "black", "improved" = "black")) +
        add_costum_colours(levels(sil_df_3D_comb$method), "fill") +
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Silhouette_3D.pdf"), silhouette_3D, units = "cm", width = 16, height = 12, dpi = 400)
      
      # 2D
      df_means_silhouette_2D <- as.data.frame(sil_df_2D_comb %>%
                                                group_by(label, method) %>%
                                                summarize(Mean = mean(sil_width, na.rm = TRUE),
                                                          se = sd(sil_width, na.rm = TRUE) / sqrt(n()),
                                                          sd = sd(sil_width)))
      
      
      silhouette_2D <- ggplot(sil_df_2D_comb, aes(x=label, y = sil_width, fill = method)) +
        geom_violin(scale = "width", width = 0.4, position = position_dodge(width = 0.5), alpha = 0.25, trim = TRUE) +
        geom_point(data = df_means_silhouette_2D, aes(y = Mean, x = label, fill = method), size = 2, position = position_dodge(width = 0.5)) +
        geom_errorbar(data = df_means_silhouette_2D, aes(y = Mean, ymin = Mean - sd, ymax = Mean + sd), position = position_dodge(width = 0.5), width = 0.2) + 
        th + th_transparent + 
        ggtitle("Silhouette Scores - 2D") +
        labs(x = "genetic ancestry", y = "silhouette scores", fill = "method") +
        scale_colour_manual(values = c("old" = "black", "improved" = "black")) +
        add_costum_colours(levels(sil_df_2D_comb$method), "fill") + 
        theme(legend.position = "top",
              legend.box.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
      
      ggsave(paste0(output_folder, "Figures/Silhouette_2D.pdf"), silhouette_2D, units = "cm", width = 16, height = 12, dpi = 400)
      
      
      #WITHIN SILHOUETTE
      if (multiple_samples_per_indiv == "yes"){
        within_silhouette_3D <- mean_sil_df_within_comb[which(mean_sil_df_within_comb$dimension == "3D"), ]
        within_silhouette_3D$method <- factor(within_silhouette_3D$method, levels = SNP_list_names)
        selected_colours <- my_colours[1:length(within_silhouette_3D$method)]
        
        within_silhouette_3D_plot <- ggplot(within_silhouette_3D, aes(x = method, y = mean_sil)) + 
          geom_bar(stat = "identity", position = position_dodge(width = 0.8), fill = selected_colours, alpha = 0.5, width = 0.8, color = "black") +
          labs(title = "Silhouette within indiv", 
               x = "method",
               y = "average silhouette score across samples") + 
          th + th_transparent + 
          theme(legend.position = "top",
                legend.box.background = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_blank(),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
          geom_hline(yintercept = 0)
        within_silhouette_3D_plot
        
        ggsave(paste0(output_folder, "Figures/Silhouette_within_3D.pdf"), within_silhouette_3D_plot, units = "cm", width = 7.5, height = 11, dpi = 400)
        
        
        
        within_silhouette_2D <- mean_sil_df_within_comb[which(mean_sil_df_within_comb$dimension == "2D"), ]
        within_silhouette_2D$method <- factor(within_silhouette_2D$method, levels = SNP_list_names)
        selected_colours <- my_colours[1:length(within_silhouette_2D$method)]
        
        within_silhouette_2D_plot <- ggplot(within_silhouette_2D, aes(x = method, y = mean_sil)) + 
          geom_bar(stat = "identity", position = position_dodge(width = 0.8), fill = selected_colours, alpha = 0.5, width = 0.8, color = "black") +
          labs(title = "Silhouette within indiv", 
               x = "method",
               y = "average silhouette score across samples") + 
          th + th_transparent + 
          theme(legend.position = "top",
                legend.box.background = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_blank(),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
          geom_hline(yintercept = 0)
        within_silhouette_2D_plot
        
        ggsave(paste0(output_folder, "Figures/Silhouette_within_2D.pdf"), within_silhouette_2D_plot, units = "cm", width = 7.5, height = 11, dpi = 400)
      }
      
      
      # COMBINED INDIVIDUAL
      
      # combine the two within individual plots - 3D
      if (multiple_samples_per_indiv == "yes"){
        design <- "AAABBBBB"
        within_silhouette_3D_plot <- within_silhouette_3D_plot + ggtitle("Silhouette Scores")
        within_plot_3D <- within_plot_3D + ggtitle("Distance to Centroid") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        layout <- (within_silhouette_3D_plot + within_plot_3D) + plot_layout(design = design) +
          plot_annotation(title = paste0("Individual Level"),
                          theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
        ggsave(paste0(output_folder, "Figures/Within_individual_comb.pdf"), layout, units = "cm", width = 15, height = 12, dpi = 400)
      }
      
      # combine the two within individual plots - 3D
      empty <- ggplot() + theme_void() #+ th + th_transparent
      empty2 <- ggplot() + theme_void() #+ th + th_transparent
      
      design <- 
        "AAA
         BBB
         CCC
         DDD
         EEE"
      silhouette_3D <- silhouette_3D + ggtitle("Silhouette Scores")
      between_cluster_3D <- between_cluster_3D + ggtitle("Distance Between Clusters")
      within_cluster_3D <- within_cluster_3D + ggtitle("Distance Within Clusters")
      layout <- (silhouette_3D + empty + between_cluster_3D + empty2 + within_cluster_3D) + 
        plot_layout(design = design, heights = c(1, 0.1, 1, 0.1, 1)) +
        plot_annotation(title = paste0("Ancestry Category Level"),
                        theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
      ggsave(paste0(output_folder, "Figures/Cluster_comb.pdf"), layout, units = "cm", width = 15, height = 30, dpi = 40)
      
      
}
 
  
  

