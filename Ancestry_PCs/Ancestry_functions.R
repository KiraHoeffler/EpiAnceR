
################################################################################
# OVERVIEW
################################################################################

# Author: Kira Hoeffler
# Part 1) Function to extract ancestry information (required resources are made available in Resources folder)
# Part 2) Function to run PCA from ancestry information





################################################################################
################################################################################
# PART I) EXTRACT ANCESTRY INFORMATION
################################################################################
################################################################################



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
#' @param dp if you have a df/matrix with detection p values already calculated, it can speed up the process, if not it takes longer to run
#'
#' @return df with ancestry information that can be used to calculate PCs on
#' 
ancestry_info <- function(background_corr_RGset, samplesheet, array_type, path_SNP_folder, cell_cols, dp = NA){
  
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
  
  suppressMessages({
    library(minfi)
    library(wateRmelon)
    library(ChAMP)
  })
  
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
    
    if (length(unique(Sex)) > 1){
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
    } else {
      if (nr_cell_type_prop_PCs == 2){
        model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + 
                      Age + 
                      Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
      }
      
      if (nr_cell_type_prop_PCs == 3){
        model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 +
                      Age + 
                      Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
      }
      
      if (nr_cell_type_prop_PCs == 4){
        model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 + cellcounts_PC4 + 
                      Age + 
                      Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
      }
      
      if (nr_cell_type_prop_PCs >= 5){
        model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 + cellcounts_PC4 + cellcounts_PC5 + 
                      Age + 
                      Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
      }
    }
    
    
    
    
    
    
    residuals <- resid(model)
    resid_beta[i, ] <- residuals
  }
  rownames(resid_beta) <- rownames(beta_SNP0bp_imp)
  colnames(resid_beta) <- colnames(beta_SNP0bp_imp)
  
  rm(model)
  
  
  ### COMBINE RESIDUALISED SNP0bp with rs beta ###
  
  rs_genotypes <- rs_genotypes[, colnames(rs_genotypes) %in% colnames(resid_beta)]
  comb_SNPs <- rbind(resid_beta, rs_genotypes)
  
  return(comb_SNPs)
  
}







################################################################################
################################################################################
# PART II) PCA
################################################################################
################################################################################

#' RUN ANCESTRY PCA
#'
#' @param SNP_info output of ancestry_info function
#' df; cpgs/snps as rows, sample names as columns:
#' structure example:
#'                123456789_R01C01    123456789_R02C01    123456789_R03C01
#'  cg03105556          0.02389115          0.01232904          0.03520888
#'  cg06248701         -0.03864776         -0.03530882         -0.09446918
#'  cg09122588          0.05765300          0.13806254          0.10284979
#' 
#' @param samplesheet samplesheet containing the samples for which ancestry PCs will be calculated;
#' df; rows as samples, variables as columns; needs to contain the column Basename that must overlap with
#' the column names of the SNP_info df
#'
#' @return PCA_result; df; samples as rows, PCs as columns
#' example structure:
#'                          PC1       PC2          PC3
#'  123456789_R01C01 -0.1519075 0.7865262 -0.686906584
#'  123456789_R02C01 -0.1886758 0.7915884 -0.736685232
#'  123456789_R03C01 -0.6935017 0.8249382  0.003974756
#' 
#' 
ancestry_PCA <- function(SNP_info, samplesheet){
  
  if (class(SNP_info) != "data.frame"){
    stop("Check that the SNP_info object is a data.frame.")
  }
  
  if (class(samplesheet) != "data.frame"){
    stop("Check that the samplesheet object is a data.frame.")
  }
  
  
  if (!"Basename" %in% colnames(samplesheet)){
    stop("The <Basename> column is missing in the samplesheet. It refers to the sample name (and the colnames of the RGset)")
  }
  
  overlapping_samples <- intersect(colnames(SNP_info), samplesheet$Basename)
  
  if (length(overlapping_samples) < 5){
    stop("The samplenames do not overlap: colnames(SNP_info) and samplesheet$Basename.")
  }
  
  if (length(overlapping_samples) < length(SNP_info)){
    print("No exact overlap between the samples of the two input objects. Only overlapping samples are used")
  }
  
  if (length(overlapping_samples) < length(samplesheet$Basename)){
    print("No exact overlap between the samples of the two input objects. Only overlapping samples are used")
  }
  
  SNP_info <- SNP_info[, overlapping_samples]
  SNP_info <- SNP_info[, order(colnames(SNP_info))]
  samplesheet <- samplesheet[which(samplesheet$Basename %in% overlapping_samples), ]
  samplesheet <- samplesheet[order(samplesheet$Basename),]
  
  pc <- prcomp(t(SNP_info))
  PCA_result <- as.data.frame(pc$x)
  PCA_result$Basename <- samplesheet$Basename
  
  return(PCA_result)
}



