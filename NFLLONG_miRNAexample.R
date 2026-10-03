###############################################################################################################
#                       Example Code for original article submitted to Brain Communications                   #
#     Cumulative Neurotrauma and Neurodegeneration-Related MicroRNA in Former Professional Football Players   #
#                             Corresponding Author: Benjamin Brett (bbrett@mcw.edu)                           #
#                                           Code Author: Colin Hansen                                         #
###############################################################################################################


################################################################################
############################# DATA PREP ########################################
################################################################################

# Load Packacges
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(multiMiR)
library(miRNAtap)
library(miRNAtap.db)
library(AnnotationDbi)
library(pathview)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(patchwork)
library(ggtext)
library(rmarkdown)
library(corrplot)

head(mirna_data[1:5, 1:5])
head(metadata)
colnames(metadata)

mirna_data <- read.csv(mirna_file, row.names = 1)

metadata <- read.csv(metadata_file,stringsAsFactors = FALSE)

head(mirna_data[1:5, 1:5])
head(metadata)
colnames(metadata)

expression_samples <- colnames(mirna_data)

sample_info <- metadata[match(expression_samples, metadata$miRNAID), ]

cat("Expression samples:", length(expression_samples), "\n")
cat("Matched metadata rows:", nrow(sample_info), "\n")
cat("Missing matches:", sum(is.na(sample_info$miRNAID)), "\n")

head(data.frame(
  expression  = expression_samples,
  metadata_id = sample_info$miRNAID,
  YrsPart     = sample_info$TotalYrsPlayed,
  LowHighCx   = sample_info$LowHighCx
))
export_path <- export_path
if (!dir.exists(export_path)) dir.create(export_path, recursive = TRUE)

################################################################################
################DIFFERENTIAL EXPRESSION YEARS OF PLAY ##########################
################################################################################

YrsPart       <- sample_info$TotalYrsPlayed
design_yrs    <- model.matrix(~ YrsPart)
colnames(design_yrs) <- c("Intercept", "YrsPart")
dge_yrs  <- DGEList(counts = mirna_data)
keep_yrs <- filterByExpr(dge_yrs, design = design_yrs)
dge_yrs  <- dge_yrs[keep_yrs, , keep.lib.sizes = FALSE]
dge_yrs  <- calcNormFactors(dge_yrs)
dge_yrs  <- estimateDisp(dge_yrs, design_yrs)
fit_yrs  <- glmQLFit(dge_yrs, design_yrs)
qlf_yrs      <- glmQLFTest(fit_yrs, coef = 2)
results_yrs  <- topTags(qlf_yrs, n = Inf)

cat("FDR < 0.10 :", sum(results_yrs$table$FDR   < 0.1),   "\n")
cat("P   < 0.01 :", sum(results_yrs$table$PValue < 0.01),  "\n")

sig_yrs_fdr10  <- results_yrs$table[results_yrs$table$FDR    < 0.10, ]
sig_yrs_p01    <- results_yrs$table[results_yrs$table$PValue < 0.01, ]

cat("Significant miRNAs (p < 0.01):", nrow(sig_yrs_p01), "\n")

write.csv(sig_yrs_fdr10,  file.path(desktop_path, "yrs_significant_miRNAs_FDR_0.10.csv"),  row.names = TRUE)
write.csv(sig_yrs_p01,    file.path(desktop_path, "yrs_significant_miRNAs_P_0.01.csv"),    row.names = TRUE)

################################################################################
##############DIFFERENTIAL EXPRESSION CONCUSSION HISTORY #######################
################################################################################

condition    <- factor(sample_info$LowHighCx, levels = c(0, 1), labels = c("Low", "High"))
design_cx    <- model.matrix(~ 0 + condition)
colnames(design_cx) <- levels(condition)
cat("Group counts:\n")
print(table(condition))
dge_cx  <- DGEList(counts = mirna_data)
keep_cx <- filterByExpr(dge_cx, design = design_cx)
dge_cx  <- dge_cx[keep_cx, , keep.lib.sizes = FALSE]
dge_cx  <- calcNormFactors(dge_cx)
dge_cx  <- estimateDisp(dge_cx, design_cx)
fit_cx  <- glmQLFit(dge_cx, design_cx)
contrast_cx <- makeContrasts(High - Low, levels = design_cx)
qlf_cx      <- glmQLFTest(fit_cx, contrast = contrast_cx)
results_cx  <- topTags(qlf_cx, n = Inf)

cat("FDR < 0.10 :", sum(results_cx$table$FDR   < 0.1),   "\n")
cat("P   < 0.01 :", sum(results_cx$table$PValue < 0.01),  "\n")

sig_cx_fdr10  <- results_cx$table[results_cx$table$FDR    < 0.10, ]
sig_cx_p01    <- results_cx$table[results_cx$table$PValue < 0.01, ]
cat("Significant miRNAs (p < 0.01):", nrow(sig_cx_p01), "\n")


write.csv(sig_cx_fdr10,  file.path(desktop_path, "cx_significant_miRNAs_FDR_0.10.csv"),  row.names = TRUE)
write.csv(sig_cx_p01,    file.path(desktop_path, "cx_significant_miRNAs_P_0.01.csv"),    row.names = TRUE)

cat("\n=== miRNA Count Summary ===\n")
cat("Raw miRNAs in dataset:              ", nrow(mirna_data), "\n")
cat("After filtering (YoP analysis):     ", nrow(dge_yrs), "\n")
cat("After filtering (ConHx analysis):   ", nrow(dge_cx), "\n")
cat("Tested in YoP results:              ", nrow(results_yrs$table), "\n")
cat("Tested in ConHx results:            ", nrow(results_cx$table), "\n")

################################################################################
####################PATHWAY ANALYSIS YEARS OF PLAY #############################
################################################################################


sig_mirna_names <- rownames(sig_yrs_p01)
head(sig_mirna_names)

target_genes <- get_multimir(
  org = "hsa",
  mirna = sig_mirna_names,
  table = "validated",  
  summary = TRUE
)
target_gene_symbols <- unique(target_genes@data$target_symbol)
cat("Number of target genes:", length(target_gene_symbols), "\n")
entrez_ids <- bitr(target_gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",  
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                           organism = "hsa",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
head(go_results@result)
head(kegg_results@result)
export_path <- export_path
cat("\n=== Exporting Pathway Analysis Results ===\n")
if (!is.null(go_results) && nrow(go_results@result) > 0) {
  write.csv(go_results@result, 
            file.path(export_path, "GO_enrichment_results.csv"), 
            row.names = FALSE)
  cat("✓ Exported: GO_enrichment_results.csv\n")
  cat("  Number of enriched GO terms:", nrow(go_results@result), "\n")
} else {
  cat("⚠ No GO enrichment results to export (no significant terms found)\n")
}
if (!is.null(kegg_results) && nrow(kegg_results@result) > 0) {
  write.csv(kegg_results@result, 
            file.path(export_path, "KEGG_pathway_results.csv"), 
            row.names = FALSE)
  cat("✓ Exported: KEGG_pathway_results.csv\n")
  cat("  Number of enriched KEGG pathways:", nrow(kegg_results@result), "\n")
} else {
  cat("⚠ No KEGG pathway results to export (no significant pathways found)\n")
}
write.csv(target_genes@data, 
          file.path(export_path, "miRNA_target_genes.csv"), 
          row.names = FALSE)
cat("✓ Exported: miRNA_target_genes.csv\n")
write.csv(data.frame(Gene_Symbol = target_gene_symbols), 
          file.path(export_path, "unique_target_genes.csv"), 
          row.names = FALSE)
cat("✓ Exported: unique_target_genes.csv\n")
cat("\n=== Saving Plots ===\n")

################################################################################
##################PATHWAY ANALYSIS CONCUSSION HISTORY###########################
################################################################################

sig_mirna_names <- rownames(sig_cx_p01)
head(sig_mirna_names)
target_genes <- get_multimir(
  org = "hsa",
  mirna = sig_mirna_names,
  table = "validated",  
  summary = TRUE
)
target_gene_symbols <- unique(target_genes@data$target_symbol)
cat("Number of target genes:", length(target_gene_symbols), "\n")
entrez_ids <- bitr(target_gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                           organism = "hsa",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
head(go_results@result)
head(kegg_results@result)
export_path <- export_path
cat("\n=== Exporting Pathway Analysis Results ===\n")
if (!is.null(go_results) && nrow(go_results@result) > 0) {
  write.csv(go_results@result, 
            file.path(export_path, "GO_enrichment_results.csv"), 
            row.names = FALSE)
  cat("✓ Exported: GO_enrichment_results.csv\n")
  cat("  Number of enriched GO terms:", nrow(go_results@result), "\n")
} else {
  cat("⚠ No GO enrichment results to export (no significant terms found)\n")
}
if (!is.null(kegg_results) && nrow(kegg_results@result) > 0) {
  write.csv(kegg_results@result, 
            file.path(export_path, "KEGG_pathway_results.csv"), 
            row.names = FALSE)
  cat("✓ Exported: KEGG_pathway_results.csv\n")
  cat("  Number of enriched KEGG pathways:", nrow(kegg_results@result), "\n")
} else {
  cat("⚠ No KEGG pathway results to export (no significant pathways found)\n")
}
write.csv(target_genes@data, 
          file.path(export_path, "miRNA_target_genes.csv"), 
          row.names = FALSE)
cat("✓ Exported: miRNA_target_genes.csv\n")
write.csv(data.frame(Gene_Symbol = target_gene_symbols), 
          file.path(export_path, "unique_target_genes.csv"), 
          row.names = FALSE)
cat("✓ Exported: unique_target_genes.csv\n")
cat("\n=== Saving Plots ===\n")

################################################################################
####################PET & MIRNA SPEARMAN CORRELATION ###########################
################################################################################

cor_matrix <- cor(sample_info[, c("WholeBrainPiB2", "WholeBrainTau2", "FrontalPiB2", "FrontalTau2", "hsamiR206", "hsamiR133a3p", "hsamiR13p", 
                           "hsamiR199b3p", "hsamiR1505p", "hsamiR34a5p", "hsamiR34c5p", "hsamiR1291",
                           "hsamiR27a5p", "hsamiR1248", "hsamiR146a5p", "hsamiR1915p", "hsamiR30e3p",
                           "hsamiR1425p", "hsamiR2235p", "hsamiR1403p", "hsamiR3383p", "hsamiR253p",
                           "hsamiR30655p")], 
                  method = "spearman", 
                  use = "complete.obs")

# Print the matrix
print(round(cor_matrix, 2))

################################################################################
##########SENSITIVITY ANALYSES W/ ADDITIONAL CORRELATIONS ######################
################################################################################
# Years of Participation

export_path <- export_path
if (!dir.exists(export_path)) dir.create(export_path, recursive = TRUE)

mirna_data <- read.csv(mirna_file, row.names = 1)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

expression_samples <- colnames(mirna_data)
sample_info <- metadata[match(expression_samples, metadata$miRNAID), ]

cat("Expression samples :", length(expression_samples), "\n")
cat("Matched metadata   :", nrow(sample_info), "\n")
cat("Missing matches    :", sum(is.na(sample_info$miRNAID)), "\n")

YrsPart <- sample_info$TotalYrsPlayed
Age     <- sample_info$AgeatGHS        
Race    <- factor(sample_info$Race) 

cat("\nCovariate missingness:\n")
cat("  YrsPart missing:", sum(is.na(YrsPart)), "\n")
cat("  Age missing    :", sum(is.na(Age)),     "\n")
cat("  Race missing   :", sum(is.na(Race)),    "\n")

run_de_model <- function(counts, design, model_label) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("MODEL:", model_label, "\n")
  cat(strrep("=", 60), "\n")
  
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design = design)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)
  dge  <- estimateDisp(dge, design)
  
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)          
  
  results <- topTags(qlf, n = Inf)$table
  
  cat("\nThreshold summary:\n")
  cat("  FDR < 0.10 :", sum(results$FDR    < 0.10),  "\n")
  cat("  p   < 0.01 :", sum(results$PValue < 0.01),  "\n")
  
  sigs <- list(
    fdr10  = results[results$FDR    < 0.10,  ],
    p01    = results[results$PValue < 0.01,   ]
  )
  
  return(list(full = results, sigs = sigs))
}

save_model_results <- function(model_results, model_label, export_path) {
  
  threshold_map <- list(
    fdr10  = "FDR_0.10",
    p01    = "P_0.01"
  )
  
  safe_label <- gsub("[^A-Za-z0-9]", "_", model_label)
  
  for (thresh_key in names(threshold_map)) {
    fname <- file.path(
      export_path,
      paste0("sig_miRNAs_", safe_label, "_", threshold_map[[thresh_key]], ".csv")
    )
    write.csv(model_results$sigs[[thresh_key]], fname, row.names = TRUE)
    cat("  Saved:", basename(fname),
        "(n =", nrow(model_results$sigs[[thresh_key]]), ")\n")
  }
}

design_M1 <- model.matrix(~ YrsPart)                    
design_M2 <- model.matrix(~ YrsPart + Age)              
design_M3 <- model.matrix(~ YrsPart + Race)             
design_M4 <- model.matrix(~ YrsPart + Age + Race)      

cat("\nDesign matrix column names:\n")
cat("  M1:", colnames(design_M1), "\n")
cat("  M2:", colnames(design_M2), "\n")
cat("  M3:", colnames(design_M3), "\n")
cat("  M4:", colnames(design_M4), "\n")

res_M1 <- run_de_model(mirna_data, design_M1, "YrsPart_only")
res_M2 <- run_de_model(mirna_data, design_M2, "YrsPart_plus_Age")
res_M3 <- run_de_model(mirna_data, design_M3, "YrsPart_plus_Race")
res_M4 <- run_de_model(mirna_data, design_M4, "YrsPart_plus_Age_plus_Race")

export_path <- export_path

cat("\n", strrep("=", 60), "\n")
cat("SAVING RESULTS TO:", export_path, "\n")
cat(strrep("=", 60), "\n\n")

cat("--- Model 1: YrsPart only ---\n")
save_model_results(res_M1, "YrsPart_only",               export_path)

cat("--- Model 2: YrsPart + Age ---\n")
save_model_results(res_M2, "YrsPart_plus_Age",            export_path)

cat("--- Model 3: YrsPart + Race ---\n")
save_model_results(res_M3, "YrsPart_plus_Race",           export_path)

cat("--- Model 4: YrsPart + Age + Race ---\n")
save_model_results(res_M4, "YrsPart_plus_Age_plus_Race",  export_path)

write.csv(res_M1$full, file.path(export_path, "full_results_YrsPart_only.csv"),               row.names = TRUE)
write.csv(res_M2$full, file.path(export_path, "full_results_YrsPart_plus_Age.csv"),            row.names = TRUE)
write.csv(res_M3$full, file.path(export_path, "full_results_YrsPart_plus_Race.csv"),           row.names = TRUE)
write.csv(res_M4$full, file.path(export_path, "full_results_YrsPart_plus_Age_plus_Race.csv"),  row.names = TRUE)

cat("\nFull results tables also saved (all miRNAs, unfiltered).\n")


# ==============================================================================
#Concussion History 

head(data.frame(
  expression = expression_samples,
  metadata   = sample_info$miRNAID,
  LowHighCx  = sample_info$LowHighCx
))

condition <- factor(sample_info$LowHighCx, levels = c(0, 1), labels = c("Low", "High"))
Age       <- sample_info$AgeatGHS          

cat("\nCovariate missingness:\n")
cat("  condition missing:", sum(is.na(condition)), "\n")
cat("  Age missing:      ", sum(is.na(Age)),       "\n")
cat("  Race missing:     ", sum(is.na(Race)),       "\n")

cat("\nCondition breakdown:\n")
print(table(condition))
cat("\nRace levels (sanitised):\n")
print(levels(Race))

run_de_model <- function(counts, design, model_label) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("MODEL:", model_label, "\n")
  cat(strrep("=", 60), "\n")
  cat("Design matrix columns:", colnames(design), "\n")
  
  dge  <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design = design)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)
  dge  <- estimateDisp(dge, design)
  
  fit             <- glmQLFit(dge, design)
  contrast_matrix <- makeContrasts(High - Low, levels = design)
  qlf             <- glmQLFTest(fit, contrast = contrast_matrix)
  results         <- topTags(qlf, n = Inf)$table
  
  cat("\nThreshold summary:\n")
  cat("  FDR < 0.10 :", sum(results$FDR    < 0.10),  "\n")
  cat("  p   < 0.01 :", sum(results$PValue < 0.01),  "\n")
  
  sigs <- list(
    fdr10  = results[results$FDR    < 0.10,  ],
    p01    = results[results$PValue < 0.01,   ]
  )
  
  return(list(full = results, sigs = sigs))
}

save_model_results <- function(model_results, model_label, out_path) {
  
  threshold_map <- list(
    fdr10  = "FDR_0.10",
    p01    = "P_0.01"
  )
  
  safe_label <- gsub("[^A-Za-z0-9]", "_", model_label)
  
  for (thresh_key in names(threshold_map)) {
    fname <- file.path(
      out_path,
      paste0("sig_miRNAs_", safe_label, "_", threshold_map[[thresh_key]], ".csv")
    )
    write.csv(model_results$sigs[[thresh_key]], fname, row.names = TRUE)
    cat("  Saved:", basename(fname),
        "(n =", nrow(model_results$sigs[[thresh_key]]), ")\n")
  }
}

design_M1 <- model.matrix(~ 0 + condition)
colnames(design_M1)[1:2] <- c("Low", "High")

design_M2 <- model.matrix(~ 0 + condition + Age)
colnames(design_M2)[1:2] <- c("Low", "High")

design_M3 <- model.matrix(~ 0 + condition + Race)
colnames(design_M3)[1:2] <- c("Low", "High")

design_M4 <- model.matrix(~ 0 + condition + Age + Race)
colnames(design_M4)[1:2] <- c("Low", "High")

res_M1 <- run_de_model(mirna_data, design_M1, "ConHx_only")
res_M2 <- run_de_model(mirna_data, design_M2, "ConHx_plus_Age")
res_M3 <- run_de_model(mirna_data, design_M3, "ConHx_plus_Race")
res_M4 <- run_de_model(mirna_data, design_M4, "ConHx_plus_Age_plus_Race")

cat("\n", strrep("=", 60), "\n")
cat("SAVING RESULTS TO:", export_path, "\n")
cat(strrep("=", 60), "\n\n")

cat("--- Model 1: ConHx only ---\n")
save_model_results(res_M1, "ConHx_only",               export_path)

cat("--- Model 2: ConHx + Age ---\n")
save_model_results(res_M2, "ConHx_plus_Age",            export_path)

cat("--- Model 3: ConHx + Race ---\n")
save_model_results(res_M3, "ConHx_plus_Race",           export_path)

cat("--- Model 4: ConHx + Age + Race ---\n")
save_model_results(res_M4, "ConHx_plus_Age_plus_Race",  export_path)

write.csv(res_M1$full, file.path(export_path, "full_results_ConHx_only.csv"),               row.names = TRUE)
write.csv(res_M2$full, file.path(export_path, "full_results_ConHx_plus_Age.csv"),            row.names = TRUE)
write.csv(res_M3$full, file.path(export_path, "full_results_ConHx_plus_Race.csv"),           row.names = TRUE)
write.csv(res_M4$full, file.path(export_path, "full_results_ConHx_plus_Age_plus_Race.csv"),  row.names = TRUE)

cat("\nFull results tables saved (all miRNAs, unfiltered).\n")
