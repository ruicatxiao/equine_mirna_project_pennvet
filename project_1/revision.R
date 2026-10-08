counts_path   <- "mirtop_mature_4da_counts.tsv"
samples_path  <- "smallrnaseq_samplesheet_4da.csv"


# Keep a miRNA if it has >=10 counts in >=min_samples of the 20 samples.
#   min_samples = 3  -> 117 genes   (primary analysis)
min_counts   <- 10
min_samples  <- 3


table_dir <- "table"
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
gene_name <- setNames(counts$gene_name, rownames(counts))  
counts <- as.matrix(counts[, -1])                   
storage.mode(counts) <- "integer"

# Sample sheet (condition + horse)
samplesheet <- read.csv(samples_path, stringsAsFactors = FALSE)
coldata <- data.frame(
  sample    = samplesheet$sample,
  condition = factor(samplesheet$condition),
  horse     = factor(samplesheet$horse),
  row.names = samplesheet$sample,
  stringsAsFactors = FALSE
)
coldata <- coldata[colnames(counts), ]              # same order as counts
stopifnot(all(rownames(coldata) == colnames(counts)))


keep <- rowSums(counts >= min_counts) >= min_samples
counts_f <- counts[keep, ]
cat("Filtered genes (>=10 counts in >=", min_samples, " samples):", sum(keep),
    "of", nrow(counts), "\n")


dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata,
                              design = ~ horse + condition)
dds <- DESeq(dds)

# Shrinkage of log2FC.
shrink_type <- "normal"


get_res <- function(treat, ref) {
  res <- results(dds, contrast = c("condition", treat, ref))
  shr <- lfcShrink(dds, contrast = c("condition", treat, ref), type = shrink_type)
  df <- data.frame(
    gene_id     = rownames(res),
    gene_name   = gene_name[rownames(res)],
    baseMean    = res$baseMean,
    log2FC      = shr$log2FoldChange,   # shrunken
    lfcSE       = shr$lfcSE,            # shrunken standard error
    pvalue      = res$pvalue,
    padj        = res$padj,
    stringsAsFactors = FALSE
  )
  df[order(df$padj), ]
}

res_ctrl_vs_IFN    <- get_res("IFN",     "control")
res_ctrl_vs_TNF    <- get_res("TNF",     "control")
res_ctrl_vs_IFNTNF <- get_res("IFN_TNF", "control")
res_IFN_vs_TNF     <- get_res("TNF",     "IFN")
res_IFNTNF_vs_IFN  <- get_res("IFN_TNF", "IFN")
res_IFNTNF_vs_TNF  <- get_res("IFN_TNF", "TNF")


# ---- volcano control vs IFN ------------------------------------------------
df <- res_ctrl_vs_IFN
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log2 FoldChange (shrunken)",
       y = "-log10(adjusted p-value)",
       title = "control vs IFN") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- volcano control vs TNF ------------------------------------------------
df <- res_ctrl_vs_TNF
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log2 FoldChange (shrunken)",
       y = "-log10(adjusted p-value)",
       title = "control vs TNF") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- volcano control vs IFN+TNF -------------------------------------------
df <- res_ctrl_vs_IFNTNF
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log2 FoldChange (shrunken)",
       y = "-log10(adjusted p-value)",
       title = "control vs IFN+TNF") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)



# ---- MA control vs IFN ------------------------------------------------
df <- res_ctrl_vs_IFN
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log10 baseMean", y = "log2 FoldChange (shrunken)",
       title = "MA: control vs IFN") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- MA control vs TNF ------------------------------------------------
df <- res_ctrl_vs_TNF
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log10 baseMean", y = "log2 FoldChange (shrunken)",
       title = "MA: control vs TNF") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- MA control vs IFN+TNF -------------------------------------------
df <- res_ctrl_vs_IFNTNF
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log10 baseMean", y = "log2 FoldChange (shrunken)",
       title = "MA: control vs IFN+TNF") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- MA IFN vs TNF ---------------------------------------------------
df <- res_IFN_vs_TNF
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log10 baseMean", y = "log2 FoldChange (shrunken)",
       title = "MA: IFN vs TNF") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- MA IFN+TNF vs IFN ----------------------------------------------
df <- res_IFNTNF_vs_IFN
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log10 baseMean", y = "log2 FoldChange (shrunken)",
       title = "MA: IFN+TNF vs IFN") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- MA IFN+TNF vs TNF ----------------------------------------------
df <- res_IFNTNF_vs_TNF
df$sig <- !is.na(df$padj) & df$padj <= 0.05
p <- ggplot(df, aes(x = log10(baseMean), y = log2FC)) +
  geom_point(aes(color = sig), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "padj <= 0.05"),
                     name = "") +
  geom_hline(yintercept = 0, color = "black") +
  geom_text_repel(data = subset(df, sig), aes(label = gene_name),
                  size = 3.5, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey50") +
  labs(x = "log10 baseMean", y = "log2 FoldChange (shrunken)",
       title = "MA: IFN+TNF vs TNF") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- PCA ----------------------------------------------
dds_pca <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata,
                                  design = ~ horse + condition)
vsd <- varianceStabilizingTransformation(dds_pca, blind = TRUE)
pca <- prcomp(t(assay(vsd)))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(pca$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df$condition <- coldata[pca_df$sample, "condition"]
pca_df$horse <- coldata[pca_df$sample, "horse"]

cond_lab <- c(control = "Control", IFN = "IFN-g", TNF = "TNF-a",
              IFN_TNF = "IFN-g+TNF-a")
pca_df$cond_label <- factor(cond_lab[as.character(pca_df$condition)],
                            levels = c("Control", "IFN-g", "TNF-a", "IFN-g+TNF-a"))

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cond_label)) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_manual(values = c("Control" = "#0072B2", "IFN-g" = "#D55E00",
                                "TNF-a" = "#009E73", "IFN-g+TNF-a" = "#CC79A7"),
                     name = "") +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf,
                  box.padding = 0.4, segment.color = "grey60") +
  labs(x = paste0("PC1: ", pct[1], "%"), y = paste0("PC2: ", pct[2], "%"),
       title = "PCA (VST, filtered)") +
  theme_bw() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)


# ---- Venn ----------------------------------------------
sig_IFN    <- res_ctrl_vs_IFN$gene_id[res_ctrl_vs_IFN$padj <= 0.05 & !is.na(res_ctrl_vs_IFN$padj)]
sig_TNF    <- res_ctrl_vs_TNF$gene_id[res_ctrl_vs_TNF$padj <= 0.05 & !is.na(res_ctrl_vs_TNF$padj)]
sig_IFNTNF <- res_ctrl_vs_IFNTNF$gene_id[res_ctrl_vs_IFNTNF$padj <= 0.05 & !is.na(res_ctrl_vs_IFNTNF$padj)]


if (requireNamespace("VennDiagram", quietly = TRUE)) {
  library(VennDiagram)
  venn_list <- list("IFN-g" = sig_IFN, "TNF-a" = sig_TNF, "IFN-g+TNF-a" = sig_IFNTNF)
  venn.plot <- venn.diagram(
    x = venn_list,
    fill = c("#D55E00", "#009E73", "#CC79A7"),
    alpha = 0.35,
    col = c("#D55E00", "#009E73", "#CC79A7"),
    cat.col = c("#D55E00", "#009E73", "#CC79A7"),
    cex = 1.5, cat.cex = 1.3,
    main = "Overlap of DE miRNAs (padj <= 0.05)",
    filename = NULL                       # prints to screen, no file
  )
  grid.newpage()
  grid.draw(venn.plot)

} else {
  cat("VennDiagram package not installed.\n")
  cat("  install.packages('VennDiagram')  then re-run this step.\n")
  cat("  DE gene counts:\n")
  cat("    IFN-g    :", length(sig_IFN), "\n")
  cat("    TNF-a    :", length(sig_TNF), "\n")
  cat("    IFN-g+TNF-a:", length(sig_IFNTNF), "\n")
}

# ---- Table ----------------------------------------------

# Columns for every table: gene_id, gene_name, baseMean, log2FC (shrunken),
# lfcSE, pvalue, padj
prepare_table <- function(res_df) {
  res_df[order(res_df$padj, res_df$pvalue), c("gene_id", "gene_name",
                                             "baseMean", "log2FC", "lfcSE",
                                             "pvalue", "padj")]
}


write.csv(prepare_table(res_ctrl_vs_IFN),    file.path(table_dir, "DE_control_vs_IFN_full.csv"),    row.names = FALSE)
write.csv(prepare_table(res_ctrl_vs_TNF),    file.path(table_dir, "DE_control_vs_TNF_full.csv"),    row.names = FALSE)
write.csv(prepare_table(res_ctrl_vs_IFNTNF), file.path(table_dir, "DE_control_vs_IFNTNF_full.csv"), row.names = FALSE)

# ---- Table 1 replacement ----------------------------------------------
sig_rows <- function(res_df, contrast_label) {
  s <- res_df[!is.na(res_df$padj) & res_df$padj <= 0.05, ]
  if (nrow(s) == 0) return(NULL)
  data.frame(
    Contrast  = contrast_label,
    gene_name = s$gene_name,
    baseMean  = round(s$baseMean, 1),
    log2FC    = round(s$log2FC, 2),
    padj      = signif(s$padj, 3),
    Direction = ifelse(s$log2FC > 0, "Up in treatment", "Down in treatment"),
    stringsAsFactors = FALSE
  )
}

table1 <- rbind(
  sig_rows(res_ctrl_vs_IFN,    "control vs IFN-g"),
  sig_rows(res_ctrl_vs_TNF,    "control vs TNF-a"),
  sig_rows(res_ctrl_vs_IFNTNF, "control vs IFN-g+TNF-a"),
  sig_rows(res_IFN_vs_TNF,     "IFN-g vs TNF-a"),
  sig_rows(res_IFNTNF_vs_IFN,  "IFN-g+TNF-a vs IFN-g"),
  sig_rows(res_IFNTNF_vs_TNF,  "IFN-g+TNF-a vs TNF-a")
)
write.csv(table1, file.path(table_dir, "Table1_significant_DE.csv"), row.names = FALSE)

# --- Table Summery Genes -----------------
counts_summary <- data.frame(
  Contrast = c("control vs IFN-g", "control vs TNF-a", "control vs IFN-g+TNF-a",
               "IFN-g vs TNF-a", "IFN-g+TNF-a vs IFN-g", "IFN-g+TNF-a vs TNF-a"),
  n_sig    = c(
    sum(!is.na(res_ctrl_vs_IFN$padj)    & res_ctrl_vs_IFN$padj    <= 0.05),
    sum(!is.na(res_ctrl_vs_TNF$padj)    & res_ctrl_vs_TNF$padj    <= 0.05),
    sum(!is.na(res_ctrl_vs_IFNTNF$padj) & res_ctrl_vs_IFNTNF$padj <= 0.05),
    sum(!is.na(res_IFN_vs_TNF$padj)     & res_IFN_vs_TNF$padj     <= 0.05),
    sum(!is.na(res_IFNTNF_vs_IFN$padj)  & res_IFNTNF_vs_IFN$padj  <= 0.05),
    sum(!is.na(res_IFNTNF_vs_TNF$padj)  & res_IFNTNF_vs_TNF$padj  <= 0.05)
  ),
  stringsAsFactors = FALSE
)
write.csv(counts_summary, file.path(table_dir, "Table_significant_counts.csv"), row.names = FALSE)

# done
