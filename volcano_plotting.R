# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)


res_table <- read.table("INPUT.TSV", header = TRUE, sep = "\t")


res_table <- res_table %>%
  mutate(
    diffexpressed = case_when(
      padj <= 0.05 & abs(log2FoldChange) > 1 ~ "padj<=0.05 & |LFC|>1",
      padj <= 0.05 ~ "padj<=0.05",
      TRUE ~ "Not significant"
    )
  )

res_table <- res_table %>%
  mutate(
    label = if_else(padj <= 0.05, as.character(gene_id), "")
  )

res_table$diffexpressed <- factor(res_table$diffexpressed, 
                                  levels = c("padj<=0.05 & |LFC|>1", "padj<=0.05", "Not significant"))


my_colors <- c("padj<=0.05 & |LFC|>1" = "red", 
               "padj<=0.05" = "blue", 
               "Not significant" = "grey")


ggplot(res_table, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = diffexpressed), size = 2, alpha = 0.9) + 
  scale_color_manual(values = my_colors, name = "Expression Status") +
  coord_cartesian(ylim = c(0, 10), xlim = c(-30, 30)) +
  scale_x_continuous(breaks = seq(-30, 30, by = 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  labs(x = expression("log"[2]*"FoldChange"), 
       y = expression("-log"[10]*"(adjusted p-value)")) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 18), axis.text= element_text(size = 12)) +
  geom_text_repel(
    data = subset(res_table, padj <= 0.05), 
    aes(label = label),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = Inf, 
    segment.color = 'grey50'
  )


# Set to 900 by 900 in size
# Saved to svg, png and pdf
# Set to control_vs_ifn, control_vs_tnf, control_vs_ifntnf
