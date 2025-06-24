# Bash and R script used for additional analysis
# This is not an automated script, R and Bash portion needs to do executed separately
# Scripts below consists post analysis for both project_1 and project_2

library(tidyverse)
library(isomiRs)
library(pheatmap)
library(viridis) 




# isoMIR R scripts
# isoMIR project_1
mirtop <- read_tsv("joined_samples_mirtop.tsv") 
sample_metadata <- tibble(
  samplename = c("H1IIFNTNF", "H1TNF", "H3IFN",   "H3IIFNTNF", "H3TNF",
                 "H4Control", "H4IIFNTNF", "H4IFN", "H3Control", "H5IIFNTNF",
                 "H4TNF",      "H5Control", "H5TNF", "H2IIFNTNF", "H5IFN",
                 "H2IFN",      "H1IFN",      "H1Control", "H2Control", "H2TNF"),
  condition   = c("IIFNTNF", "TNF", "IFN", "IIFNTNF", "TNF",
                  "Control", "IIFNTNF", "IFN", "Control", "IIFNTNF",
                  "TNF",     "Control",  "TNF", "IIFNTNF",   "IFN",
                  "IFN",     "IFN",      "Control", "Control", "TNF")
) 

sample_metadata <- tibble(
  samplename = c("H1IIFNTNF", "H1TNF", "H3IFN",   "H3IIFNTNF", "H3TNF",
                 "H4Control", "H4IIFNTNF", "H4IFN", "H3Control", "H5IIFNTNF",
                 "H4TNF",      "H5Control", "H5TNF", "H2IIFNTNF", "H5IFN",
                 "H2IFN",      "H1IFN",      "H1Control", "H2Control", "H2TNF"),
  condition   = c("Horse_1", "Horse_1", "Horse_3", "Horse_3", "Horse_3",
                  "Horse_4", "Horse_4", "Horse_4", "Horse_3", "Horse_5",
                  "Horse_4",     "Horse_5",  "Horse_5", "Horse_2",   "Horse_5",
                  "Horse_2",     "Horse_1",      "Horse_1", "Horse_2", "Horse_2")
) 

features <- mirtop %>%
  select(
    seq  = Read,
    mir  = miRNA,
    mism = iso_snp_nt,
    add  = iso_add3p_nt,
    t5   = iso_5p_nt,
    t3   = iso_3p_nt,
    uid  = UID
  )  

samples <- names(mirtop)[13:ncol(mirtop)]
input <- bind_cols(features, mirtop[, samples])  
input_raw <- input %>%
  select(-uid) %>%
  mutate(change = as.numeric(mism != 0))

coldata <- tibble(sample = samples) %>%
  left_join(sample_metadata, by = c("sample" = "samplename")) %>%  
  column_to_rownames(var = "sample")

iso_data <- IsomirDataSeqFromRawData(
  input_raw,
  coldata = coldata
)

ids <- isoCounts(iso_data)

isoPlot(ids, type="all")


ids = isoNorm(ids, formula = ~ condition)

mat <- counts(ids, norm = TRUE)

coldata <- as.data.frame(colData(ids))
coldata$sample    <- rownames(coldata)
coldata$condition <- factor(
  coldata$condition,
  levels = c("Control", "IFN", "TNF", "IIFNTNF")
)

ord <- with(coldata, order(condition, sample))
mat_ord        <- mat[, ord]
annotation_ord <- coldata[ord, , drop = FALSE]
ann_col <- annotation_ord["condition", drop = FALSE]


pheatmap(
  mat_ord,
  annotation_col = ann_col,
  cluster_cols   = FALSE,   # disable column clustering :contentReference[oaicite:0]{index=0}
  cluster_rows   = TRUE,    # keep row clustering :contentReference[oaicite:1]{index=1}
  show_rownames  = TRUE,    # display row names :contentReference[oaicite:2]{index=2}
  scale          = "row",
  border_color   = NA
)
head(isoAnnotate(ids))
View(isoAnnotate(ids))
write.csv(isoAnnotate(ids),"ec_project1_isomir_importance_annotation.csv", row.names = TRUE)


# isoMIR project_2
mirtop <- read_tsv("joined_samples_mirtop_2.tsv") 

sample_metadata <- tibble(
  samplename = c("H3Restigen","H2BMA","H2PPP","H4Prostride","H2Prostride","H2Plasma","H5Prostride","H3Prostride","H6Plasma","H1Restigen","H6Restigen","H3CellPP","H1Prostride","H5PPP","H3Plasma","H1PPP","H4PPP","H1Plasma","H6Prostride","H1Centrate","H4Restigen","H5Centrate","H5Restigen","H4Plasma","H4BMA","H3PPP","H6BMA","H2Restigen","H1BMA","H6PPP","H1CellPP","H5Plasma","H5BMA","H5CellPP","H6CellPP","H2Centrate","H4Centrate","H2CellPP","H3BMA","H3Centrate","H6Centrate","H4CellPP"),
  condition   = c("Restigen","BMA","PPP","Prostride","Prostride","Plasma","Prostride","Prostride","Plasma","Restigen","Restigen","CellPP","Prostride","PPP","Plasma","PPP","PPP","Plasma","Prostride","Centrate","Restigen","Centrate","Restigen","Plasma","BMA","PPP","BMA","Restigen","BMA","PPP","CellPP","Plasma","BMA","CellPP","CellPP","Centrate","Centrate","CellPP","BMA","Centrate","Centrate","CellPP")
) 

sample_metadata <- tibble(
  samplename = c("H3Restigen","H2BMA","H2PPP","H4Prostride","H2Prostride","H2Plasma","H5Prostride","H3Prostride","H6Plasma","H1Restigen","H6Restigen","H3CellPP","H1Prostride","H5PPP","H3Plasma","H1PPP","H4PPP","H1Plasma","H6Prostride","H1Centrate","H4Restigen","H5Centrate","H5Restigen","H4Plasma","H4BMA","H3PPP","H6BMA","H2Restigen","H1BMA","H6PPP","H1CellPP","H5Plasma","H5BMA","H5CellPP","H6CellPP","H2Centrate","H4Centrate","H2CellPP","H3BMA","H3Centrate","H6Centrate","H4CellPP"),
  condition   = c("Horse_3","Horse_2","Horse_2","Horse_4","Horse_2","Horse_2","Horse_5","Horse_3","Horse_6","Horse_1","Horse_6","Horse_3","Horse_1","Horse_5","Horse_3","Horse_1","Horse_4","Horse_1","Horse_6","Horse_1","Horse_4","Horse_5","Horse_5","Horse_4","Horse_4","Horse_3","Horse_6","Horse_2","Horse_1","Horse_6","Horse_1","Horse_5","Horse_5","Horse_5","Horse_6","Horse_2","Horse_4","Horse_2","Horse_3","Horse_3","Horse_6","Horse_4")
) 

features <- mirtop %>%
  select(
    seq  = Read,
    mir  = miRNA,
    mism = iso_snp_nt,
    add  = iso_add3p_nt,
    t5   = iso_5p_nt,
    t3   = iso_3p_nt,
    uid  = UID
  )  

samples <- names(mirtop)[13:ncol(mirtop)]
input <- bind_cols(features, mirtop[, samples])  
input_raw <- input %>%
  select(-uid) %>%
  mutate(change = as.numeric(mism != 0))

coldata <- tibble(sample = samples) %>%
  left_join(sample_metadata, by = c("sample" = "samplename")) %>%  
  column_to_rownames(var = "sample")

iso_data <- IsomirDataSeqFromRawData(
  input_raw,
  coldata = coldata
)


ids <- isoCounts(iso_data)
isoPlot(ids, type="all")
isoPlot(ids, type="snv")
ids = isoNorm(ids, formula = ~ condition)
mat <- counts(ids, norm = TRUE)
coldata <- as.data.frame(colData(ids))
coldata$sample    <- rownames(coldata)
coldata$condition <- factor(
  coldata$condition,
  levels = c("BMA","CellPP","Centrate","PPP","Plasma","Prostride","Restigen")
)

ord <- with(coldata, order(condition, sample))
mat_ord        <- mat[, ord]
annotation_ord <- coldata[ord, , drop = FALSE]
ann_col <- annotation_ord["condition", drop = FALSE]

pheatmap(
  mat_ord,
  annotation_col = ann_col,
  cluster_cols   = FALSE,   # disable column clustering :contentReference[oaicite:0]{index=0}
  cluster_rows   = TRUE,    # keep row clustering :contentReference[oaicite:1]{index=1}
  show_rownames  = TRUE,    # display row names :contentReference[oaicite:2]{index=2}
  scale          = "row",
  border_color   = NA
)


head(isoAnnotate(ids))
View(isoAnnotate(ids))
write.csv(isoAnnotate(ids),"ec_project2_isomir_importance_annotation.csv", row.names = TRUE)



# Need to run custom Python script to transform mirdeep2 report before it can be used

python process_mirdeep2.py



# mirDeep2 analysis R scrtipts
# Annotated miRNA expressed vs non expressed
# Project_1
known_ids <- readLines("known_id.txt") %>% trimws()
report_files <- list.files(pattern = "^report_.*\\.txt$")
process_report <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_known <- df %>% filter(type == "known") %>% select(sample, ID)
  return(df_known)
}
df_list <- lapply(report_files, process_report)
df_combined <- bind_rows(df_list)

samples <- unique(df_combined$sample)
presence_matrix <- matrix(0, nrow = length(known_ids), ncol = length(samples),
                          dimnames = list(known_ids, samples))

for(i in 1:nrow(df_combined)) {
  id <- df_combined$ID[i]
  sample <- df_combined$sample[i]
  if(id %in% known_ids) {
    presence_matrix[id, sample] <- 1
  }
}

presence_df <- as.data.frame(presence_matrix) %>% 
  rownames_to_column(var = "miRNA") %>% 
  gather(key = "sample", value = "presence", -miRNA)

presence_df$presence <- factor(presence_df$presence, levels = c(0, 1))
presence_df$miRNA <- factor(presence_df$miRNA, levels = sort(unique(presence_df$miRNA)))
presence_df$sample <- factor(presence_df$sample, levels = sort(unique(presence_df$sample)))

ggplot(presence_df, aes(x = sample, y = miRNA, fill = presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  labs(title = "Known miRNA Presence Heatmap") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
hc_rows <- hclust(dist(presence_matrix))
row_order <- rownames(presence_matrix)[hc_rows$order]
hc_cols <- hclust(dist(t(presence_matrix)))
col_order <- colnames(presence_matrix)[hc_cols$order]
presence_df$miRNA <- factor(presence_df$miRNA, levels = row_order)
presence_df$sample <- factor(presence_df$sample, levels = col_order)
ggplot(presence_df, aes(x = sample, y = miRNA, fill = presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  labs(title = "Known miRNA Presence Heatmap (Clustered)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# Project_2
process_report <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_known <- df %>% filter(type == "known") %>% select(sample, ID)
  return(df_known)
}

known_ids <- readLines("known_id.txt") %>% trimws()
report_files <- list.files(pattern = "^report_.*\\.txt$")

df_list <- lapply(report_files, process_report)
df_combined <- bind_rows(df_list)

samples <- unique(df_combined$sample)
presence_matrix <- matrix(0, nrow = length(known_ids), ncol = length(samples),
                          dimnames = list(known_ids, samples))

for(i in 1:nrow(df_combined)) {
  id <- df_combined$ID[i]
  sample <- df_combined$sample[i]
  if(id %in% known_ids) {
    presence_matrix[id, sample] <- 1
  }
}

presence_df <- as.data.frame(presence_matrix) %>% 
  rownames_to_column(var = "miRNA") %>% 
  gather(key = "sample", value = "presence", -miRNA)


presence_df$presence <- factor(presence_df$presence, levels = c(0, 1))


presence_df$miRNA <- factor(presence_df$miRNA, levels = sort(unique(presence_df$miRNA)))
presence_df$sample <- factor(presence_df$sample, levels = sort(unique(presence_df$sample)))

ggplot(presence_df, aes(x = sample, y = miRNA, fill = presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  labs(title = "Known miRNA Presence Heatmap") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


hc_rows <- hclust(dist(presence_matrix))
row_order <- rownames(presence_matrix)[hc_rows$order]


hc_cols <- hclust(dist(t(presence_matrix)))
col_order <- colnames(presence_matrix)[hc_cols$order]

presence_df$miRNA <- factor(presence_df$miRNA, levels = row_order)
presence_df$sample <- factor(presence_df$sample, levels = col_order)

ggplot(presence_df, aes(x = sample, y = miRNA, fill = presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  labs(title = "Known miRNA Presence Heatmap (Clustered)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# mirDeep2 analysis R scrtipts
# Novel miRNA expressed vs non expressed in bar chart
# Project_1
report_files <- list.files(pattern = "^report_.*\\.txt$")
df_all <- map_df(report_files, ~ read.delim(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
df_counts <- df_all %>%
  group_by(sample, type) %>%
  summarise(count = n(), .groups = "drop")
desired_order <- c("H1Control", "H2Control", "H3Control", "H4Control", "H5Control", "H1IFN", "H2IFN", "H3IFN", "H4IFN", "H5IFN", "H1TNF", "H2TNF", "H3TNF", "H4TNF", "H5TNF","H1IIFNTNF", "H2IIFNTNF", "H3IIFNTNF", "H4IIFNTNF","H5IIFNTNF")

df_counts$sample <- factor(df_counts$sample, levels = desired_order)
ggplot(df_counts, aes(x = sample, y = count, fill = forcats::fct_rev(type))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Project_2
report_files <- list.files(pattern = "^report_.*\\.txt$")
df_all <- map_df(report_files, ~ read.delim(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE))

df_counts <- df_all %>%
  group_by(sample, type) %>%
  summarise(count = n(), .groups = "drop")

desired_order <- c("H1BMA","H2BMA","H3BMA","H4BMA","H5BMA","H6BMA","H1CellPP","H2CellPP","H3CellPP","H4CellPP","H5CellPP","H6CellPP","H1Centrate","H2Centrate","H3Centrate","H4Centrate","H5Centrate","H6Centrate","H1PPP","H2PPP","H3PPP","H4PPP","H5PPP","H6PPP","H1Plasma","H2Plasma","H3Plasma","H4Plasma","H5Plasma","H6Plasma","H1Prostride","H2Prostride","H3Prostride","H4Prostride","H5Prostride","H6Prostride","H1Restigen","H2Restigen","H3Restigen","H4Restigen","H5Restigen","H6Restigen")
df_counts$sample <- factor(df_counts$sample, levels = desired_order)

#    The bars will show counts of known and novel miRNAs for each sample
ggplot(df_counts, aes(x = sample, y = count, fill = forcats::fct_rev(type))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# mirDeep2 analysis Bash scrtipts
# Novel miRNA similarity among different samples

# Project_1
cat report*.txt \
    | awk '$2=="novel"' | cut -f 4 \
    | awk -F ':' '{print $1"_"$2}' \
    | sort -V | uniq -c | awk '{print $2,$1}' OFS='\t' \
    | sort -k2,2rn

cat report*.txt \
    | awk '$2=="novel"' | cut -f 4 \
    | awk -F ':' '{print $1"_"$2}' \
    | sed 's/\.\./\_/g' \
    | sed 's/\_/\t/g' \
    | sort -k1,1 -k2,2n > novel.bed


bedtools merge -d 10 -i novel.bed | awk '{print $1,$2,$3,"novel_"NR}' OFS='\t' > novel_merged.bed
for file in report_*.txt; do
    echo
    name1=$(basename "$file" .txt)
    name2=${name1/#report_}
    cat ${file} \
        | awk '$2=="novel"' | cut -f 4 \
        | awk -F ':' '{print $1"_"$2}' \
        | sed 's/\.\./\_/g' \
        | sed 's/\_/\t/g' \
        | sort -k1,1 -k2,2n \
        | awk -v var=${name2} '{print $1,$2,$3,var"_"NR}' OFS='\t' > "samples_${name2}_novel.bed"
    echo
done

bedtools intersect -a novel_merged.bed -b samples_H1Control_novel.bed \
-c  -F 1 | cut -f 4,5


for region in samples_*_novel.bed; do
    echo
    name1=$(basename "$region" _novel.bed)
    name2=${name1/#samples_} 
    bedtools intersect \
        -a novel_merged.bed \
        -b "${region}" -c -F 1 \
            | cut -f 4,5 | awk -F '/s' -v var=${name2} '{print var,$1.$2}' OFS='\t' \
            | awk '{if ($3 >=1) print $1,$2,1; else print $1,$2,0}' OFS='\t' >> novel_miRNA_all_samples_count.txt
    echo
done

# Project_2
cat report*.txt \
    | awk '$2=="novel"' | cut -f 4 \
    | awk -F ':' '{print $1"_"$2}' \
    | sort -V | uniq -c | awk '{print $2,$1}' OFS='\t' \
    | sort -k2,2rn

cat report*.txt \
    | awk '$2=="novel"' | cut -f 4 \
    | awk -F ':' '{print $1"_"$2}' \
    | sed 's/\.\./\_/g' \
    | sed 's/\_/\t/g' \
    | sort -k1,1 -k2,2n > novel.bed


bedtools merge -d 10 -i novel.bed \
    | awk '{print $1,$2,$3,"novel_"NR}' OFS='\t' > novel_merged.bed

for file in report_*.txt; do
    echo
    name1=$(basename "$file" .txt)
    name2=${name1/#report_}
    cat ${file} \
        | awk '$2=="novel"' | cut -f 4 \
        | awk -F ':' '{print $1"_"$2}' \
        | sed 's/\.\./\_/g' \
        | sed 's/\_/\t/g' \
        | sort -k1,1 -k2,2n \
        | awk -v var=${name2} '{print $1,$2,$3,var"_"NR}' OFS='\t' > "samples_${name2}_novel.bed"
    echo
done

bedtools intersect -a novel_merged.bed -b samples_H1Control_novel.bed \
-c  -F 1 | cut -f 4,5
for region in samples_*_novel.bed; do
    echo

    name1=$(basename "$region" _novel.bed)
    name2=${name1/#samples_}
    bedtools intersect \
        -a novel_merged.bed \
        -b "${region}" -c -F 1 \
            | cut -f 4,5 | awk -F '/s' -v var=${name2} '{print var,$1.$2}' OFS='\t' \
            | awk '{if ($3 >=1) print $1,$2,1; else print $1,$2,0}' OFS='\t' >> novel_miRNA_all_samples_count.txt
    echo
done



# mirDeep2 analysis R scrtipts
# Novel miRNA similarity among different samples
# Project_1
df <- read.delim("novel_miRNA_all_samples_count.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
wide_df <- df %>%
  pivot_wider(names_from = sample, values_from = count)
region_names <- wide_df$region
wide_matrix <- as.matrix(wide_df[,-1])
rownames(wide_matrix) <- region_names
hc_rows <- hclust(dist(wide_matrix))
hc_cols <- hclust(dist(t(wide_matrix)))
row_order <- rownames(wide_matrix)[hc_rows$order]
col_order <- colnames(wide_matrix)[hc_cols$order]
df_long <- as.data.frame(wide_matrix) %>%
  rownames_to_column(var = "region") %>%
  pivot_longer(-region, names_to = "sample", values_to = "count")

df_long$region <- factor(df_long$region, levels = row_order)
df_long$sample <- factor(df_long$sample, levels = col_order)
ggplot(df_long, aes(x = sample, y = region, fill = factor(count))) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue"),
                    name = "Presence",
                    labels = c("0" = "Absent", "1" = "Present")) +
  labs(title = "Heatmap of Novel miRNA Presence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# Project_2
df <- read.delim("novel_miRNA_all_samples_count.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
wide_df <- df %>%
  pivot_wider(names_from = sample, values_from = count)
region_names <- wide_df$region
wide_matrix <- as.matrix(wide_df[,-1])
rownames(wide_matrix) <- region_names

hc_rows <- hclust(dist(wide_matrix))
hc_cols <- hclust(dist(t(wide_matrix)))
row_order <- rownames(wide_matrix)[hc_rows$order]
col_order <- colnames(wide_matrix)[hc_cols$order]


df_long <- as.data.frame(wide_matrix) %>%
  rownames_to_column(var = "region") %>%
  pivot_longer(-region, names_to = "sample", values_to = "count")

df_long$region <- factor(df_long$region, levels = row_order)
df_long$sample <- factor(df_long$sample, levels = col_order)
ggplot(df_long, aes(x = sample, y = region, fill = factor(count))) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue"),
                    name = "Presence",
                    labels = c("0" = "Absent", "1" = "Present")) +
  labs(title = "Heatmap of Novel miRNA Presence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


# Intersection plot R script
# Project_1 intersect plot
df1 <- read.delim("ec_mature_4upset.txt", header = T, row.names="miRNA")
upset(df1, sets = c("IFN_vs_Control_UP","IFN_vs_Control_DOWN","TNF_vs_Control_UP","TNF_vs_Control_DOWN","INFTNF_vs_Control_UP","INFTNF_vs_Control_DOWN","TNF_vs_INF_UP","TNF_vs_INF_DOWN","IFN_vs_IFNTNF_UP" ,"IFN_vs_IFNTNF_DOWN","TNF_vs_IFNTNF_UP" ,"TNF_vs_IFNTNF_DOWN"), order.by = "degree", keep.order = FALSE, point.size = 5, line.size = 1, text.scale=2, sets.x.label="miRNA", nintersects = NA, show.numbers = TRUE)

# Project_2 intersect plot
df1 <- read.delim("ec_mature_4upset_2.txt", header = T, row.names="miRNA")
upset(df1, sets = c("CellPP_vs_BMA_DOWN","Centrate_vs_BMA_DOWN","Prostride_vs_Plasma_UP","Prostride_vs_Plasma_DOWN","Restigen_vs_Plasma_UP","Restigen_vs_Plasma_DOWN","Plasma_vs_PPP_UP","Plasma_vs_PPP_DOWN","Restigen_vs_PPP_UP","Restigen_vs_PPP_DOWN"), order.by = "degree", keep.order = FALSE, point.size = 5, line.size = 1, text.scale=3, sets.x.label="miRNA", nintersects = NA, show.numbers = TRUE)








