# Load TCGAbiolinks
library(TCGAbiolinks)

# Create a query for prostate cancer
query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
library(TCGAbiolinks)
GDCdownload(query)
data <- GDCprepare(query)
dim(data)
data
library(SummarizedExperiment)
assay(data)[1:5, 1:5]
colnames(data)
table(data$sample_type)
counts <- assay(data)
dim(counts)
counts[1:5, 1:5]
colData(data)
group <- data$sample_type
table(group)
length(group)
ncol(counts)
BiocManager::install("DESeq2")
n
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(group = group),
  design = ~ group
)
dds
group <- factor(data$sample_type)
levels(group)
group <- relevel(group, ref = "Solid Tissue Normal")
levels(group)
keep <- group %in% c("Solid Tissue Normal", "Primary Tumor")
table(group[keep])
dim(counts_filtered)
counts_filtered <- counts[, keep]
dim(counts_filtered)
keep_genes <- rowSums(counts_filtered >= 10) >= 10
table(keep_genes)
sum(keep_genes)
counts_filtered <- counts_filtered[keep_genes, ]
dim(counts_filtered)
col_data <- data.frame(
  group = group[keep]
)

rownames(col_data) <- colnames(counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = col_data,
  design = ~ group
)
col_data <- data.frame(
  group = group[keep]
)
rownames(col_data) <- colnames(counts_filtered)
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = col_data,
  design = ~ group
)
table(dds$group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Primary Tumor", "Solid Tissue Normal"))
res
summary(res)
sig <- res[
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) >= 1,
]
nrow(sig)
upregulated <- sig[sig$log2FoldChange > 1, ]
downregulated <- sig[sig$log2FoldChange < -1, ]
nrow(upregulated)
nrow(downregulated)
head(
  upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ],
  20
)
head(
  downregulated[order(downregulated$log2FoldChange), ],
  20
)
range(upregulated$log2FoldChange)
range(downregulated$log2FoldChange)
top_genes <- head(
  res[order(res$padj), ],
  20
)
top_genes
library(org.Hs.eg.db)
BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
library(org.Hs.eg.db)
res$ensembl_id <- sub("\\..*", "", rownames(res))
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = res$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res$gene_symbol <- gene_symbols
head(
  res[order(res$padj), 
      c("gene_symbol", "baseMean", "log2FoldChange", "padj")],
  20
)
upregulated$gene_symbol <- res$gene_symbol[
  match(rownames(upregulated), rownames(res))
]
head(
  upregulated[
    order(upregulated$log2FoldChange, decreasing = TRUE),
    c("gene_symbol", "baseMean", "log2FoldChange", "padj")
  ],
  20
)
downregulated$gene_symbol <- res$gene_symbol[
  match(rownames(downregulated), rownames(res))
]
head(
  downregulated[
    order(downregulated$log2FoldChange),
    c("gene_symbol", "baseMean", "log2FoldChange", "padj")
  ],
  20
)
summary(res$log2FoldChange)
summary(-log10(res$padj))
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 20,
  cex = 0.5,
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "TCGA Prostate Cancer: Tumor vs Normal"
)
res$significant <- !is.na(res$padj) &
  res$padj < 0.05 &
  abs(res$log2FoldChange) >= 1
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 20,
  cex = 0.5,
  col = ifelse(res$significant, "red", "grey"),
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Differential Gene Expression: Prostate Cancer vs Normal"
)

abline(
  v = c(-1, 1),
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)
res$significant <- !is.na(res$padj) &
  res$padj < 0.05 &
  abs(res$log2FoldChange) >= 1
table(res$significant)
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 20,
  cex = 0.5,
  col = ifelse(res$significant, "red", "grey"),
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Differential Gene Expression: Prostate Cancer vs Normal"
)

abline(
  v = c(-1, 1),
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 20,
  cex = 0.5,
  col = ifelse(res$significant, "red", "grey"),
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Differential Gene Expression: Prostate Cancer vs Normal"
)

abline(
  v = c(-1, 1),
  lty = 2
)
up_genes <- rownames(res)[
  !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange >= 1
]

down_genes <- rownames(res)[
  !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange <= -1
]

length(up_genes)
length(down_genes)
# Remove Ensembl version numbers
up_ensembl <- sub("\\..*$", "", up_genes)
down_ensembl <- sub("\\..*$", "", down_genes)

# Check the first few IDs
head(up_ensembl)
head(down_ensembl)
library(org.Hs.eg.db)
library(AnnotationDbi)

up_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = up_ensembl,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL"
)

down_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = down_ensembl,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL"
)
head(up_map)
up_entrez <- unique(
  na.omit(up_map$ENTREZID)
)

down_entrez <- unique(
  na.omit(down_map$ENTREZID)
)

length(up_entrez)
length(down_entrez)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

library(clusterProfiler)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
library(clusterProfiler)
ego_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(ego_up)
ego_up_df <- as.data.frame(ego_up)
head(
  ego_up_df[
    order(ego_up_df$p.adjust),
    c("Description", "GeneRatio", "Count", "p.adjust")
  ],
  20
)
ego_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
ego_down_df <- as.data.frame(ego_down)
head(
  ego_down_df[
    order(ego_down_df$p.adjust),
    c("Description", "GeneRatio", "Count", "p.adjust")
  ],
  20
)
ego_down_df[
  ego_down_df$Description == "muscle contraction",
  c("Description", "Count", "p.adjust", "geneID")
]
dotplot(
  ego_up,
  showCategory = 15,
  title = "GO Biological Processes - Upregulated Genes"
)
dotplot(
  ego_up,
  showCategory = 10,
  title = "GO Biological Processes - Upregulated Genes"
)
dotplot(
  ego_down,
  showCategory = 10,
  title = "GO Biological Processes - Downregulated Genes"
)
# Remove Ensembl version numbers
ensembl_ids <- sub("\\..*$", "", rownames(res))

# Convert Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Create ranked vector using DESeq2 Wald statistic
gene_list <- res$stat

names(gene_list) <- entrez_ids

# Remove genes without Entrez IDs
gene_list <- gene_list[!is.na(names(gene_list))]

# Remove duplicated Entrez IDs
gene_list <- gene_list[!duplicated(names(gene_list))]

# Sort from most upregulated to most downregulated
gene_list <- sort(gene_list, decreasing = TRUE)

# Check
head(gene_list)
tail(gene_list)
length(gene_list)
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0
)
library(BiocParallel)

register(SerialParam())
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0
)
gsea_df <- as.data.frame(gsea_go)

head(
  gsea_df[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")],
  20
)
sum(gsea_df$p.adjust < 0.05, na.rm = TRUE)
dotplot(
  gsea_go,
  showCategory = 20,
  title = "GSEA - GO Biological Processes"
)
library(ggplot2)

gsea_up_plot <- head(
  gsea_up[order(gsea_up$p.adjust), ],
  15
)

ggplot(
  gsea_up_plot,
  aes(
    x = NES,
    y = reorder(Description, NES),
    size = setSize,
    color = -log10(p.adjust)
  )
) +
  geom_point() +
  labs(
    title = "GSEA - Pathways Enriched in Primary Tumor",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "Gene Set Size",
    color = "-log10 Adjusted P-value"
  ) +
  theme_minimal()
gsea_up <- gsea_df[gsea_df$NES > 0 & gsea_df$p.adjust < 0.05, ]

gsea_down <- gsea_df[gsea_df$NES < 0 & gsea_df$p.adjust < 0.05, ]

nrow(gsea_up)
nrow(gsea_down)
gsea_up_plot <- head(
  gsea_up[order(gsea_up$p.adjust), ],
  15
)

gsea_down_plot <- head(
  gsea_down[order(gsea_down$p.adjust), ],
  15
)
gsea_up_plot[, c("ID", "Description", "NES", "p.adjust")]
gsea_down_plot[, c("ID", "Description", "NES", "p.adjust")]
library(ggplot2)

ggplot(gsea_up_plot,
       aes(x = NES,
           y = reorder(Description, NES),
           size = -log10(p.adjust),
           color = NES)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "GSEA - Upregulated Biological Processes",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "-log10 adjusted p-value"
  )
ggplot(gsea_down_plot,
       aes(x = NES,
           y = reorder(Description, NES),
           size = -log10(p.adjust),
           color = NES)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "GSEA - Downregulated Biological Processes",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = "-log10 adjusted p-value"
  )
library(clusterProfiler)

gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_kegg_df <- as.data.frame(gsea_kegg)

head(
  gsea_kegg_df[
    order(gsea_kegg_df$p.adjust),
    c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")
  ],
  20
)
kegg_up <- gsea_kegg_df[
  gsea_kegg_df$NES > 0 & gsea_kegg_df$p.adjust < 0.05,
]

kegg_down <- gsea_kegg_df[
  gsea_kegg_df$NES < 0 & gsea_kegg_df$p.adjust < 0.05,
]
head(
  kegg_up[
    order(kegg_up$p.adjust),
    c("ID", "Description", "setSize", "NES", "p.adjust")
  ],
  15
)
head(
  kegg_down[
    order(kegg_down$p.adjust),
    c("ID", "Description", "setSize", "NES", "p.adjust")
  ],
  15
)
ribosome_genes <- gsea_kegg@result[
  gsea_kegg@result$ID == "hsa03010",
  "core_enrichment"
]

ribosome_genes
ribosome_genes <- unlist(strsplit(ribosome_genes, "/"))

length(ribosome_genes)
head(ribosome_genes)
ribosome_symbols <- bitr(
  ribosome_genes,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

head(ribosome_symbols, 30)
ribosome_df <- data.frame(
  ENTREZID = names(gene_list),
  statistic = as.numeric(gene_list)
)

ribosome_df <- merge(
  ribosome_df,
  ribosome_symbols,
  by = "ENTREZID"
)

ribosome_df <- ribosome_df[
  order(-abs(ribosome_df$statistic)),
]

head(ribosome_df, 30)
focal_genes <- gsea_kegg@result[
  gsea_kegg@result$ID == "hsa04510",
  "core_enrichment"
]

focal_genes
focal_genes <- unlist(strsplit(focal_genes, "/"))

length(focal_genes)
head(focal_genes)
focal_symbols <- bitr(
  focal_genes,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)
focal_df <- data.frame(
  ENTREZID = names(gene_list),
  statistic = as.numeric(gene_list)
)

focal_df <- merge(
  focal_df,
  focal_symbols,
  by = "ENTREZID"
)

focal_df <- focal_df[
  order(focal_df$statistic),
]

head(focal_df, 30)
library(enrichplot)
library(ggplot2)

# Top upregulated pathway
gseaplot2(
  gsea_go,
  geneSetID = "GO:0022613",
  title = "Ribonucleoprotein complex biogenesis"
)
gseaplot2(
  gsea_go,
  geneSetID = "GO:0006936",
  title = "Muscle contraction"
)
gseaplot2(
  kegg_go,
  geneSetID = "hsa03010",
  title = "Ribosome"
)
gseaplot2(
  gsea_kegg,
  geneSetID = "hsa03010",
  title = "KEGG GSEA - Ribosome"
)
gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04971",
  title = "KEGG GSEA - Gastric Acid Secretion"
)
gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04510",
  title = "KEGG GSEA - Focal Adhesion"
)
gsea_kegg@result[
  gsea_kegg@result$ID == "hsa03010",
  c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
gsea_kegg@result[
  gsea_kegg@result$ID == "hsa04510",
  c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
gsea_kegg@result[
  gsea_kegg@result$ID == "hsa04971",
  c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
library(AnnotationDbi)
library(org.Hs.eg.db)

ribosome_genes <- strsplit(
  gsea_kegg@result[
    gsea_kegg@result$ID == "hsa03010",
    "core_enrichment"
  ],
  "/"
)[[1]]

ribosome_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ribosome_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

head(ribosome_symbols, 30)
focal_genes <- strsplit(
  gsea_kegg@result[
    gsea_kegg@result$ID == "hsa04510",
    "core_enrichment"
  ],
  "/"
)[[1]]

focal_symbols <- mapIds(
  org.Hs.eg.db,
  keys = focal_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

head(focal_symbols, 30)
gastric_genes <- strsplit(
  gsea_kegg@result[
    gsea_kegg@result$ID == "hsa04971",
    "core_enrichment"
  ],
  "/"
)[[1]]

gastric_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gastric_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

head(gastric_symbols, 30)
library(enrichplot)
library(ggplot2)

dotplot(
  gsea_kegg,
  showCategory = 15,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  ggtitle("KEGG Gene Set Enrichment Analysis")
library(ggplot2)
library(dplyr)

# Select the top 10 upregulated and top 10 downregulated pathways
kegg_plot <- gsea_kegg@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10) %>%
  bind_rows(
    gsea_kegg@result %>%
      filter(p.adjust < 0.05) %>%
      arrange(NES) %>%
      slice_head(n = 10)
  ) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    Description = factor(
      Description,
      levels = Description[order(NES)]
    )
  )

ggplot(kegg_plot,
       aes(x = NES,
           y = Description,
           size = setSize,
           color = p.adjust)) +
  
  geom_point(alpha = 0.9) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  scale_color_continuous(
    trans = "reverse",
    name = "Adjusted\np-value"
  ) +
  
  scale_size_continuous(
    name = "Gene set size",
    range = c(3, 9)
  ) +
  
  labs(
    title = "KEGG Gene Set Enrichment Analysis",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.text.x = element_text(
      size = 10
    ),
    legend.title = element_text(
      face = "bold"
    ),
    panel.grid.major.y = element_line(
      linewidth = 0.2
    )
  )
ggplot(
  kegg_plot,
  aes(
    x = NES,
    y = Description,
    size = setSize,
    color = NES
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(alpha = 0.9) +
  
  scale_size_continuous(
    name = "Gene set size",
    range = c(3, 8)
  ) +
  
  scale_color_gradient2(
    name = "NES",
    midpoint = 0
  ) +
  
  labs(
    title = "KEGG Gene Set Enrichment Analysis",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    plot.margin = margin(
      10, 20, 10, 10
    )
  )
library(org.Hs.eg.db)
library(AnnotationDbi)

convert_to_symbols <- function(entrez_ids) {
  
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  data.frame(
    ENTREZID = entrez_ids,
    SYMBOL = unname(symbols),
    stringsAsFactors = FALSE
  )
}

ribosome_df <- convert_to_symbols(ribosome_genes)
focal_df <- convert_to_symbols(focal_genes)
gastric_df <- convert_to_symbols(gastric_genes)
calcium_df <- convert_to_symbols(calcium_genes)
calcium_genes <- get_core_genes("hsa04020")

length(calcium_genes)
calcium_genes <- get_core_genes("hsa04020")
length(calcium_genes)
calcium_genes <- strsplit(
  gsea_kegg@result$core_enrichment[
    gsea_kegg@result$ID == "hsa04020"
  ],
  "/"
)[[1]]

length(calcium_genes)
calcium_df <- convert_to_symbols(calcium_genes)

head(calcium_df, 30)
calcium_ranked <- data.frame(
  ENTREZID = calcium_genes,
  statistic = gene_list[calcium_genes]
)

calcium_ranked$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = calcium_ranked$ENTREZID,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

calcium_ranked <- calcium_ranked[
  order(calcium_ranked$statistic),
]

head(calcium_ranked, 20)
# Convert the leading-edge gene lists to symbols

focal_symbols_all <- mapIds(
  org.Hs.eg.db,
  keys = focal_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

calcium_symbols_all <- mapIds(
  org.Hs.eg.db,
  keys = calcium_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

ribosome_symbols_all <- mapIds(
  org.Hs.eg.db,
  keys = ribosome_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

# Find genes shared by focal adhesion and calcium signaling

focal_calcium_overlap <- intersect(
  na.omit(focal_symbols_all),
  na.omit(calcium_symbols_all)
)

focal_calcium_overlap
overlap_stats <- data.frame(
  SYMBOL = focal_calcium_overlap
)

overlap_stats$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = overlap_stats$SYMBOL,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

overlap_stats$statistic <- gene_list[
  as.character(overlap_stats$ENTREZID)
]

overlap_stats <- overlap_stats[
  order(overlap_stats$statistic),
]

overlap_stats
focal_calcium_ribosome <- Reduce(
  intersect,
  list(
    na.omit(focal_symbols_all),
    na.omit(calcium_symbols_all),
    na.omit(ribosome_symbols_all)
  )
)

focal_calcium_ribosome
muscle_genes <- strsplit(
  gsea_kegg@result$core_enrichment[
    gsea_kegg@result$ID == "hsa04820"
  ],
  "/"
)[[1]]

muscle_symbols_all <- mapIds(
  org.Hs.eg.db,
  keys = muscle_genes,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)
down_overlap <- Reduce(
  intersect,
  list(
    na.omit(focal_symbols_all),
    na.omit(calcium_symbols_all),
    na.omit(muscle_symbols_all)
  )
)

down_overlap
ribosome_ranked <- data.frame(
  ENTREZID = ribosome_genes,
  statistic = gene_list[ribosome_genes]
)

ribosome_ranked$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = ribosome_ranked$ENTREZID,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

ribosome_ranked <- ribosome_ranked[
  order(-ribosome_ranked$statistic),
]

head(ribosome_ranked, 20)
down_candidates <- merge(
  focal_ranked[, c("ENTREZID", "statistic", "SYMBOL")],
  calcium_ranked[, c("ENTREZID", "statistic")],
  by = "ENTREZID",
  suffixes = c("_focal", "_calcium")
)

down_candidates$mean_statistic <- rowMeans(
  down_candidates[, c("statistic_focal", "statistic_calcium")],
  na.rm = TRUE
)

down_candidates <- down_candidates[
  order(down_candidates$mean_statistic),
]

head(down_candidates, 20)
focal_ranked <- data.frame(
  ENTREZID = focal_genes,
  statistic = gene_list[focal_genes]
)

focal_ranked$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = focal_ranked$ENTREZID,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

focal_ranked <- focal_ranked[
  order(focal_ranked$statistic),
]
head(focal_ranked, 20)
down_candidates <- merge(
  focal_ranked[, c("ENTREZID", "statistic", "SYMBOL")],
  calcium_ranked[, c("ENTREZID", "statistic")],
  by = "ENTREZID",
  suffixes = c("_focal", "_calcium")
)

down_candidates$mean_statistic <- rowMeans(
  down_candidates[, c("statistic_focal", "statistic_calcium")],
  na.rm = TRUE
)

down_candidates <- down_candidates[
  order(down_candidates$mean_statistic),
]

down_candidates
colnames(colData(dds))
head(as.data.frame(colData(dds)))
save.image("TCGA_PRAD_GSEA_progress.RData")
file.exists("TCGA_PRAD_GSEA_progress.RData")
head(colnames(dds), 10)
sample_ids <- colnames(dds)

patient_ids <- substr(
  sample_ids,
  1,
  12
)

head(patient_ids, 10)
length(unique(patient_ids))
table(table(patient_ids))
patient_counts <- table(patient_ids)

patient_counts[patient_counts > 1]
sample_info <- data.frame(
  sample_id = colnames(dds),
  patient_id = patient_ids,
  group = as.character(colData(dds)$group)
)

sample_info[
  sample_info$patient_id %in% names(patient_counts[patient_counts == 3]),
]
tumor_samples <- sample_info[
  sample_info$group == "Primary Tumor",
]

tumor_counts <- table(tumor_samples$patient_id)

table(tumor_counts)
tumor_counts[tumor_counts > 1]
tumor_samples_unique <- tumor_samples[
  !duplicated(tumor_samples$patient_id),
]

dim(tumor_samples_unique)
length(unique(tumor_samples_unique$patient_id))
ls()
save.image("TCGA_PRAD_GSEA_progress.RData")
file.exists("TCGA_PRAD_GSEA_progress.RData")
"TCGAbiolinks" %in% rownames(installed.packages())
library(TCGAbiolinks)
clinical <- GDCquery_clinic(
  project = "TCGA-PRAD",
  type = "clinical"
)
dim(clinical)
colnames(clinical)
clinical_check <- clinical[, c(
  "submitter_id",
  "age_at_diagnosis",
  "gleason_score",
  "ajcc_pathologic_stage",
  "ajcc_pathologic_t",
  "ajcc_pathologic_n",
  "ajcc_pathologic_m",
  "progression_or_recurrence",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up"
)]

colSums(is.na(clinical_check))
table(clinical$gleason_score, useNA = "ifany")
table(clinical$ajcc_pathologic_stage, useNA = "ifany")
head(clinical$submitter_id)
head(tumor_samples_unique$patient_id)
clinical_prad <- clinical[
  clinical$submitter_id %in% tumor_samples_unique$patient_id,
  c(
    "submitter_id",
    "age_at_diagnosis",
    "gleason_score",
    "ajcc_pathologic_t",
    "ajcc_pathologic_n",
    "vital_status",
    "days_to_last_follow_up"
  )
]
nrow(clinical_prad)
sum(
  tumor_samples_unique$patient_id %in% clinical_prad$submitter_id
)
tumor_patient_ids <- tumor_samples_unique$patient_id
tumor_sample_ids <- tumor_samples_unique$sample_id
head(tumor_patient_ids)
head(tumor_sample_ids)
dim(counts_filtered)
vsd <- vst(dds, blind = FALSE)
vsd
vsd_tumor <- assay(vsd)[, tumor_sample_ids]
dim(vsd_tumor)
all(colnames(vsd_tumor) == tumor_sample_ids)
candidate_genes <- c(
  "MYLK",
  "PRKCB",
  "PDGFD",
  "PRKCA",
  "PDGFC",
  "VEGFD",
  "MET",
  "MYLK4",
  "KDR",
  "EGFR",
  "PDGFRA"
)
candidate_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = candidate_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
)

candidate_map
candidate_ensembl <- candidate_map$ENSEMBL

candidate_expr <- vsd_tumor[
  rownames(vsd_tumor) %in% candidate_ensembl,
  ,
  drop = FALSE
]

rownames(candidate_expr)
candidate_expr <- candidate_expr[
  match(candidate_ensembl, rownames(candidate_expr)),
  ,
  drop = FALSE
]

rownames(candidate_expr) <- candidate_map$SYMBOL

dim(candidate_expr)
head(candidate_expr)
# Remove Ensembl version numbers from VSD rownames
vsd_ensembl <- sub("\\..*$", "", rownames(vsd_tumor))

# Find the 11 candidate genes
candidate_idx <- match(candidate_map$ENSEMBL, vsd_ensembl)

# Check that all 11 were found
candidate_idx
sum(is.na(candidate_idx))
candidate_expr <- vsd_tumor[candidate_idx, , drop = FALSE]

rownames(candidate_expr) <- candidate_map$SYMBOL

dim(candidate_expr)
head(candidate_expr)
sum(is.na(candidate_idx))
rownames(candidate_expr)
# Create clinical data for the 497 tumor patients
gleason_df <- clinical_prad[, c("submitter_id", "gleason_score")]

# Match clinical data to the tumor samples
gleason_df <- gleason_df[
  match(tumor_patient_ids, gleason_df$submitter_id),
]

# Check that the order matches
all(gleason_df$submitter_id == tumor_patient_ids)

# Add expression for each candidate gene
gleason_expr <- t(candidate_expr)

# Make sure samples are in the same order
all(rownames(gleason_expr) == tumor_sample_ids)

# Combine
gleason_analysis <- data.frame(
  patient_id = tumor_patient_ids,
  sample_id = tumor_sample_ids,
  gleason_score = gleason_df$gleason_score,
  gleason_expr
)

head(gleason_analysis)
table(gleason_analysis$gleason_score, useNA = "ifany")
sum(is.na(gleason_analysis$gleason_score))
cor_results <- data.frame(
  Gene = rownames(candidate_expr),
  Spearman_rho = NA_real_,
  p_value = NA_real_
)

for (i in seq_len(nrow(candidate_expr))) {
  
  test <- cor.test(
    candidate_expr[i, ],
    gleason_analysis$gleason_score,
    method = "spearman",
    exact = FALSE
  )
  
  cor_results$Spearman_rho[i] <- unname(test$estimate)
  cor_results$p_value[i] <- test$p.value
}

cor_results$p_adjusted <- p.adjust(
  cor_results$p_value,
  method = "BH"
)

cor_results <- cor_results[
  order(cor_results$p_adjusted),
]

cor_results
cor_results[, c(
  "Gene",
  "Spearman_rho",
  "p_value",
  "p_adjusted"
)]
library(ggplot2)

ggplot(
  gleason_analysis,
  aes(
    x = factor(gleason_score),
    y = MYLK
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.35) +
  labs(
    x = "Gleason Score",
    y = "VST Expression",
    title = "MYLK Expression Across Gleason Scores"
  ) +
  theme_classic()
gleason_analysis$age_at_diagnosis <- clinical_prad$age_at_diagnosis[
  match(gleason_analysis$patient_id, clinical_prad$submitter_id)
]

summary(gleason_analysis$age_at_diagnosis)
gleason_analysis$age_years <- gleason_analysis$age_at_diagnosis / 365.25

summary(gleason_analysis$age_years)
mylk_model <- lm(
  MYLK ~ gleason_score + age_years,
  data = gleason_analysis
)

summary(mylk_model)
sig_genes <- c(
  "MYLK",
  "PRKCB",
  "MYLK4",
  "PDGFC",
  "VEGFD"
)

adjusted_results <- data.frame(
  Gene = sig_genes,
  Beta_Gleason = NA_real_,
  p_value = NA_real_
)

for (i in seq_along(sig_genes)) {
  
  gene <- sig_genes[i]
  
  model <- lm(
    as.formula(paste(gene, "~ gleason_score + age_years")),
    data = gleason_analysis
  )
  
  adjusted_results$Beta_Gleason[i] <-
    coef(summary(model))["gleason_score", "Estimate"]
  
  adjusted_results$p_value[i] <-
    coef(summary(model))["gleason_score", "Pr(>|t|)"]
}

adjusted_results$p_adjusted <- p.adjust(
  adjusted_results$p_value,
  method = "BH"
)

adjusted_results[
  order(adjusted_results$p_adjusted),
]
gleason_analysis$gleason_group <- ifelse(
  gleason_analysis$gleason_score <= 7,
  "Gleason 6-7",
  "Gleason 8-10"
)

table(gleason_analysis$gleason_group)
wilcox.test(
  MYLK ~ gleason_group,
  data = gleason_analysis
)
group_results <- data.frame(
  Gene = sig_genes,
  p_value = NA_real_
)

for (i in seq_along(sig_genes)) {
  
  gene <- sig_genes[i]
  
  test <- wilcox.test(
    as.formula(paste(gene, "~ gleason_group")),
    data = gleason_analysis
  )
  
  group_results$p_value[i] <- test$p.value
}

group_results$p_adjusted <- p.adjust(
  group_results$p_value,
  method = "BH"
)

group_results[
  order(group_results$p_adjusted),
]
library(ggplot2)
library(tidyr)
library(dplyr)

plot_df <- gleason_analysis %>%
  select(
    gleason_group,
    MYLK,
    PRKCB,
    PDGFC,
    MYLK4,
    VEGFD
  ) %>%
  pivot_longer(
    cols = c(MYLK, PRKCB, PDGFC, MYLK4, VEGFD),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_df$Gene <- factor(
  plot_df$Gene,
  levels = c("MYLK", "PRKCB", "PDGFC", "MYLK4", "VEGFD")
)

ggplot(
  plot_df,
  aes(
    x = gleason_group,
    y = Expression
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.25,
    size = 1
  ) +
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title = "Candidate Gene Expression by Gleason Grade",
    x = NULL,
    y = "VST Expression"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
    
    
    
    
    plot_df <- gleason_analysis %>%plot_df <- gleason_analysis %>%
      select(
        gleason_group,
        MYLK,
        PRKCB,
        PDGFC,
        MYLK4,
        VEGFD
      ) %>%
      pivot_longer(
        cols = c(MYLK, PRKCB, PDGFC, MYLK4, VEGFD),
        names_to = "Gene",
        values_to = "Expression"
      )
    head(plot_df)
    ggplot(
      plot_df,
      aes(
        x = gleason_group,
        y = Expression
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.55
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.25,
        size = 1.2
      ) +
      facet_wrap(
        ~ Gene,
        scales = "free_y",
        ncol = 3
      ) +
      labs(
        title = "Candidate Gene Expression by Gleason Grade",
        x = NULL,
        y = "VST Expression"
      ) +
      theme_classic(base_size = 13) +
      theme(
        strip.text = element_text(face = "bold"),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.text.x = element_text(
          angle = 0,
          size = 10
        )
      )
    ggplot(
      plot_df,
      aes(
        x = gleason_group,
        y = Expression
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.55
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.20,
        size = 1.1
      ) +
      facet_wrap(
        ~ Gene,
        scales = "free_y",
        ncol = 3
      ) +
      scale_x_discrete(
        labels = c(
          "Gleason 6-7" = "6–7",
          "Gleason 8-10" = "8–10"
        )
      ) +
      labs(
        title = "Candidate Gene Expression by Gleason Grade",
        x = "Gleason group",
        y = "VST expression"
      ) +
      theme_classic(base_size = 13) +
      theme(
        strip.text = element_text(
          face = "bold",
          size = 12
        ),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 16
        ),
        axis.text.x = element_text(
          size = 10
        ),
        axis.title = element_text(
          size = 12
        )
      )
    plot_df$Gene <- factor(
      plot_df$Gene,
      levels = c("MYLK", "PRKCB", "VEGFD", "PDGFC", "MYLK4")
    )
    ggplot(
      plot_df,
      aes(
        x = gleason_group,
        y = Expression
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.55
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.20,
        size = 1.1
      ) +
      facet_wrap(
        ~ Gene,
        scales = "free_y",
        ncol = 3
      ) +
      scale_x_discrete(
        labels = c(
          "Gleason 6-7" = "6–7",
          "Gleason 8-10" = "8–10"
        )
      ) +
      labs(
        title = "Candidate Gene Expression by Gleason Grade",
        x = "Gleason group",
        y = "VST expression"
      ) +
      theme_classic(base_size = 13) +
      theme(
        strip.text = element_text(
          face = "bold",
          size = 12
        ),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 16
        ),
        axis.text.x = element_text(
          size = 10
        ),
        axis.title = element_text(
          size = 12
        )
      )
    gleason_model_df <- gleason_analysis %>%
      mutate(
        age_years = age_at_diagnosis / 365.25
      ) %>%
      select(
        gleason_score,
        age_years,
        MYLK,
        PRKCB,
        VEGFD,
        PDGFC,
        MYLK4
      ) %>%
      na.omit()
    
    dim(gleason_model_df)
    multivariable_model <- lm(
      gleason_score ~ age_years + MYLK + PRKCB + VEGFD + PDGFC + MYLK4,
      data = gleason_model_df
    )
    
    summary(multivariable_model)
    coef_results <- as.data.frame(summary(multivariable_model)$coefficients)
    
    coef_results
    coef_results[, c("Estimate", "Pr(>|t|)")]
    cor(
      gleason_model_df[, c("MYLK", "PRKCB", "VEGFD", "PDGFC", "MYLK4")],
      method = "spearman"
    )  
    library(car)
    
    vif(multivariable_model)
    predictors <- c("age_years", "MYLK", "PRKCB", "VEGFD", "PDGFC", "MYLK4")
    
    vif_values <- sapply(predictors, function(x) {
      
      other_predictors <- setdiff(predictors, x)
      
      formula_vif <- as.formula(
        paste(x, "~", paste(other_predictors, collapse = " + "))
      )
      
      r2 <- summary(lm(formula_vif, data = gleason_model_df))$r.squared
      
      1 / (1 - r2)
    })
    
    vif_values
    ggplot(
      gleason_analysis,
      aes(
        x = gleason_score,
        y = VEGFD
      )
    ) +
      geom_jitter(
        width = 0.15,
        alpha = 0.25,
        size = 1
      ) +
      geom_smooth(
        method = "lm",
        se = TRUE
      ) +
      labs(
        title = "VEGFD Expression and Gleason Score",
        x = "Gleason Score",
        y = "VST Expression"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        )
      )
    table(clinical_prad$vital_status, useNA = "ifany")
    
    summary(clinical_prad$days_to_death)
    
    summary(clinical_prad$days_to_last_follow_up)
    
    "survival" %in% rownames(installed.packages())
    head(clinical_prad[, c(
      "submitter_id",
      "vital_status",
      "days_to_death",
      "days_to_last_follow_up"
    )])
    grep(
      "death|follow|vital",
      colnames(clinical_prad),
      ignore.case = TRUE,
      value = TRUE
    )
    colnames(clinical_prad)
    clinical_prad[clinical_prad$vital_status == "Dead", ]    
    clinical[
      clinical$submitter_id %in% clinical_prad$submitter_id &
        clinical$vital_status == "Dead",
      c(
        "submitter_id",
        "vital_status",
        "days_to_death",
        "days_to_last_follow_up"
      )
    ] 
    survival_clinical <- clinical_prad %>%
      mutate(
        survival_time = ifelse(
          vital_status == "Dead",
          clinical$days_to_death[match(submitter_id, clinical$submitter_id)],
          days_to_last_follow_up
        ),
        event = ifelse(vital_status == "Dead", 1, 0)
      )
    table(survival_clinical$event)
    
    summary(survival_clinical$survival_time)
    
    sum(is.na(survival_clinical$survival_time))
    vegfd_expr <- candidate_expr["VEGFD", ]
    
    length(vegfd_expr)
    vegfd_df <- data.frame(
      sample_id = colnames(candidate_expr),
      VEGFD = as.numeric(vegfd_expr),
      patient_id = substr(colnames(candidate_expr), 1, 12)
    )
    head(vegfd_df)
    nrow(vegfd_df)
    length(unique(vegfd_df$patient_id))
    survival_df <- survival_clinical %>%
      select(
        submitter_id,
        survival_time,
        event
      ) %>%
      left_join(
        vegfd_df %>%
          select(patient_id, VEGFD),
        by = c("submitter_id" = "patient_id")
      )
    dim(survival_df)
    
    sum(is.na(survival_df$VEGFD))
    
    table(survival_df$event)
    library(survival)
    library(survminer)
    survival_df$VEGFD_group <- ifelse(
      survival_df$VEGFD >= median(survival_df$VEGFD),
      "High VEGFD",
      "Low VEGFD"
    )
    
    table(survival_df$VEGFD_group)
    surv_object <- Surv(
      time = survival_df$survival_time,
      event = survival_df$event
    )
    km_vegfd <- survfit(
      surv_object ~ VEGFD_group,
      data = survival_df
    )
    
    summary(km_vegfd)
    logrank_vegfd <- survdiff(
      surv_object ~ VEGFD_group,
      data = survival_df
    )
    
    logrank_vegfd
    p_logrank <- 1 - pchisq(
      logrank_vegfd$chisq,
      df = length(logrank_vegfd$n) - 1
    )
    
    p_logrank
    table(
      gleason_analysis$ajcc_pathologic_t,
      useNA = "ifany"
    )
    table(
      gleason_analysis$ajcc_pathologic_n,
      useNA = "ifany"
    )
    table(
      gleason_analysis$gleason_score,
      useNA = "ifany"
    )
    ordinal_df <- gleason_analysis %>%
      select(
        submitter_id,
        gleason_score,
        age_at_diagnosis,
        VEGFD
      ) %>%
      left_join(
        clinical_prad %>%
          select(
            submitter_id,
            ajcc_pathologic_t,
            ajcc_pathologic_n
          ),
        by = "submitter_id"
      )
    colnames(gleason_analysis)
    head(gleason_analysis)
    ordinal_df <- gleason_analysis %>%
      select(
        patient_id,
        gleason_score,
        age_years,
        VEGFD
      ) %>%
      left_join(
        clinical_prad %>%
          select(
            submitter_id,
            ajcc_pathologic_t,
            ajcc_pathologic_n
          ),
        by = c("patient_id" = "submitter_id")
      )    
    dim(ordinal_df)
    
    table(ordinal_df$ajcc_pathologic_t, useNA = "ifany")
    
    table(ordinal_df$ajcc_pathologic_n, useNA = "ifany")
    
    colSums(is.na(ordinal_df))
    ordinal_df$gleason_score <- ordered(
      ordinal_df$gleason_score,
      levels = c(6, 7, 8, 9, 10)
    )    
    ordinal_t <- ordinal_df %>%
      filter(
        !is.na(age_years),
        !is.na(ajcc_pathologic_t)
      )
    dim(ordinal_t)
    
    table(ordinal_t$gleason_score)
    
    table(ordinal_t$ajcc_pathologic_t)
    "ordinal" %in% rownames(installed.packages())
    "nnet" %in% rownames(installed.packages())
    "MASS" %in% rownames(installed.packages())    
    library(MASS)
    
    ordinal_model_t <- polr(
      gleason_score ~ VEGFD + age_years + ajcc_pathologic_t,
      data = ordinal_t,
      Hess = TRUE
    )
    
    summary(ordinal_model_t)    
    coef_table <- coef(summary(ordinal_model_t))
    
    coef_table
    p_values <- 2 * pnorm(
      abs(coef_table[, "t value"]),
      lower.tail = FALSE
    )
    
    cbind(
      coef_table,
      p_value = p_values
    )
    VEGFD_OR <- exp(coef(ordinal_model_t)["VEGFD"])
    
    VEGFD_OR
    VEGFD_CI <- exp(
      confint(ordinal_model_t, parm = "VEGFD")
    )
    
    VEGFD_CI
    summary(
      lm(
        as.numeric(ajcc_pathologic_t) ~ VEGFD,
        data = ordinal_t
      )
    )
    aggregate(
      VEGFD ~ ajcc_pathologic_t,
      data = ordinal_t,
      FUN = median
    )
    aggregate(
      VEGFD ~ ajcc_pathologic_t,
      data = ordinal_t,
      FUN = mean
    )
    table(ordinal_t$ajcc_pathologic_t)
    kruskal.test(
      VEGFD ~ ajcc_pathologic_t,
      data = ordinal_t
    )   
    ordinal_t$T_stage_numeric <- match(
      ordinal_t$ajcc_pathologic_t,
      c("T2a", "T2b", "T2c", "T3a", "T3b", "T4")
    )
    cor.test(
      ordinal_t$VEGFD,
      ordinal_t$T_stage_numeric,
      method = "spearman"
    )
    t_model_df <- ordinal_df %>%
      filter(
        !is.na(age_years),
        !is.na(ajcc_pathologic_t)
      )
    t_model_df$T_stage <- ordered(
      t_model_df$ajcc_pathologic_t,
      levels = c("T2a", "T2b", "T2c", "T3a", "T3b", "T4")
    )
    table(t_model_df$T_stage)
    
    sum(is.na(t_model_df$VEGFD))
    sum(is.na(t_model_df$gleason_score))
    sum(is.na(t_model_df$age_years))
    library(MASS)
    
    t_stage_model <- polr(
      T_stage ~ VEGFD + gleason_score + age_years,
      data = t_model_df,
      Hess = TRUE
    )
    
    summary(t_stage_model)
    t_coef <- coef(summary(t_stage_model))
    
    t_pvalues <- 2 * pnorm(
      abs(t_coef[, "t value"]),
      lower.tail = FALSE
    )
    
    cbind(
      t_coef,
      p_value = t_pvalues
    )
    VEGFD_T_OR <- exp(coef(t_stage_model)["VEGFD"])
    
    VEGFD_T_OR
    VEGFD_T_CI <- exp(
      confint(t_stage_model, parm = "VEGFD")
    )
    
    VEGFD_T_CI
    "GEOquery" %in% rownames(installed.packages())
    "limma" %in% rownames(installed.packages())    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    BiocManager::install(
      c("GEOquery", "limma"),
      ask = FALSE,
      update = FALSE
    )
    "GEOquery" %in% rownames(installed.packages())
    "limma" %in% rownames(installed.packages())
    library(GEOquery)
    library(limma)
    
    gse70770 <- getGEO(
      "GSE70770",
      GSEMatrix = TRUE,
      AnnotGPL = TRUE
    )
    gse70768 <- getGEO(
      "GSE70768",
      GSEMatrix = TRUE,
      AnnotGPL = TRUE
    )
    length(gse70768)
    lapply(gse70768, dim)  
    gse70768_eset <- gse70768[[1]]
    
    dim(exprs(gse70768_eset))    
    colnames(pData(gse70768_eset))
    head(pData(gse70768_eset), 10)    
    grep(
      "gleason|grade|stage|path|tumor|patient",
      colnames(pData(gse70768_eset)),
      ignore.case = TRUE,
      value = TRUE
    )   
    table(
      pData(gse70768_eset)$`sample type:ch1`,
      useNA = "ifany"
    )
    table(
      pData(gse70768_eset)$`pathology stage:ch1`,
      useNA = "ifany"
    )
    table(
      pData(gse70768_eset)$`tumour gleason:ch1`,
      useNA = "ifany"
    )
    table(
      pData(gse70768_eset)$`clinical stage:ch1`,
      useNA = "ifany"
    )
    table(
      grepl(
        "castrate-resistant|CRPC",
        pData(gse70768_eset)$title,
        ignore.case = TRUE
      )
    )
    head(
      pData(gse70768_eset)[
        ,
        c(
          "title",
          "age at diag:ch1",
          "tumour gleason:ch1",
          "pathology stage:ch1",
          "sample type:ch1"
        )
      ],
      20
    )
    geo_pd <- pData(gse70768_eset)
    
    geo_pd$CRPC <- grepl(
      "castrate-resistant|CRPC",
      geo_pd$title,
      ignore.case = TRUE
    )
    
    geo_validation <- geo_pd[
      geo_pd$`sample type:ch1` == "Tumour" &
        !geo_pd$CRPC &
        !is.na(geo_pd$`pathology stage:ch1`) &
        geo_pd$`pathology stage:ch1` != "UNKNOWN",
    ]
    
    dim(geo_validation)
    table(geo_validation$`pathology stage:ch1`)
    table(geo_validation$`tumour gleason:ch1`)    
    nrow(geo_validation)   
    sum(geo_validation$CRPC)    
    dim(fData(gse70768_eset))    
    colnames(fData(gse70768_eset))   
    vegfd_probes <- fData(gse70768_eset)[
      toupper(fData(gse70768_eset)$`Gene symbol`) == "VEGFD",
      c("ID", "Gene symbol", "Gene title")
    ]
    
    vegfd_probes    
    nrow(vegfd_probes)
    vegfd_probes    
    vegfd_hits <- which(
      apply(
        fData(gse70768_eset),
        1,
        function(x) any(grepl("VEGFD", x, ignore.case = TRUE))
      )
    )
    
    length(vegfd_hits)   
    fData(gse70768_eset)[
      vegfd_hits,
      c(
        "ID",
        "Gene symbol",
        "Gene title",
        "Gene ID",
        "UniGene symbol",
        "UniGene ID"
      )
    ]
    figf_hits <- which(
      apply(
        fData(gse70768_eset),
        1,
        function(x) any(grepl("FIGF", x, ignore.case = TRUE))
      )
    )
    
    length(figf_hits)
    fData(gse70768_eset)[
      figf_hits,
      c(
        "ID",
        "Gene symbol",
        "Gene title",
        "Gene ID",
        "UniGene symbol",
        "UniGene ID"
      )
    ]
    vegfd_probe <- "ILMN_1707612"
    
    vegfd_geo <- exprs(gse70768_eset)[
      vegfd_probe,
      rownames(geo_validation)
    ]
    
    head(vegfd_geo)
    summary(vegfd_geo)
    all(names(vegfd_geo) == rownames(geo_validation))   
    geo_validation$T_group <- ifelse(
      grepl("^p?T2", geo_validation$`pathology stage:ch1`, ignore.case = TRUE) |
        grepl("^T2", geo_validation$`pathology stage:ch1`, ignore.case = TRUE),
      "T2",
      "T3/T4"
    )    
    table(geo_validation$T_group)
    geo_validation$VEGFD <- as.numeric(vegfd_geo)
    
    wilcox_vegfd <- wilcox.test(
      VEGFD ~ T_group,
      data = geo_validation,
      exact = FALSE
    )
    
    wilcox_vegfd 
    aggregate(
      VEGFD ~ T_group,
      data = geo_validation,
      FUN = median
    )
    aggregate(
      VEGFD ~ T_group,
      data = geo_validation,
      FUN = mean
    )
    tapply(
      geo_validation$VEGFD,
      geo_validation$T_group,
      sd
    )
    candidate_genes <- c(
      "MYLK",
      "PRKCB",
      "PDGFD",
      "PRKCA",
      "PDGFC",
      "VEGFD",
      "MET",
      "MYLK4",
      "KDR",
      "EGFR",
      "PDGFRA"
    )
    candidate_hits <- lapply(
      candidate_genes,
      function(gene) {
        
        synonyms <- switch(
          gene,
          VEGFD = c("VEGFD", "FIGF"),
          gene
        )
        
        idx <- which(
          apply(
            fData(gse70768_eset),
            1,
            function(x) {
              any(
                sapply(
                  synonyms,
                  function(s) grepl(
                    paste0("^", s, "$"),
                    x,
                    ignore.case = TRUE
                  )
                )
              )
            }
          )
        )
        
        data.frame(
          gene = gene,
          probe = fData(gse70768_eset)$ID[idx],
          symbol = fData(gse70768_eset)$`Gene symbol`[idx]
        )
      }
    )
    
    candidate_map_geo <- do.call(rbind, candidate_hits)
    
    candidate_map_geo
    table(candidate_map_geo$gene)
    candidate_probes <- candidate_map_geo$probe
    
    candidate_expr_geo <- exprs(gse70768_eset)[
      candidate_probes,
      rownames(geo_validation)
    ]
    
    dim(candidate_expr_geo)   
    candidate_expr_geo[1:5, 1:5]
    cor(
      t(candidate_expr_geo[
        candidate_map_geo$gene == "MYLK",
        ,
        drop = FALSE
      ]),
      method = "spearman"
    )    
    cor(
      t(candidate_expr_geo[
        candidate_map_geo$gene == "PRKCB",
        ,
        drop = FALSE
      ]),
      method = "spearman"
    )
    cor(
      t(candidate_expr_geo[
        candidate_map_geo$gene == "EGFR",
        ,
        drop = FALSE
      ]),
      method = "spearman"
    )
    multi_genes <- c(
      "MYLK",
      "PRKCB",
      "PDGFD",
      "MYLK4",
      "EGFR",
      "PDGFRA"
    )
    
    candidate_map_geo[
      candidate_map_geo$gene %in% multi_genes,
    ]
    fData(gse70768_eset)[
      candidate_map_geo$probe[
        candidate_map_geo$gene %in% multi_genes
      ],
      c(
        "ID",
        "Gene symbol",
        "Gene title",
        "Gene ID",
        "Chromosome location",
        "Platform_CLONEID",
        "Platform_SEQUENCE"
      )
    ]
    candidate_gene_expr <- sapply(
      candidate_genes,
      function(gene) {
        
        probes <- candidate_map_geo$probe[
          candidate_map_geo$gene == gene
        ]
        
        if (length(probes) == 1) {
          as.numeric(candidate_expr_geo[probes, ])
        } else {
          apply(
            candidate_expr_geo[probes, , drop = FALSE],
            2,
            median,
            na.rm = TRUE
          )
        }
      }
    )
    
    candidate_gene_expr <- t(candidate_gene_expr)
    
    rownames(candidate_gene_expr) <- candidate_genes
    
    dim(candidate_gene_expr)
    head(candidate_gene_expr)
    colnames(candidate_gene_expr)[1:6]    
    sum(is.na(candidate_gene_expr))    
    geo_validation_genes <- geo_validation
    
    geo_validation_genes[
      ,
      candidate_genes
    ] <- t(candidate_gene_expr)    
    dim(geo_validation_genes)
    
    head(
      geo_validation_genes[
        ,
        c("T_group", candidate_genes)
      ]
    )
    validation_results <- lapply(
      candidate_genes,
      function(gene) {
        
        test <- wilcox.test(
          geo_validation_genes[[gene]] ~ geo_validation_genes$T_group,
          exact = FALSE
        )
        
        data.frame(
          Gene = gene,
          p_value = test$p.value,
          T2_median = median(
            geo_validation_genes[
              geo_validation_genes$T_group == "T2",
              gene
            ]
          ),
          T3_T4_median = median(
            geo_validation_genes[
              geo_validation_genes$T_group == "T3/T4",
              gene
            ]
          )
        )
      }
    )
    
    validation_results <- do.call(
      rbind,
      validation_results
    )
    
    validation_results$p_adjusted <- p.adjust(
      validation_results$p_value,
      method = "BH"
    )
    
    validation_results <- validation_results[
      order(validation_results$p_adjusted),
    ]
    
    validation_results
    geo_ordinal$T_stage <- sub(
      ".*?(T[234][a-z]?).*",
      "\\1",
      geo_ordinal$`pathology stage:ch1`
    )
    geo_ordinal <- geo_validation_genes[
      !is.na(geo_validation_genes$`pathology stage:ch1`),
    ]
    geo_ordinal$T_stage <- sub(
      ".*?(T[234][a-z]?).*",
      "\\1",
      geo_ordinal$`pathology stage:ch1`
    )
    table(geo_ordinal$T_stage, useNA = "ifany")
    geo_ordinal$T_stage_simple <- ifelse(
      grepl("^T2", geo_ordinal$T_stage),
      "T2",
      ifelse(
        grepl("^T3", geo_ordinal$T_stage),
        "T3",
        "T4"
      )
    )
    
    geo_ordinal$T_stage_simple <- factor(
      geo_ordinal$T_stage_simple,
      levels = c("T2", "T3", "T4"),
      ordered = TRUE
    )
    
    table(geo_ordinal$T_stage_simple)    
    geo_ordinal$T_advanced <- ifelse(
      geo_ordinal$T_stage_simple == "T2",
      0,
      1
    )
    
    table(geo_ordinal$T_advanced)
    colSums(
      is.na(
        geo_ordinal[
          ,
          c(
            "MET",
            "PDGFRA",
            "age at diag:ch1",
            "tumour gleason:ch1"
          )
        ]
      )
    )
    geo_ordinal$gleason_score <- as.numeric(
      sub(
        "=.*",
        "",
        geo_ordinal$`tumour gleason:ch1`
      )
    )
    
    table(geo_ordinal$gleason_score, useNA = "ifany")
    geo_ordinal$age_years <- as.numeric(
      geo_ordinal$`age at diag:ch1`
    )
    summary(
      geo_ordinal[
        ,
        c(
          "age_years",
          "gleason_score",
          "MET",
          "PDGFRA"
        )
      ]
    )
    geo_model <- glm(
      T_advanced ~ MET + PDGFRA + gleason_score + age_years,
      data = geo_ordinal,
      family = binomial
    )
    
    summary(geo_model)
    geo_coef <- summary(geo_model)$coefficients
    
    geo_OR <- exp(coef(geo_model))
    
    geo_CI <- exp(
      confint.default(geo_model)
    )
    
    geo_results <- data.frame(
      OR = geo_OR,
      CI_lower = geo_CI[, 1],
      CI_upper = geo_CI[, 2],
      p_value = geo_coef[, 4]
    )
    
    geo_results
    cor.test(
      geo_ordinal$MET,
      geo_ordinal$PDGFRA,
      method = "spearman"
    )
    model_MET <- glm(
      T_advanced ~ MET + gleason_score + age_years,
      data = geo_ordinal,
      family = binomial
    )
    
    model_PDGFRA <- glm(
      T_advanced ~ PDGFRA + gleason_score + age_years,
      data = geo_ordinal,
      family = binomial
    )
    
    summary(model_MET)
    summary(model_PDGFRA)
    table(
      geo_ordinal$T_advanced,
      cut(
        geo_ordinal$MET,
        breaks = quantile(
          geo_ordinal$MET,
          probs = c(0, 0.5, 1),
          na.rm = TRUE
        ),
        include.lowest = TRUE
      )
    )
    table(
      geo_ordinal$T_advanced,
      cut(
        geo_ordinal$PDGFRA,
        breaks = quantile(
          geo_ordinal$PDGFRA,
          probs = c(0, 0.5, 1),
          na.rm = TRUE
        ),
        include.lowest = TRUE
      )
    )
    geo_ordinal$T_numeric <- ifelse(
      grepl("^T2", geo_ordinal$T_stage), 2,
      ifelse(
        grepl("^T3a", geo_ordinal$T_stage), 3,
        ifelse(
          grepl("^T3b", geo_ordinal$T_stage), 3,
          ifelse(
            grepl("^T4", geo_ordinal$T_stage), 4,
            NA
          )
        )
      )
    )
    table(geo_ordinal$T_numeric, useNA = "ifany")
    cor.test(
      geo_ordinal$MET,
      geo_ordinal$T_numeric,
      method = "spearman"
    )  
    cor.test(
      geo_ordinal$PDGFRA,
      geo_ordinal$T_numeric,
      method = "spearman"
    )
    table(
      gleason_analysis$ajcc_pathologic_t,
      useNA = "ifany"
    )
    cor.test(
      gleason_analysis$MET,
      as.numeric(gleason_analysis$ajcc_pathologic_t),
      method = "spearman"
    )
    head(clinical_prad[, c(
      "submitter_id",
      "ajcc_pathologic_t",
      "ajcc_pathologic_n"
    )])
    sum(
      gleason_analysis$patient_id %in%
        clinical_prad$submitter_id
    )
    t_stage_map <- clinical_prad[
      ,
      c(
        "submitter_id",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n"
      )
    ]
    
    gleason_analysis_t <- merge(
      gleason_analysis,
      t_stage_map,
      by.x = "patient_id",
      by.y = "submitter_id",
      all.x = TRUE,
      sort = FALSE
    )
    dim(gleason_analysis_t)
    table(
      gleason_analysis_t$ajcc_pathologic_t,
      useNA = "ifany"
    )    
    colSums(
      is.na(
        gleason_analysis_t[
          ,
          c(
            "MET",
            "PDGFRA",
            "gleason_score",
            "age_years",
            "ajcc_pathologic_t"
          )
        ]
      )
    )   
    tcga_t <- gleason_analysis_t[
      !is.na(gleason_analysis_t$ajcc_pathologic_t),
    ]
    tcga_t$T_advanced <- ifelse(
      grepl("^T2", tcga_t$ajcc_pathologic_t),
      0,
      ifelse(
        grepl("^T3|^T4", tcga_t$ajcc_pathologic_t),
        1,
        NA
      )
    )
    table(
      tcga_t$T_advanced,
      useNA = "ifany"
    )
    tcga_t_model <- tcga_t[
      complete.cases(
        tcga_t[
          ,
          c(
            "T_advanced",
            "MET",
            "PDGFRA",
            "gleason_score",
            "age_years"
          )
        ]
      ),
    ]
    dim(tcga_t_model)
    table(tcga_t_model$T_advanced)   
    tcga_model <- glm(
      T_advanced ~ MET + PDGFRA + gleason_score + age_years,
      data = tcga_t_model,
      family = binomial
    )
    
    summary(tcga_model)    
    tcga_coef <- summary(tcga_model)$coefficients
    
    tcga_OR <- exp(coef(tcga_model))
    
    tcga_CI <- exp(
      confint.default(tcga_model)
    )
    
    tcga_results <- data.frame(
      OR = tcga_OR,
      CI_lower = tcga_CI[, 1],
      CI_upper = tcga_CI[, 2],
      p_value = tcga_coef[, 4]
    )
    
    tcga_results
    tcga_MET <- glm(
      T_advanced ~ MET + gleason_score + age_years,
      data = tcga_t_model,
      family = binomial
    )
    
    tcga_PDGFRA <- glm(
      T_advanced ~ PDGFRA + gleason_score + age_years,
      data = tcga_t_model,
      family = binomial
    )
    summary(tcga_MET)
    summary(tcga_PDGFRA)
    table(
      geo_validation$`tumour gleason:ch1`,
      useNA = "ifany"
    )
    summary(
      geo_validation$`age at diag:ch1`
    )
    geo_validation$gleason_score <- as.numeric(
      sub(
        "=.*",
        "",
        geo_validation$`tumour gleason:ch1`
      )
    )
    
    table(
      geo_validation$gleason_score,
      useNA = "ifany"
    )    
    geo_validation$age_years <- as.numeric(
      geo_validation$`age at diag:ch1`
    )
    summary(geo_validation$age_years)
    
    sum(is.na(geo_validation$age_years))
    gleason_validation_genes <- c(
      "MYLK",
      "PDGFC",
      "PRKCB",
      "VEGFD",
      "MYLK4"
    )
    
    gleason_validation_results <- lapply(
      gleason_validation_genes,
      function(gene) {
        
        test <- cor.test(
          geo_validation[[gene]],
          geo_validation$gleason_score,
          method = "spearman",
          exact = FALSE
        )
        
        data.frame(
          Gene = gene,
          Spearman_rho = unname(test$estimate),
          p_value = test$p.value
        )
      }
    )
    
    gleason_validation_results <- do.call(
      rbind,
      gleason_validation_results
    )
    
    gleason_validation_results$p_adjusted <- p.adjust(
      gleason_validation_results$p_value,
      method = "BH"
    )
    
    gleason_validation_results[
      order(gleason_validation_results$p_adjusted),
    ]
    sapply(
      geo_validation[
        ,
        c("MYLK", "PDGFC", "PRKCB", "VEGFD", "MYLK4")
      ],
      class
    )
    colnames(geo_validation)
    sapply(geo_validation, class)   
    genes_to_validate <- c(
      "MYLK",
      "PDGFC",
      "PRKCB",
      "VEGFD",
      "MYLK4"
    )
    
    gene_expr_geo <- sapply(
      genes_to_validate,
      function(gene) {
        
        probes <- candidate_map_geo$probe[
          candidate_map_geo$gene == gene
        ]
        
        apply(
          candidate_expr_geo[probes, , drop = FALSE],
          2,
          median,
          na.rm = TRUE
        )
      }
    )
    
    gene_expr_geo <- as.data.frame(gene_expr_geo)
    
    head(gene_expr_geo)
    dim(gene_expr_geo)    
    all(
      rownames(gene_expr_geo) == rownames(geo_validation)
    )
    geo_gleason <- cbind(
      geo_validation[
        ,
        c(
          "geo_accession",
          "gleason_score",
          "age_years"
        )
      ],
      gene_expr_geo
    )    
    dim(geo_gleason)
    
    colSums(
      is.na(
        geo_gleason[
          ,
          c(
            "gleason_score",
            "age_years",
            genes_to_validate
          )
        ]
      )
    )
    gleason_validation_results <- lapply(
      genes_to_validate,
      function(gene) {
        
        test <- cor.test(
          geo_gleason[[gene]],
          geo_gleason$gleason_score,
          method = "spearman",
          exact = FALSE
        )
        
        data.frame(
          Gene = gene,
          Spearman_rho = unname(test$estimate),
          p_value = test$p.value
        )
      }
    )
    
    gleason_validation_results <- do.call(
      rbind,
      gleason_validation_results
    )
    
    gleason_validation_results$p_adjusted <- p.adjust(
      gleason_validation_results$p_value,
      method = "BH"
    )
    
    gleason_validation_results[
      order(gleason_validation_results$p_adjusted),
    ]
    mylk_model_geo <- lm(
      gleason_score ~ MYLK + age_years,
      data = geo_gleason
    )
    
    summary(mylk_model_geo)
    summary(mylk_model_geo)$coefficients
    mylk_beta <- coef(mylk_model_geo)["MYLK"]
    
    mylk_ci <- confint(
      mylk_model_geo,
      "MYLK"
    )
    
    mylk_beta
    mylk_ci    
    geo_gleason$gleason_ordinal <- factor(
      geo_gleason$gleason_score,
      levels = c(6, 7, 8, 9),
      ordered = TRUE
    )
    table(geo_gleason$gleason_ordinal)
    library(MASS)
    
    mylk_ordinal_geo <- polr(
      gleason_ordinal ~ MYLK + age_years,
      data = geo_gleason,
      Hess = TRUE
    )
    
    summary(mylk_ordinal_geo)    
    mylk_coef <- coef(summary(mylk_ordinal_geo))
    
    mylk_p <- 2 * pnorm(
      abs(mylk_coef["MYLK", "t value"]),
      lower.tail = FALSE
    )
    
    mylk_p
    mylk_OR <- exp(
      coef(mylk_ordinal_geo)["MYLK"]
    )
    
    mylk_CI <- exp(
      confint(
        mylk_ordinal_geo,
        parm = "MYLK"
      )
    )
    
    mylk_OR
    mylk_CI
    aggregate(
      MYLK ~ gleason_score,
      data = geo_gleason,
      FUN = median
    )
    aggregate(
      MYLK ~ gleason_score,
      data = geo_gleason,
      FUN = mean
    )
    tapply(
      geo_gleason$MYLK,
      geo_gleason$gleason_score,
      sd
    )
    tcga_mylk_model <- lm(
      gleason_score ~ MYLK + age_years,
      data = tcga_t_model
    )
    
    summary(tcga_mylk_model)
    summary(tcga_mylk_model)$coefficients
    tcga_mylk_beta <- coef(
      tcga_mylk_model
    )["MYLK"]
    
    tcga_mylk_CI <- confint(
      tcga_mylk_model,
      "MYLK"
    )
    
    tcga_mylk_beta
    tcga_mylk_CI    
    tcga_mylk_full <- lm(
      gleason_score ~ MYLK + PRKCB + PDGFC + VEGFD + MYLK4 + age_years,
      data = tcga_t_model
    )
    
    summary(tcga_mylk_full)
    summary(tcga_mylk_full)$coefficients
    cor(
      tcga_t_model[
        ,
        c("MYLK", "PRKCB", "PDGFC", "VEGFD", "MYLK4")
      ],
      method = "spearman"
    )
    vif_manual <- function(data, vars) {
      
      sapply(vars, function(v) {
        
        others <- setdiff(vars, v)
        
        formula <- as.formula(
          paste(v, "~", paste(others, collapse = " + "))
        )
        
        r2 <- summary(
          lm(formula, data = data)
        )$r.squared
        
        1 / (1 - r2)
      })
    }
    
    candidate_vars <- c(
      "MYLK",
      "PRKCB",
      "PDGFC",
      "VEGFD",
      "MYLK4"
    )
    
    vif_manual(
      tcga_t_model,
      candidate_vars
    )
    cor(
      geo_gleason[
        ,
        c("MYLK", "PRKCB", "PDGFC")
      ],
      method = "spearman"
    )
    cor.test(
      geo_gleason$MYLK,
      geo_gleason$PRKCB,
      method = "spearman",
      exact = FALSE
    )
    
    cor.test(
      geo_gleason$MYLK,
      geo_gleason$PDGFC,
      method = "spearman",
      exact = FALSE
    )
    
    cor.test(
      geo_gleason$PRKCB,
      geo_gleason$PDGFC,
      method = "spearman",
      exact = FALSE
    )
    tcga_network <- tcga_t_model[
      ,
      c("MYLK", "PRKCB", "PDGFC")
    ]
    
    tcga_pca <- prcomp(
      tcga_network,
      scale. = TRUE
    )
    
    summary(tcga_pca)
    tcga_pca_df <- data.frame(
      patient_id = rownames(tcga_t_model),
      Gleason = tcga_t_model$gleason_score,
      PC1 = tcga_pca$x[, 1]
    )
    
    cor.test(
      tcga_pca_df$PC1,
      tcga_pca_df$Gleason,
      method = "spearman",
      exact = FALSE
    )
    tcga_pca$rotation
    geo_network <- geo_gleason[
      ,
      c("MYLK", "PRKCB", "PDGFC")
    ]
    
    geo_pca <- prcomp(
      geo_network,
      scale. = TRUE
    )
    
    summary(geo_pca)    
    geo_pca_df <- data.frame(
      patient_id = rownames(geo_gleason),
      Gleason = geo_gleason$gleason_score,
      PC1 = geo_pca$x[, 1]
    )
    
    cor.test(
      geo_pca_df$PC1,
      geo_pca_df$Gleason,
      method = "spearman",
      exact = FALSE
    )
    geo_pca$rotation
    colnames(tcga_t_model)    
    head(rownames(tcga_t_model))   
    head(tcga_pca_df)    
    tcga_pca_df <- data.frame(
      patient_id = tcga_t_model$patient_id,
      Gleason = tcga_t_model$gleason_score,
      age_years = tcga_t_model$age_years,
      PC1 = tcga_pca$x[, 1]
    )
    
    head(tcga_pca_df)   
    tcga_pc1_model <- lm(
      Gleason ~ PC1 + age_years,
      data = tcga_pca_df
    )
    
    summary(tcga_pc1_model)
    summary(tcga_pc1_model)$coefficients
    tcga_pc1_beta <- coef(tcga_pc1_model)["PC1"]
    
    tcga_pc1_CI <- confint(
      tcga_pc1_model,
      "PC1"
    )
    
    tcga_pc1_beta
    tcga_pc1_CI  
    geo_pca_df <- data.frame(
      patient_id = geo_gleason$geo_accession,
      Gleason = geo_gleason$gleason_score,
      age_years = geo_gleason$age_years,
      PC1 = geo_pca$x[, 1]
    )
    
    geo_pc1_model <- lm(
      Gleason ~ PC1 + age_years,
      data = geo_pca_df
    )
    
    summary(geo_pc1_model)
    summary(geo_pc1_model)$coefficients
    tcga_gene_means <- colMeans(
      tcga_network,
      na.rm = TRUE
    )
    
    tcga_gene_sds <- apply(
      tcga_network,
      2,
      sd,
      na.rm = TRUE
    )
    
    tcga_gene_means
    tcga_gene_sds   
    tcga_weights <- c(
      MYLK = 0.5896508,
      PRKCB = 0.5686800,
      PDGFC = 0.5735111
    )
    
    tcga_signature <- as.matrix(
      scale(
        tcga_network,
        center = tcga_gene_means,
        scale = tcga_gene_sds
      )
    ) %*% tcga_weights
    
    tcga_signature <- as.numeric(tcga_signature)
    
    summary(tcga_signature)
    cor(
      tcga_signature,
      -tcga_pca$x[, 1]
    )
    geo_network <- geo_gleason[
      ,
      c("MYLK", "PRKCB", "PDGFC")
    ]
    
    geo_signature <- as.matrix(
      scale(
        geo_network,
        center = tcga_gene_means,
        scale = tcga_gene_sds
      )
    ) %*% tcga_weights
    
    geo_signature <- as.numeric(geo_signature)
    
    summary(geo_signature)
    geo_signature_model <- lm(
      gleason_score ~ geo_signature + age_years,
      data = geo_gleason
    )
    
    summary(geo_signature_model)
    summary(geo_signature_model)$coefficients
    head(rownames(geo_gleason))
    head(rownames(geo_ordinal))
    
    all(rownames(geo_gleason) == rownames(geo_ordinal))    
    table(geo_ordinal$T_stage_simple)
    geo_ordinal$signature <- geo_signature    
    summary(geo_ordinal$signature)
    
    sum(is.na(geo_ordinal$signature))   
    geo_ordinal$T_advanced <- ifelse(
      geo_ordinal$T_stage_simple == "T2",
      0,
      1
    )
    
    table(geo_ordinal$T_advanced)
    geo_signature_T_model <- glm(
      T_advanced ~ signature + gleason_score + age_years,
      data = geo_ordinal,
      family = binomial
    )
    
    summary(geo_signature_T_model)
    geo_T_coef <- summary(
      geo_signature_T_model
    )$coefficients
    
    geo_T_OR <- exp(
      coef(geo_signature_T_model)
    )
    
    geo_T_CI <- exp(
      confint.default(geo_signature_T_model)
    )
    
    geo_T_results <- data.frame(
      OR = geo_T_OR,
      CI_lower = geo_T_CI[, 1],
      CI_upper = geo_T_CI[, 2],
      p_value = geo_T_coef[, 4]
    )
    
    geo_T_results
    ls()
    class(vsd_tumor)
    dim(vsd_tumor)
    
    class(vsd)
    dim(vsd)
    
    class(counts_filtered)
    dim(counts_filtered)
    
    class(tcga_network)
    dim(tcga_network)
    
    class(gene_expr_geo)
    dim(gene_expr_geo) 
    head(rownames(vsd_tumor))
    head(colnames(vsd_tumor))
    
    head(rownames(vsd))
    head(colnames(vsd))
    head(gene_symbols)
    head(convert_to_symbols)
    grep(
      "MYLK|PRKCB|PDGFC",
      gene_symbols,
      value = TRUE
    )
    signature_ensembl <- c(
      MYLK = "ENSG00000065534",
      PRKCB = "ENSG00000166501",
      PDGFC = "ENSG00000145431"
    )
    
    vsd_gene_ids <- sub(
      "\\..*$",
      "",
      rownames(vsd_tumor)
    )
    
    signature_rows <- match(
      signature_ensembl,
      vsd_gene_ids
    )
    
    signature_rows
    tcga_signature_expr <- t(
      vsd_tumor[
        signature_rows,
        ,
        drop = FALSE
      ]
    )
    
    colnames(tcga_signature_expr) <- names(signature_ensembl)
    
    dim(tcga_signature_expr)
    head(tcga_signature_expr)
    summary(tcga_signature_expr)
    cor(
      tcga_signature_expr,
      method = "spearman"
    ) 
    
    tcga_weights <- c(
      MYLK = 0.5896508,
      PRKCB = 0.5686800,
      PDGFC = 0.5735111
    )
    tcga_signature_means <- colMeans(
      tcga_signature_expr
    )
    
    tcga_signature_sds <- apply(
      tcga_signature_expr,
      2,
      sd
    )
    
    tcga_signature_means
    tcga_signature_sds
    tcga_signature_final <- as.numeric(
      scale(
        tcga_signature_expr,
        center = tcga_signature_means,
        scale = tcga_signature_sds
      ) %*%
        tcga_weights
    )
    summary(tcga_signature_final)
    cor(
      tcga_signature_final,
      -tcga_pca$x[, 1]
    )   
    length(tcga_signature_final)
    length(tcga_pca$x[, 1])
    
    dim(tcga_signature_expr)
    dim(tcga_network)
    summary(tcga_signature_final)
    tcga_final_df <- data.frame(
      patient_id = rownames(tcga_signature_expr),
      signature = tcga_signature_final,
      gleason_score = tcga_t_model$gleason_score,
      age_years = tcga_t_model$age_years
    )   
    tcga_signature_means <- colMeans(
      tcga_signature_expr
    )
    
    tcga_signature_sds <- apply(
      tcga_signature_expr,
      2,
      sd
    )
    
    tcga_signature_final <- as.numeric(
      scale(
        tcga_signature_expr,
        center = tcga_signature_means,
        scale = tcga_signature_sds
      ) %*%
        tcga_weights
    )
    
    length(tcga_signature_final) 
    summary(tcga_signature_final)
    tcga_t_model$patient_id    
    tcga_signature_df <- data.frame(
      patient_id = sub(
        "-01A.*$",
        "",
        rownames(tcga_signature_expr)
      ),
      signature = tcga_signature_final,
      stringsAsFactors = FALSE
    )
    
    dim(tcga_signature_df)
    head(tcga_signature_df)    
    tcga_final_df <- merge(
      tcga_t_model[
        ,
        c(
          "patient_id",
          "gleason_score",
          "age_years"
        )
      ],
      tcga_signature_df,
      by = "patient_id",
      all = FALSE
    )
    
    dim(tcga_final_df)
    colSums(
      is.na(
        tcga_final_df[
          ,
          c(
            "signature",
            "gleason_score",
            "age_years"
          )
        ]
      )
    )
    tcga_final_gleason_model <- lm(
      gleason_score ~ signature + age_years,
      data = tcga_final_df
    )
    
    summary(tcga_final_gleason_model)
    setdiff(
      tcga_t_model$patient_id,
      tcga_signature_df$patient_id
    )
    length(
      setdiff(
        tcga_t_model$patient_id,
        tcga_signature_df$patient_id
      )
    )
    length(
      setdiff(
        tcga_signature_df$patient_id,
        tcga_t_model$patient_id
      )
    )
    missing_patients <- setdiff(
      tcga_t_model$patient_id,
      tcga_signature_df$patient_id
    )
    
    missing_patients
    tcga_signature_df[
      tcga_signature_df$patient_id %in% missing_patients,
    ]
    colnames(vsd_tumor)[
      sapply(
        colnames(vsd_tumor),
        function(x) {
          any(
            sapply(
              missing_patients,
              function(p) startsWith(x, p)
            )
          )
        }
      )
    ]
    tcga_signature_df <- data.frame(
      patient_id = sub(
        "^((TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})).*$",
        "\\1",
        rownames(tcga_signature_expr)
      ),
      signature = tcga_signature_final,
      stringsAsFactors = FALSE
    )
    dim(tcga_signature_df)
    
    length(
      intersect(
        tcga_t_model$patient_id,
        tcga_signature_df$patient_id
      )
    )
    
    setdiff(
      tcga_t_model$patient_id,
      tcga_signature_df$patient_id
    )
    tcga_final_df <- merge(
      tcga_t_model[
        ,
        c(
          "patient_id",
          "gleason_score",
          "age_years"
        )
      ],
      tcga_signature_df,
      by = "patient_id",
      all = FALSE
    )
    
    dim(tcga_final_df)
    
    colSums(
      is.na(
        tcga_final_df[
          ,
          c(
            "signature",
            "gleason_score",
            "age_years"
          )
        ]
      )
    )
    tcga_final_gleason_model <- lm(
      gleason_score ~ signature + age_years,
      data = tcga_final_df
    )
    
    summary(tcga_final_gleason_model)
    tcga_final_beta <- coef(
      tcga_final_gleason_model
    )["signature"]
    
    tcga_final_CI <- confint(
      tcga_final_gleason_model,
      "signature"
    )
    
    tcga_final_beta
    tcga_final_CI
    cor.test(
      tcga_final_df$signature,
      tcga_final_df$gleason_score,
      method = "spearman",
      exact = FALSE
    )
    table(
      tcga_t_model$T_advanced,
      useNA = "ifany"
    )
    
    table(
      tcga_t_model$ajcc_pathologic_t,
      useNA = "ifany"
    )
    length(tcga_t_model$patient_id)
    
    length(tcga_signature_df$patient_id)
    
    length(
      intersect(
        tcga_t_model$patient_id,
        tcga_signature_df$patient_id
      )
    )
    tcga_T_validation <- merge(
      tcga_t_model[
        ,
        c(
          "patient_id",
          "T_advanced",
          "gleason_score",
          "age_years"
        )
      ],
      tcga_signature_df,
      by = "patient_id",
      all = FALSE
    )
    
    dim(tcga_T_validation)
    
    colSums(
      is.na(
        tcga_T_validation[
          ,
          c(
            "signature",
            "T_advanced",
            "gleason_score",
            "age_years"
          )
        ]
      )
    )
    tcga_signature_T_model <- glm(
      T_advanced ~ signature + gleason_score + age_years,
      data = tcga_T_validation,
      family = binomial
    )
    
    summary(tcga_signature_T_model)
    tcga_T_coef <- summary(
      tcga_signature_T_model
    )$coefficients
    
    tcga_T_OR <- exp(
      coef(tcga_signature_T_model)
    )
    
    tcga_T_CI <- exp(
      confint.default(tcga_signature_T_model)
    )
    
    tcga_T_results <- data.frame(
      OR = tcga_T_OR,
      CI_lower = tcga_T_CI[, 1],
      CI_upper = tcga_T_CI[, 2],
      p_value = tcga_T_coef[, 4]
    )
    
    tcga_T_results
    wilcox.test(
      signature ~ T_advanced,
      data = tcga_T_validation
    )
    aggregate(
      signature ~ T_advanced,
      data = tcga_T_validation,
      FUN = median
    )
    table(
      tcga_final_df$gleason_score,
      useNA = "ifany"
     )
    tcga_final_df$gleason_ordinal <- factor(
      tcga_final_df$gleason_score,
      levels = sort(unique(tcga_final_df$gleason_score)),
      ordered = TRUE
    )
    
    table(tcga_final_df$gleason_ordinal)
    library(MASS)
    
    tcga_signature_ordinal <- polr(
      gleason_ordinal ~ signature + age_years,
      data = tcga_final_df,
      Hess = TRUE
    )
    
    summary(tcga_signature_ordinal)
    tcga_ord_coef <- coef(
      summary(tcga_signature_ordinal)
    )
    
    tcga_ord_p <- 2 * pnorm(
      abs(tcga_ord_coef["signature", "t value"]),
      lower.tail = FALSE
    )
    
    tcga_ord_OR <- exp(
      coef(tcga_signature_ordinal)["signature"]
    )
    
    tcga_ord_CI <- exp(
      confint(
        tcga_signature_ordinal,
        parm = "signature"
      )
    )
    
    tcga_ord_OR
    tcga_ord_CI
    tcga_ord_p
    tcga_final_df_no10 <- subset(
      tcga_final_df,
      gleason_score < 10
    )
    
    tcga_signature_ordinal_no10 <- polr(
      factor(
        gleason_score,
        levels = sort(unique(gleason_score)),
        ordered = TRUE
      ) ~ signature + age_years,
      data = tcga_final_df_no10,
      Hess = TRUE
    )
    
    summary(tcga_signature_ordinal_no10)
    coef_no10 <- coef(
      summary(tcga_signature_ordinal_no10)
    )
    
    p_no10 <- 2 * pnorm(
      abs(coef_no10["signature", "t value"]),
      lower.tail = FALSE
    )
    
    OR_no10 <- exp(
      coef(tcga_signature_ordinal_no10)["signature"]
    )
    
    CI_no10 <- exp(
      confint(
        tcga_signature_ordinal_no10,
        parm = "signature"
      )
    )
    
    OR_no10
    CI_no10
    p_no10
    tcga_signature_sd <- sd(
      tcga_final_df$signature,
      na.rm = TRUE
    )
    
    tcga_signature_sd
    tcga_beta_1SD <- coef(
      tcga_final_gleason_model
    )["signature"] * tcga_signature_sd
    
    tcga_CI_1SD <- confint(
      tcga_final_gleason_model,
      "signature"
    ) * tcga_signature_sd
    
    tcga_beta_1SD
    tcga_CI_1SD
    geo_ordinal$signature
    geo_signature_sd <- sd(
      geo_ordinal$signature,
      na.rm = TRUE
    )
    
    geo_signature_sd    
    geo_beta_1SD <- coef(
      geo_signature_model
    )["geo_signature"] * geo_signature_sd
    
    geo_CI_1SD <- confint(
      geo_signature_model,
      "geo_signature"
    ) * geo_signature_sd
    
    geo_beta_1SD
    geo_CI_1SD
    standardized_results <- data.frame(
      Cohort = c("GEO", "TCGA"),
      N = c(
        nrow(geo_ordinal),
        nrow(tcga_final_df)
      ),
      Beta_per_SD = c(
        geo_beta_1SD,
        tcga_beta_1SD
      ),
      CI_lower = c(
        geo_CI_1SD[1],
        tcga_CI_1SD[1]
      ),
      CI_upper = c(
        geo_CI_1SD[2],
        tcga_CI_1SD[2]
      ),
      P_value = c(
        summary(geo_signature_model)$coefficients[
          "geo_signature",
          "Pr(>|t|)"
        ],
        summary(tcga_final_gleason_model)$coefficients[
          "signature",
          "Pr(>|t|)"
        ]
      )
    )
    
    standardized_results
    tcga_mylk_model_final <- lm(
      gleason_score ~ MYLK + age_years,
      data = tcga_t_model
    )
    
    summary(tcga_mylk_model_final)
    summary(tcga_final_gleason_model)
    anova(
      tcga_mylk_model_final,
      tcga_final_gleason_model
    )
    library(ggplot2)
    
    # GEO
    geo_plot <- data.frame(
      Cohort = "GEO",
      Gleason = geo_ordinal$gleason_score,
      signature = geo_ordinal$signature,
      age_years = geo_ordinal$age_years
    )
    
    # TCGA
    tcga_plot <- data.frame(
      Cohort = "TCGA",
      Gleason = tcga_final_df$gleason_score,
      signature = tcga_final_df$signature,
      age_years = tcga_final_df$age_years
    )
    
    # Standardize within each cohort
    geo_plot$signature_z <- as.numeric(
      scale(geo_plot$signature)
    )
    
    tcga_plot$signature_z <- as.numeric(
      scale(tcga_plot$signature)
    )
    
    gleason_plot <- rbind(
      geo_plot,
      tcga_plot
    )
    
    table(gleason_plot$Cohort)
    p_gleason <- ggplot(
      gleason_plot,
      aes(
        x = signature_z,
        y = Gleason
      )
    ) +
      geom_jitter(
        width = 0.05,
        height = 0.08,
        alpha = 0.45,
        size = 1.8
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.9
      ) +
      facet_wrap(
        ~ Cohort,
        nrow = 1
      ) +
      scale_y_continuous(
        breaks = 6:10
      ) +
      labs(
        x = "MYLK–PRKCB–PDGFC signature (Z-score)",
        y = "Gleason score"
      ) +
      theme_classic(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(
          face = "bold"
        )
      )
    
    p_gleason
    annotation_df <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(-2.7, -2.7),
      y = c(9.8, 9.8),
      label = c(
        "β/SD = −0.158\nP = 0.000993",
        "β/SD = −0.167\nP = 0.000277"
      )
    )
    
    p_gleason_final <- ggplot(
      gleason_plot,
      aes(
        x = signature_z,
        y = Gleason
      )
    ) +
      geom_jitter(
        width = 0.05,
        height = 0.08,
        alpha = 0.45,
        size = 1.8
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.9
      ) +
      facet_wrap(
        ~ Cohort,
        nrow = 1
      ) +
      geom_text(
        data = annotation_df,
        aes(
          x = x,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 4
      ) +
      scale_y_continuous(
        breaks = 6:10
      ) +
      labs(
        x = "MYLK–PRKCB–PDGFC signature (Z-score)",
        y = "Gleason score"
      ) +
      theme_classic(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(
          face = "bold"
        )
      )
    
    p_gleason_final
    ggsave(
      "Figure_2_GEO_TCGA_Gleason_signature.png",
      p_gleason_final,
      width = 10,
      height = 5.5,
      dpi = 600
    )
    
    ggsave(
      "Figure_2_GEO_TCGA_Gleason_signature.pdf",
      p_gleason_final,
      width = 10,
      height = 5.5
    )
    annotation_df <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(-2.7, -2.7),
      y = c(9.85, 9.85),
      label = c(
        "N = 111\nβ/SD = −0.158\nP < 0.001",
        "N = 448\nβ/SD = −0.167\nP < 0.001"
      )
    )
    
    p_gleason_final <- ggplot(
      gleason_plot,
      aes(x = signature_z, y = Gleason)
    ) +
      geom_jitter(
        width = 0.05,
        height = 0.06,
        alpha = 0.30,
        size = 1.4
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.9
      ) +
      facet_wrap(
        ~ Cohort,
        nrow = 1
      ) +
      geom_text(
        data = annotation_df,
        aes(
          x = x,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 4
      ) +
      scale_y_continuous(
        breaks = 6:10
      ) +
      labs(
        x = "MYLK–PRKCB–PDGFC signature (Z-score)",
        y = "Gleason score"
      ) +
      theme_classic(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(
          face = "bold"
        )
      )
    
    p_gleason_final
    annotation_df <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(-2.7, 0.4),
      y = c(9.85, 9.85),
      label = c(
        "N = 111\nβ/SD = −0.158\nP < 0.001",
        "N = 448\nβ/SD = −0.167\nP < 0.001"
      )
    )
    # ================================
    # FIGURE 3A — SIGNATURE vs GLEASON
    # ================================
    
    # GEO
    geo_gleason_plot <- geo_ordinal[
      ,
      c("signature", "gleason_score")
    ]
    
    geo_gleason_plot$Cohort <- "GEO"
    
    # TCGA
    tcga_gleason_plot <- tcga_final_df[
      ,
      c("signature", "gleason_score")
    ]
    
    tcga_gleason_plot$Cohort <- "TCGA"
    
    # Combine
    fig3_gleason <- rbind(
      geo_gleason_plot,
      tcga_gleason_plot
    )
    
    fig3_gleason$Gleason <- factor(
      fig3_gleason$gleason_score,
      levels = c(6, 7, 8, 9, 10)
    )
    
    table(
      fig3_gleason$Cohort,
      fig3_gleason$Gleason,
      useNA = "ifany"
    )
    # GEO
    geo_gleason_kw <- kruskal.test(
      signature ~ Gleason,
      data = subset(fig3_gleason, Cohort == "GEO")
    )
    
    # TCGA
    tcga_gleason_kw <- kruskal.test(
      signature ~ Gleason,
      data = subset(fig3_gleason, Cohort == "TCGA")
    )
    
    geo_gleason_kw
    tcga_gleason_kw
    aggregate(
      signature ~ Cohort + Gleason,
      data = fig3_gleason,
      FUN = median
    )
    table(tcga_t_model$ajcc_pathologic_t)
    T2 = 0
    T3/T4 = 1    
    # ================================
    # FIGURE 3B — SIGNATURE vs T STAGE
    # ================================
    
    # GEO
    geo_T_plot <- data.frame(
      signature = geo_ordinal$signature,
      T_stage = ifelse(
        geo_ordinal$T_stage_simple == "T2",
        "T2",
        "T3/T4"
      ),
      Cohort = "GEO"
    )
    
    # TCGA
    tcga_T_plot <- data.frame(
      signature = tcga_T_validation$signature,
      T_stage = ifelse(
        tcga_T_validation$T_advanced == 0,
        "T2",
        "T3/T4"
      ),
      Cohort = "TCGA"
    )
    
    # Combine
    fig3_T <- rbind(
      geo_T_plot,
      tcga_T_plot
    )
    
    # Set order
    fig3_T$T_stage <- factor(
      fig3_T$T_stage,
      levels = c("T2", "T3/T4")
    )
    
    # Check
    table(
      fig3_T$Cohort,
      fig3_T$T_stage,
      useNA = "ifany"
    )
    
    geo_T_test <- wilcox.test(
      signature ~ T_stage,
      data = subset(fig3_T, Cohort == "GEO")
    )
    
    tcga_T_test <- wilcox.test(
      signature ~ T_stage,
      data = subset(fig3_T, Cohort == "TCGA")
    )
    
    geo_T_test
    tcga_T_test
    aggregate(
      signature ~ Cohort + T_stage,
      data = fig3_T,
      FUN = median
    )
    # Install if needed
    # install.packages("FSA")
    
    library(FSA)
    
    # GEO
    geo_gleason_dunn <- dunnTest(
      signature ~ Gleason,
      data = subset(fig3_gleason, Cohort == "GEO"),
      method = "bh"
    )
    
    geo_gleason_dunn
    # ================================
    # Dunn-style pairwise Wilcoxon tests
    # with BH correction
    # ================================
    
    # GEO
    geo_pairwise <- pairwise.wilcox.test(
      x = subset(fig3_gleason, Cohort == "GEO")$signature,
      g = subset(fig3_gleason, Cohort == "GEO")$Gleason,
      p.adjust.method = "BH",
      exact = FALSE
    )
    
    geo_pairwise
    # TCGA
    tcga_pairwise <- pairwise.wilcox.test(
      x = subset(fig3_gleason, Cohort == "TCGA")$signature,
      g = subset(fig3_gleason, Cohort == "TCGA")$Gleason,
      p.adjust.method = "BH",
      exact = FALSE
    )
    
    tcga_pairwise
    range(
      fig3_gleason$signature[
        fig3_gleason$Cohort == "GEO"
      ]
    )
    
    range(
      fig3_gleason$signature[
        fig3_gleason$Cohort == "TCGA"
      ]
    )
    library(ggplot2)
    
    fig3_gleason$Gleason <- factor(
      fig3_gleason$gleason_score,
      levels = c(6, 7, 8, 9, 10)
    )
    
    kw_annotation <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(6.15, 6.15),
      y = c(-2.65, 3.0),
      label = c(
        "Kruskal–Wallis P = 0.0176",
        "Kruskal–Wallis P < 0.001"
      )
    )
    
    p_fig3 <- ggplot(
      fig3_gleason,
      aes(
        x = Gleason,
        y = signature
      )
    ) +
      geom_boxplot(
        width = 0.55,
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.35,
        size = 1.3
      ) +
      facet_wrap(
        ~ Cohort,
        nrow = 1,
        scales = "free_y"
      ) +
      geom_text(
        data = kw_annotation,
        aes(
          x = x,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = 4
      ) +
      labs(
        x = "Gleason score",
        y = "MYLK–PRKCB–PDGFC signature"
      ) +
      theme_classic(
        base_size = 14
      ) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(
          face = "bold"
        )
      )
    
    p_fig3
    # ============================================
    # FIGURE 3 — FINAL ANNOTATION
    # ============================================
    
    library(ggplot2)
    
    # Get panel-specific y positions
    annotation_df <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(6.15, 6.15),
      y = c(
        max(fig3_gleason$signature[
          fig3_gleason$Cohort == "GEO"
        ]) - 0.05,
        
        max(fig3_gleason$signature[
          fig3_gleason$Cohort == "TCGA"
        ]) - 0.05
      ),
      label = c(
        "P = 0.0176",
        "P < 0.001"
      )
    )
    
    p_fig3 <- ggplot(
      fig3_gleason,
      aes(
        x = Gleason,
        y = signature
      )
    ) +
      geom_boxplot(
        width = 0.55,
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.35,
        size = 1.3
      ) +
      facet_wrap(
        ~ Cohort,
        nrow = 1,
        scales = "free_y"
      ) +
      geom_text(
        data = annotation_df,
        aes(
          x = x,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 4,
        fontface = "bold"
      ) +
      labs(
        x = "Gleason score",
        y = "MYLK–PRKCB–PDGFC signature"
      ) +
      theme_classic(
        base_size = 14
      ) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(
          face = "bold"
        )
      )
    
    p_fig3
    
   
    head(gsea_kegg_df, 15)
    head(gsea_df, 15)    
    head(ego_up_df, 15)
    head(ego_down_df, 15)   
    # ================================
    # FIGURE 4 — GSEA PATHWAY SUMMARY
    # ================================
    
    # GO
    gsea_go_plot <- gsea_df |>
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust)
      ) |>
      dplyr::arrange(p.adjust)
    
    # KEGG
    gsea_kegg_plot <- gsea_kegg_df |>
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust)
      ) |>
      dplyr::arrange(p.adjust)
    
    dim(gsea_go_plot)
    dim(gsea_kegg_plot)
    
    head(
      gsea_go_plot[, c(
        "Description",
        "NES",
        "p.adjust"
      )],
      20
    )
    
    head(
      gsea_kegg_plot[, c(
        "Description",
        "NES",
        "p.adjust"
      )],
      20
    )
    # ============================================================
    # FIGURE 4 — GSEA SUMMARY
    # ============================================================
    
    library(dplyr)
    library(ggplot2)
    library(stringr)
    
    # ------------------------------------------------------------
    # 1. Select top GO pathways
    # ------------------------------------------------------------
    
    go_up <- gsea_go_plot %>%
      filter(NES > 0, p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = 6)
    
    go_down <- gsea_go_plot %>%
      filter(NES < 0, p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = 6)
    
    go_fig <- bind_rows(go_up, go_down) %>%
      mutate(
        Direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
        Description = str_wrap(Description, width = 40),
        Description = factor(
          Description,
          levels = Description[order(NES)]
        ),
        neg_log10_FDR = -log10(p.adjust)
      )
    
    
    # ------------------------------------------------------------
    # 2. Select top KEGG pathways
    # ------------------------------------------------------------
    
    kegg_up <- gsea_kegg_plot %>%
      filter(NES > 0, p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = 6)
    
    kegg_down <- gsea_kegg_plot %>%
      filter(NES < 0, p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = 6)
    
    kegg_fig <- bind_rows(kegg_up, kegg_down) %>%
      mutate(
        Direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
        Description = str_wrap(Description, width = 40),
        Description = factor(
          Description,
          levels = Description[order(NES)]
        ),
        neg_log10_FDR = -log10(p.adjust)
      )
    
    
    # ------------------------------------------------------------
    # 3. GO plot
    # ------------------------------------------------------------
    
    p_GO <- ggplot(
      go_fig,
      aes(
        x = NES,
        y = Description,
        size = neg_log10_FDR,
        color = NES
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.5
      ) +
      geom_point(alpha = 0.9) +
      scale_size_continuous(
        name = expression(-log[10]("FDR"))
      ) +
      scale_color_gradient2(
        low = "steelblue4",
        mid = "grey70",
        high = "firebrick3",
        midpoint = 0,
        name = "NES"
      ) +
      labs(
        title = "GO Biological Process GSEA",
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 14
        ),
        axis.text.y = element_text(
          size = 10
        ),
        axis.text.x = element_text(
          size = 10
        ),
        legend.title = element_text(
          face = "bold"
        )
      )
    
    
    # ------------------------------------------------------------
    # 4. KEGG plot
    # ------------------------------------------------------------
    
    p_KEGG <- ggplot(
      kegg_fig,
      aes(
        x = NES,
        y = Description,
        size = neg_log10_FDR,
        color = NES
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.5
      ) +
      geom_point(alpha = 0.9) +
      scale_size_continuous(
        name = expression(-log[10]("FDR"))
      ) +
      scale_color_gradient2(
        low = "steelblue4",
        mid = "grey70",
        high = "firebrick3",
        midpoint = 0,
        name = "NES"
      ) +
      labs(
        title = "KEGG GSEA",
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 14
        ),
        axis.text.y = element_text(
          size = 10
        ),
        axis.text.x = element_text(
          size = 10
        ),
        legend.title = element_text(
          face = "bold"
        )
      )
    
    
    # Display
    p_GO
    p_KEGG
    
    library(ggplot2)
    
    go_test <- gsea_go_plot[
      order(gsea_go_plot$p.adjust),
      c("Description", "NES", "p.adjust")
    ]
    
    go_test <- go_test[
      is.finite(go_test$NES) &
        is.finite(go_test$p.adjust) &
        go_test$p.adjust > 0,
    ]
    
    go_test <- head(go_test, 15)
    
    ggplot(
      go_test,
      aes(
        x = NES,
        y = reorder(Description, NES)
      )
    ) +
      geom_point(size = 4) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      theme_classic() +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        title = "GO GSEA"
      )
    library(ggplot2)
    
    # Remove invalid values
    kegg_test <- gsea_kegg_plot[
      is.finite(gsea_kegg_plot$NES) &
        is.finite(gsea_kegg_plot$p.adjust) &
        gsea_kegg_plot$p.adjust > 0,
    ]
    
    # Select the top 8 positive and top 8 negative pathways
    kegg_up <- kegg_test[
      kegg_test$NES > 0,
    ]
    kegg_up <- kegg_up[
      order(kegg_up$p.adjust),
    ]
    kegg_up <- head(kegg_up, 8)
    
    kegg_down <- kegg_test[
      kegg_test$NES < 0,
    ]
    kegg_down <- kegg_down[
      order(kegg_down$p.adjust),
    ]
    kegg_down <- head(kegg_down, 8)
    
    kegg_final <- rbind(
      kegg_up,
      kegg_down
    )
    
    # Order pathways by NES
    kegg_final$Description <- factor(
      kegg_final$Description,
      levels = kegg_final$Description[
        order(kegg_final$NES)
      ]
    )
    
    # Plot
    ggplot(
      kegg_final,
      aes(
        x = NES,
        y = Description
      )
    ) +
      geom_point(size = 4) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      theme_classic() +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        title = "KEGG GSEA"
      )
    library(ggplot2)
    
    # ============================================================
    # FINAL GSEA FIGURE — GO + KEGG
    # ============================================================
    
    # ---------- GO ----------
    go_plot <- gsea_go_plot[
      is.finite(gsea_go_plot$NES) &
        is.finite(gsea_go_plot$p.adjust) &
        gsea_go_plot$p.adjust > 0,
    ]
    
    # Top 8 positive
    go_up <- go_plot[go_plot$NES > 0, ]
    go_up <- go_up[order(go_up$p.adjust), ]
    go_up <- head(go_up, 8)
    
    # Top 8 negative
    go_down <- go_plot[go_plot$NES < 0, ]
    go_down <- go_down[order(go_down$p.adjust), ]
    go_down <- head(go_down, 8)
    
    go_final <- rbind(go_up, go_down)
    
    go_final$logFDR <- -log10(go_final$p.adjust)
    
    go_final$Description <- factor(
      go_final$Description,
      levels = go_final$Description[
        order(go_final$NES)
      ]
    )
    
    
    # ---------- KEGG ----------
    kegg_plot <- gsea_kegg_plot[
      is.finite(gsea_kegg_plot$NES) &
        is.finite(gsea_kegg_plot$p.adjust) &
        gsea_kegg_plot$p.adjust > 0,
    ]
    
    # Top 8 positive
    kegg_up <- kegg_plot[kegg_plot$NES > 0, ]
    kegg_up <- kegg_up[order(kegg_up$p.adjust), ]
    kegg_up <- head(kegg_up, 8)
    
    # Top 8 negative
    kegg_down <- kegg_plot[kegg_plot$NES < 0, ]
    kegg_down <- kegg_down[order(kegg_down$p.adjust), ]
    kegg_down <- head(kegg_down, 8)
    
    kegg_final <- rbind(kegg_up, kegg_down)
    
    kegg_final$logFDR <- -log10(kegg_final$p.adjust)
    
    kegg_final$Description <- factor(
      kegg_final$Description,
      levels = kegg_final$Description[
        order(kegg_final$NES)
      ]
    )
    
    
    # ============================================================
    # GO PANEL
    # ============================================================
    
    p_GO <- ggplot(
      go_final,
      aes(
        x = NES,
        y = Description,
        size = logFDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      theme_classic() +
      labs(
        title = "A  GO GSEA",
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        size = expression(-log[10]~"(FDR)")
      ) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 16
        ),
        axis.title.x = element_text(
          size = 12
        ),
        axis.text.y = element_text(
          size = 10
        ),
        legend.title = element_text(
          size = 10
        )
      )
    
    
    # ============================================================
    # KEGG PANEL
    # ============================================================
    
    p_KEGG <- ggplot(
      kegg_final,
      aes(
        x = NES,
        y = Description,
        size = logFDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      theme_classic() +
      labs(
        title = "B  KEGG GSEA",
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        size = expression(-log[10]~"(FDR)")
      ) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 16
        ),
        axis.title.x = element_text(
          size = 12
        ),
        axis.text.y = element_text(
          size = 10
        ),
        legend.title = element_text(
          size = 10
        )
      )
    
    
    # Display
    p_GO
    p_KEGG
    library(patchwork)
    
    final_GSEA_figure <- p_GO + p_KEGG
    
    final_GSEA_figure
    ggsave(
      "Figure_GSEA_GO_KEGG.pdf",
      final_GSEA_figure,
      width = 12,
      height = 8
    )
    
    ggsave(
      "Figure_GSEA_GO_KEGG.png",
      final_GSEA_figure,
      width = 12,
      height = 8,
      dpi = 600
    )
    colnames(survival_clinical)
    
    dim(survival_clinical)
    
    head(survival_clinical)
    colnames(tcga_t_model)
    
    head(tcga_t_model[, c(
      "patient_id",
      "gleason_score",
      "age_years"
    )])
    summary(survival_clinical)
    
    table(
      survival_clinical$vital_status,
      useNA = "ifany"
    )
    grep(
      "days|death|survival|status|follow",
      colnames(survival_clinical),
      ignore.case = TRUE,
      value = TRUE
    )
    dim(survival_df)
    
    colnames(survival_df)
    
    head(survival_df)
    
    summary(survival_df)
    class(surv_object)
    
    surv_object
    library(survival)
    
    tcga_survival_signature <- merge(
      survival_clinical[
        ,
        c(
          "submitter_id",
          "survival_time",
          "event",
          "age_at_diagnosis",
          "gleason_score"
        )
      ],
      tcga_signature_df,
      by = "submitter_id",
      all = FALSE
    )
    
    dim(tcga_survival_signature)
    
    colSums(
      is.na(
        tcga_survival_signature[
          ,
          c(
            "survival_time",
            "event",
            "signature",
            "age_at_diagnosis",
            "gleason_score"
          )
        ]
      )
    )
    
    table(tcga_survival_signature$event)
    tcga_survival_signature <- merge(
      survival_clinical[
        ,
        c(
          "submitter_id",
          "survival_time",
          "event",
          "age_at_diagnosis",
          "gleason_score"
        )
      ],
      tcga_signature_df,
      by.x = "submitter_id",
      by.y = "patient_id",
      all = FALSE
    )
    
    dim(tcga_survival_signature)
    head(tcga_survival_signature)
    
    colSums(
      is.na(
        tcga_survival_signature[
          ,
          c(
            "survival_time",
            "event",
            "signature",
            "age_at_diagnosis",
            "gleason_score"
          )
        ]
      )
    )
    
    table(tcga_survival_signature$event)
    library(survival)
    
    tcga_signature_cox <- coxph(
      Surv(survival_time, event) ~ signature,
      data = tcga_survival_signature
    )
    
    summary(tcga_signature_cox)
    tcga_cox_coef <- summary(
      tcga_signature_cox
    )$coefficients
    
    tcga_cox_HR <- exp(
      coef(tcga_signature_cox)
    )
    
    tcga_cox_CI <- exp(
      confint(tcga_signature_cox)
    )
    
    tcga_cox_results <- data.frame(
      HR = tcga_cox_HR,
      CI_lower = tcga_cox_CI[, 1],
      CI_upper = tcga_cox_CI[, 2],
      P_value = tcga_cox_coef[, "Pr(>|z|)"]
    )
    
    tcga_cox_results
    tcga_signature_sd <- sd(
      tcga_survival_signature$signature,
      na.rm = TRUE
    )
    
    tcga_signature_sd
    tcga_survival_signature$signature_z <- 
      tcga_survival_signature$signature / tcga_signature_sd
    tcga_signature_cox_z <- coxph(
      Surv(survival_time, event) ~ signature_z,
      data = tcga_survival_signature
    )
    
    summary(tcga_signature_cox_z)
    sum(
      complete.cases(
        tcga_survival_signature[
          ,
          c(
            "survival_time",
            "event",
            "signature",
            "age_at_diagnosis",
            "gleason_score"
          )
        ]
      )
    )
    tcga_survival_complete <- tcga_survival_signature[
      complete.cases(
        tcga_survival_signature[
          ,
          c(
            "survival_time",
            "event",
            "signature",
            "age_at_diagnosis",
            "gleason_score"
          )
        ]
      ),
    ]
    dim(tcga_survival_complete)
    
    table(tcga_survival_complete$event)
    tcga_signature_cox_adj <- coxph(
      Surv(survival_time, event) ~
        signature_z +
        gleason_score +
        age_at_diagnosis,
      data = tcga_survival_complete
    )
    
    summary(tcga_signature_cox_adj)
    adj_coef <- summary(
      tcga_signature_cox_adj
    )$coefficients
    
    adj_CI <- exp(
      confint(tcga_signature_cox_adj)
    )
    
    tcga_survival_adj_results <- data.frame(
      Variable = rownames(adj_coef),
      HR = exp(adj_coef[, "coef"]),
      CI_lower = adj_CI[, 1],
      CI_upper = adj_CI[, 2],
      P_value = adj_coef[, "Pr(>|z|)"],
      row.names = NULL
    )
    
    tcga_survival_adj_results
    tcga_cox_ph_test <- cox.zph(
      tcga_signature_cox_adj
    )
    
    tcga_cox_ph_test
    plot(tcga_cox_ph_test)
    cor(
      tcga_signature,
      tcga_signature_expr[, "MYLK"],
      method = "spearman"
    )
    
    cor(
      tcga_signature,
      tcga_signature_expr[, "PRKCB"],
      method = "spearman"
    )
    
    cor(
      tcga_signature,
      tcga_signature_expr[, "PDGFC"],
      method = "spearman"
    )    
    length(tcga_signature)
    
    dim(tcga_signature_expr)
    
    head(rownames(tcga_signature_expr))
    
    head(colnames(vsd_tumor))
    head(names(tcga_signature))
    tcga_signature_check <- as.numeric(
      scale(
        tcga_signature_expr,
        center = tcga_gene_means,
        scale = tcga_gene_sds
      ) %*% tcga_weights
    )    
    length(tcga_signature_check)
    
    cor(
      tcga_signature_check,
      tcga_signature_expr[, "MYLK"],
      method = "spearman"
    )
    
    cor(
      tcga_signature_check,
      tcga_signature_expr[, "PRKCB"],
      method = "spearman"
    )
    
    cor(
      tcga_signature_check,
      tcga_signature_expr[, "PDGFC"],
      method = "spearman"
    )
    cor(
      tcga_signature_check,
      tcga_signature,
      method = "pearson"
    )
    tcga_signature_497 <- tcga_signature_check
    
    names(tcga_signature_497) <- rownames(
      tcga_signature_expr
    )
    length(tcga_signature_497)
    
    head(names(tcga_signature_497))
    
    summary(tcga_signature_497)
    tcga_signature_497_df <- data.frame(
      sample_id = names(tcga_signature_497),
      patient_id = sub(
        "^((TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})).*$",
        "\\1",
        names(tcga_signature_497)
      ),
      signature = as.numeric(tcga_signature_497),
      stringsAsFactors = FALSE
    )
    dim(tcga_signature_497_df)
    
    length(unique(tcga_signature_497_df$patient_id))
    
    head(tcga_signature_497_df)
    scale(
      tcga_signature_expr,
      center = tcga_gene_means,
      scale = tcga_gene_sds
    ) %*% tcga_weights
    nrow(tcga_pca$x)
    
    length(tcga_signature_497)
    # Sample IDs used by the PCA
    pca_samples <- colnames(vsd_tumor)[
      seq_len(ncol(vsd_tumor))
    ]
    
    length(pca_samples)
    head(pca_samples)
    dim(vsd_tumor)
    
    dim(tcga_pca$x)
    
    length(colnames(vsd_tumor))
    
    head(colnames(vsd_tumor))
    tcga_signature_pca
    cor(
      tcga_signature_pca,
      tcga_pca$x[, 1],
      method = "pearson"
    )    
    rownames(tcga_pca$x)[1:10]
    length(rownames(tcga_pca$x))   
    head(tcga_signature_497)   
    head(names(tcga_signature_497))    
    sum(
      rownames(tcga_pca$x) %in% names(tcga_signature_497)
    )
    
    sum(
      rownames(tcga_pca$x) %in% rownames(tcga_signature_expr)
    )    
    tcga_signature_pca <- tcga_signature_497[
      rownames(tcga_pca$x)
    ]
    cor(
      tcga_signature_pca,
      tcga_pca$x[, 1],
      method = "pearson"
    )
    # Patient IDs corresponding to the PCA observations
    pca_patient_ids <- tcga_t_model$patient_id[
      as.integer(rownames(tcga_pca$x))
    ]
    
    head(pca_patient_ids)
    length(pca_patient_ids)
    signature_patient_ids <- sub(
      "^((TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})).*$",
      "\\1",
      names(tcga_signature_497)
    )
    
    head(signature_patient_ids)
    tcga_signature_patient <- tcga_signature_497
    names(tcga_signature_patient) <- signature_patient_ids
    tcga_signature_pca <- tcga_signature_patient[
      pca_patient_ids
    ]
    length(tcga_signature_pca)
    
    sum(is.na(tcga_signature_pca))
    
    head(tcga_signature_pca)
    
    head(pca_patient_ids)
    
    head(names(tcga_signature_pca))
    all(names(tcga_signature_pca) == pca_patient_ids)
    sum(is.na(pca_patient_ids))
    sum(is.na(tcga_signature_pca))
    
    which(is.na(pca_patient_ids))
    which(is.na(tcga_signature_pca))    
    head(pca_patient_ids, 20)
    head(names(tcga_signature_pca), 20)
    sum(pca_patient_ids %in% names(tcga_signature_patient))
    sum(!is.na(tcga_signature_pca))    
    missing_pca_patients <- pca_patient_ids[
      is.na(tcga_signature_pca)
    ]
    
    missing_pca_patients    
    length(missing_pca_patients)
    missing_pca_samples <- colnames(vsd_tumor)[
      sapply(
        colnames(vsd_tumor),
        function(x) {
          any(
            startsWith(x, missing_pca_patients)
          )
        }
      )
    ]
    
    missing_pca_samples
    length(missing_pca_samples)    
    tcga_signature_clean <- as.numeric(
      scale(
        tcga_signature_expr,
        center = tcga_gene_means,
        scale = tcga_gene_sds
      ) %*% tcga_weights
    )
    
    names(tcga_signature_clean) <- rownames(tcga_signature_expr)
    length(tcga_signature_clean)
    
    sum(is.na(tcga_signature_clean))
    
    head(names(tcga_signature_clean))
    tcga_signature_clean_patient <- tcga_signature_clean
    
    names(tcga_signature_clean_patient) <- sub(
      "^((TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})).*$",
      "\\1",
      names(tcga_signature_clean_patient)
    )
    length(tcga_signature_clean_patient)
    
    sum(duplicated(names(tcga_signature_clean_patient)))
    
    sum(
      pca_patient_ids %in%
        names(tcga_signature_clean_patient)
    )
    missing_pca_patients <- pca_patient_ids[
      !pca_patient_ids %in% names(tcga_signature_clean_patient)
    ]
    
    missing_pca_patients
    length(missing_pca_patients)
    rownames(tcga_pca$x)[
      !pca_patient_ids %in% names(tcga_signature_clean_patient)
    ]    
    tcga_t_model[
      as.integer(
        rownames(tcga_pca$x)[
          !pca_patient_ids %in% names(tcga_signature_clean_patient)
        ]
      ),
      c("patient_id", "sample_id")
    ]
    tail(tcga_t_model[, c("patient_id", "sample_id")], 60)
    tail(colnames(vsd_tumor), 60)    
    # Match PCA row names to the actual row names of tcga_t_model
    pca_idx <- match(
      rownames(tcga_pca$x),
      rownames(tcga_t_model)
    )
    
    # Check
    length(pca_idx)
    sum(is.na(pca_idx))    
    pca_patient_ids <- tcga_t_model$patient_id[pca_idx]
    head(pca_patient_ids)
    tail(pca_patient_ids)
    length(pca_patient_ids)
    sum(is.na(pca_patient_ids))    
    tcga_signature_pca <- tcga_signature_clean_patient[
      pca_patient_ids
    ]
    length(tcga_signature_pca)
    
    sum(is.na(tcga_signature_pca))
    all(
      names(tcga_signature_pca) == pca_patient_ids
    )
    cor(
      tcga_signature_pca,
      tcga_pca$x[, 1],
      method = "pearson"
    )   
    cor(
      tcga_signature_pca,
      -tcga_pca$x[, 1],
      method = "pearson"
    )
    cor(
      tcga_signature_pca,
      tcga_pca$x[, 1],
      method = "spearman"
    )
    data.frame(
      Pearson_PC1 = cor(
        tcga_signature_pca,
        tcga_pca$x[, 1],
        method = "pearson"
      ),
      Spearman_PC1 = cor(
        tcga_signature_pca,
        tcga_pca$x[, 1],
        method = "spearman"
      ),
      Pearson_negPC1 = cor(
        tcga_signature_pca,
        -tcga_pca$x[, 1],
        method = "pearson"
      )
    )
    # Individual gene models
    tcga_prkcb_model_final <- lm(
      gleason_score ~ PRKCB + age_years,
      data = tcga_t_model
    )
    
    tcga_pdgfc_model_final <- lm(
      gleason_score ~ PDGFC + age_years,
      data = tcga_t_model
    )
    
    # Extract model statistics
    model_comparison <- data.frame(
      Model = c(
        "MYLK",
        "PRKCB",
        "PDGFC",
        "3-gene signature"
      ),
      Beta = c(
        coef(tcga_mylk_model_final)["MYLK"],
        coef(tcga_prkcb_model_final)["PRKCB"],
        coef(tcga_pdgfc_model_final)["PDGFC"],
        coef(tcga_final_gleason_model)["signature"]
      ),
      P_value = c(
        summary(tcga_mylk_model_final)$coefficients["MYLK", "Pr(>|t|)"],
        summary(tcga_prkcb_model_final)$coefficients["PRKCB", "Pr(>|t|)"],
        summary(tcga_pdgfc_model_final)$coefficients["PDGFC", "Pr(>|t|)"],
        summary(tcga_final_gleason_model)$coefficients["signature", "Pr(>|t|)"]
      ),
      Adjusted_R2 = c(
        summary(tcga_mylk_model_final)$adj.r.squared,
        summary(tcga_prkcb_model_final)$adj.r.squared,
        summary(tcga_pdgfc_model_final)$adj.r.squared,
        summary(tcga_final_gleason_model)$adj.r.squared
      ),
      AIC = c(
        AIC(tcga_mylk_model_final),
        AIC(tcga_prkcb_model_final),
        AIC(tcga_pdgfc_model_final),
        AIC(tcga_final_gleason_model)
      )
    )
    
    model_comparison
    geo_mylk_model <- lm(
      gleason_score ~ MYLK + age_years,
      data = geo_gleason
    )
    
    geo_prkcb_model <- lm(
      gleason_score ~ PRKCB + age_years,
      data = geo_gleason
    )
    
    geo_pdgfc_model <- lm(
      gleason_score ~ PDGFC + age_years,
      data = geo_gleason
    )
    
    geo_model_comparison <- data.frame(
      Model = c(
        "MYLK",
        "PRKCB",
        "PDGFC",
        "3-gene signature"
      ),
      Beta = c(
        coef(geo_mylk_model)["MYLK"],
        coef(geo_prkcb_model)["PRKCB"],
        coef(geo_pdgfc_model)["PDGFC"],
        coef(geo_signature_model)["geo_signature"]
      ),
      P_value = c(
        summary(geo_mylk_model)$coefficients["MYLK", "Pr(>|t|)"],
        summary(geo_prkcb_model)$coefficients["PRKCB", "Pr(>|t|)"],
        summary(geo_pdgfc_model)$coefficients["PDGFC", "Pr(>|t|)"],
        summary(geo_signature_model)$coefficients["geo_signature", "Pr(>|t|)"]
      ),
      Adjusted_R2 = c(
        summary(geo_mylk_model)$adj.r.squared,
        summary(geo_prkcb_model)$adj.r.squared,
        summary(geo_pdgfc_model)$adj.r.squared,
        summary(geo_signature_model)$adj.r.squared
      ),
      AIC = c(
        AIC(geo_mylk_model),
        AIC(geo_prkcb_model),
        AIC(geo_pdgfc_model),
        AIC(geo_signature_model)
      )
    )
    
    geo_model_comparison
    summary(tcga_prkcb_model_final)
    summary(tcga_pdgfc_model_final)
    summary(geo_mylk_model)
    summary(geo_prkcb_model)
    summary(geo_pdgfc_model)
    geo_gleason_ordinal <- factor(
      geo_gleason$gleason_score,
      levels = sort(unique(geo_gleason$gleason_score)),
      ordered = TRUE
    )
    
    geo_signature_ordinal <- MASS::polr(
      geo_gleason_ordinal ~ geo_signature + age_years,
      data = geo_gleason,
      Hess = TRUE
    )
    
    summary(geo_signature_ordinal)
    geo_ord_coef <- coef(
      summary(geo_signature_ordinal)
    )
    
    geo_ord_p <- 2 * pnorm(
      abs(geo_ord_coef["geo_signature", "t value"]),
      lower.tail = FALSE
    )
    
    geo_ord_OR <- exp(
      coef(geo_signature_ordinal)["geo_signature"]
    )
    
    geo_ord_CI <- exp(
      confint(
        geo_signature_ordinal,
        parm = "geo_signature"
      )
    )
    
    geo_ord_OR
    geo_ord_CI
    geo_ord_p
    geo_signature_sd <- sd(
      geo_gleason$geo_signature,
      na.rm = TRUE
    )
    
    geo_gleason$geo_signature_z <- 
      geo_gleason$geo_signature / geo_signature_sd
    
    geo_signature_ordinal_z <- MASS::polr(
      geo_gleason_ordinal ~ geo_signature_z + age_years,
      data = geo_gleason,
      Hess = TRUE
    )
    
    summary(geo_signature_ordinal_z)
    colnames(geo_gleason)
    colnames(geo_ordinal)   
    head(geo_gleason)
    head(geo_ordinal)   
    # Add the existing GEO signature to geo_gleason
    geo_gleason$geo_signature <- geo_ordinal$signature[
      match(
        geo_gleason$geo_accession,
        geo_ordinal$geo_accession
      )
    ]
    
    # Check
    length(geo_gleason$geo_signature)
    sum(is.na(geo_gleason$geo_signature))
    cor(
      geo_gleason$geo_signature,
      geo_ordinal$signature[
        match(
          geo_gleason$geo_accession,
          geo_ordinal$geo_accession
        )
      ],
      method = "pearson"
    )
    head(geo_gleason$geo_signature)    
    geo_signature_sd <- sd(
      geo_gleason$geo_signature,
      na.rm = TRUE
    )
    
    geo_gleason$geo_signature_z <-
      geo_gleason$geo_signature / geo_signature_sd
    
    geo_gleason_ordinal <- factor(
      geo_gleason$gleason_score,
      levels = sort(unique(geo_gleason$gleason_score)),
      ordered = TRUE
    )
    
    geo_signature_ordinal_z <- MASS::polr(
      geo_gleason_ordinal ~ geo_signature_z + age_years,
      data = geo_gleason,
      Hess = TRUE
    )
    
    summary(geo_signature_ordinal_z)
    geo_z_coef <- coef(
      summary(geo_signature_ordinal_z)
    )
    
    geo_z_p <- 2 * pnorm(
      abs(
        geo_z_coef["geo_signature_z", "t value"]
      ),
      lower.tail = FALSE
    )
    
    geo_z_OR <- exp(
      coef(geo_signature_ordinal_z)["geo_signature_z"]
    )
    
    geo_z_CI <- exp(
      confint(
        geo_signature_ordinal_z,
        parm = "geo_signature_z"
      )
    )
    
    geo_signature_sd
    geo_z_OR
    geo_z_CI
    geo_z_p
    # ============================================
    # STEP 1 — Separate positive and negative GSEA
    # ============================================
    
    gsea_go_positive <- gsea_go_plot[
      gsea_go_plot$NES > 0 &
        gsea_go_plot$p.adjust < 0.05,
    ]
    
    gsea_go_negative <- gsea_go_plot[
      gsea_go_plot$NES < 0 &
        gsea_go_plot$p.adjust < 0.05,
    ]
    
    gsea_kegg_positive <- gsea_kegg_plot[
      gsea_kegg_plot$NES > 0 &
        gsea_kegg_plot$p.adjust < 0.05,
    ]
    
    gsea_kegg_negative <- gsea_kegg_plot[
      gsea_kegg_plot$NES < 0 &
        gsea_kegg_plot$p.adjust < 0.05,
    ]
    # Top positive GO
    head(
      gsea_go_positive[
        order(gsea_go_positive$p.adjust),
        c("Description", "NES", "p.adjust")
      ],
      10
    )
    
    # Top negative GO
    head(
      gsea_go_negative[
        order(gsea_go_negative$p.adjust),
        c("Description", "NES", "p.adjust")
      ],
      10
    )
    
    # Top positive KEGG
    head(
      gsea_kegg_positive[
        order(gsea_kegg_positive$p.adjust),
        c("Description", "NES", "p.adjust")
      ],
      10
    )
    
    # Top negative KEGG
    head(
      gsea_kegg_negative[
        order(gsea_kegg_negative$p.adjust),
        c("Description", "NES", "p.adjust")
      ],
      10
    )
    # ============================================
    # SIGNATURE GENE CONTRIBUTION TO GSEA
    # ============================================
    
    signature_genes <- c(
      "MYLK",
      "PRKCB",
      "PDGFC"
    )
    
    signature_genes
    signature_genes %in% names(gene_list)
    head(names(gene_list))
    tail(names(gene_list))    
    library(org.Hs.eg.db)
    
    signature_entrez <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = signature_genes,
      columns = c("SYMBOL", "ENTREZID"),
      keytype = "SYMBOL"
    )
    
    signature_entrez
    signature_entrez$ENTREZID %in% names(gene_list)
    gene_list[
      names(gene_list) %in% signature_entrez$ENTREZID
    ]    
    match(
      signature_entrez$ENTREZID,
      names(gene_list)
    )
    signature_rank_check <- data.frame(
      SYMBOL = signature_entrez$SYMBOL,
      ENTREZID = signature_entrez$ENTREZID,
      Rank = match(
        signature_entrez$ENTREZID,
        names(gene_list)
      ),
      Statistic = gene_list[
        match(
          signature_entrez$ENTREZID,
          names(gene_list)
        )
      ]
    )
    
    signature_rank_check
    gene_list_no_signature <- gene_list[
      !names(gene_list) %in% signature_entrez$ENTREZID
    ]
    
    length(gene_list)
    length(gene_list_no_signature)
    gsea_go_no_signature <- gseGO(
      geneList = gene_list_no_signature,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
    gsea_go_no_signature_plot <- as.data.frame(
      gsea_go_no_signature
    )
    gsea_go_no_signature_plot[
      grepl(
        "ribosome|rRNA|translation|ribonucleoprotein",
        gsea_go_no_signature_plot$Description,
        ignore.case = TRUE
      ),
      c("Description", "NES", "p.adjust")
    ]
    gsea_kegg_no_signature <- gseKEGG(
      geneList = gene_list_no_signature,
      organism = "hsa",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
    
    gsea_kegg_no_signature_plot <- as.data.frame(
      gsea_kegg_no_signature
    )
    
    gsea_kegg_no_signature_plot[
      grepl(
        "ribosome|translation|oxidative phosphorylation|proteasome",
        gsea_kegg_no_signature_plot$Description,
        ignore.case = TRUE
      ),
      c("Description", "NES", "p.adjust")
    ]
    gsea_go_no_signature_plot[
      gsea_go_no_signature_plot$NES < 0 &
        gsea_go_no_signature_plot$p.adjust < 0.05 &
        grepl(
          "adhesion|cytoskeleton|motility|migration|locomotion|contraction|wound",
          gsea_go_no_signature_plot$Description,
          ignore.case = TRUE
        ),
      c("Description", "NES", "p.adjust")
    ]
    gsea_kegg_no_signature_plot[
      gsea_kegg_no_signature_plot$NES < 0 &
        gsea_kegg_no_signature_plot$p.adjust < 0.05 &
        grepl(
          "adhesion|cytoskeleton|calcium|axon|proteoglycan|motility",
          gsea_kegg_no_signature_plot$Description,
          ignore.case = TRUE
        ),
      c("Description", "NES", "p.adjust")
    ]
    library(ggplot2)
    library(dplyr)
    
    # Make sure Gleason is ordered
    fig3_gleason$Gleason <- factor(
      fig3_gleason$Gleason,
      levels = c(6, 7, 8, 9, 10),
      ordered = TRUE
    )
    
    # Plot
    p3A <- ggplot(
      fig3_gleason,
      aes(x = Gleason, y = signature)
    ) +
      geom_boxplot(
        width = 0.65,
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.35,
        size = 1.2
      ) +
      facet_wrap(~ Cohort) +
      labs(
        x = "Gleason score",
        y = "3-gene signature score"
      ) +
      theme_classic(base_size = 13)
    
    p3A
    library(ggplot2)
    library(dplyr)
    
    # Sample sizes
    gleason_n <- fig3_gleason %>%
      count(Cohort, Gleason)
    
    # Annotation positions
    gleason_annotation <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(6.15, 6.15),
      y = c(-2.0, 3.25),
      label = c(
        "Kruskal–Wallis P = 0.0176",
        "Kruskal–Wallis P < 0.001"
      )
    )
    
    p3A_final <- ggplot(
      fig3_gleason,
      aes(x = Gleason, y = signature)
    ) +
      geom_boxplot(
        width = 0.62,
        outlier.shape = NA,
        linewidth = 0.7
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.35,
        size = 1.3
      ) +
      facet_wrap(~ Cohort) +
      geom_text(
        data = gleason_annotation,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        fontface = "bold",
        size = 4.2
      ) +
      scale_x_discrete(drop = FALSE) +
      labs(
        x = "Gleason score",
        y = "MYLK–PRKCB–PDGFC signature"
      ) +
      theme_classic(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 16
        ),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")
      )
    
    p3A_final
    # ================================
    # FIGURE 3B — SIGNATURE vs T STAGE
    # ================================
    
    T_annotation <- data.frame(
      Cohort = c("GEO", "TCGA"),
      x = c(1, 1),
      y = c(-2.0, 3.25),
      label = c(
        "Wilcoxon P = 0.296",
        "Wilcoxon P = 0.052"
      )
    )
    
    p3B_final <- ggplot(
      fig3_T,
      aes(x = T_stage, y = signature)
    ) +
      geom_boxplot(
        width = 0.60,
        outlier.shape = NA,
        linewidth = 0.7
      ) +
      geom_jitter(
        width = 0.12,
        alpha = 0.35,
        size = 1.3
      ) +
      facet_wrap(~ Cohort) +
      geom_text(
        data = T_annotation,
        aes(
          x = x,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        fontface = "bold",
        size = 4.2
      ) +
      labs(
        x = "Pathologic T stage",
        y = "MYLK–PRKCB–PDGFC signature"
      ) +
      theme_classic(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold",
          size = 16
        ),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")
      )
    
    p3B_final
    # ==========================================
    # FIGURE 3C — MULTIVARIABLE ORDINAL MODEL
    # ==========================================
    
    panel_C <- data.frame(
      Variable = "3-gene signature\n(per 1-SD increase)",
      OR = geo_z_OR,
      CI_lower = geo_z_CI[1],
      CI_upper = geo_z_CI[2],
      P_value = geo_z_p
    )
    
    panel_C
    p3C <- ggplot(
      panel_C,
      aes(
        x = OR,
        y = Variable
      )
    ) +
      geom_vline(
        xintercept = 1,
        linetype = "dashed",
        linewidth = 0.7
      ) +
      geom_errorbarh(
        aes(
          xmin = CI_lower,
          xmax = CI_upper
        ),
        height = 0.15,
        linewidth = 1
      ) +
      geom_point(
        size = 4
      ) +
      scale_x_log10(
        limits = c(0.1, 2),
        breaks = c(0.1, 0.25, 0.5, 1, 2)
      ) +
      labs(
        x = "Odds ratio (95% CI)",
        y = NULL,
        title = "Multivariable ordinal regression — GEO"
      ) +
      annotate(
        "text",
        x = 1.45,
        y = 1,
        label = "OR = 0.445\n95% CI 0.274–0.703\nP < 0.001",
        hjust = 0,
        fontface = "bold",
        size = 4
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.title.x = element_text(face = "bold")
      )
    
    p3C
    # ==========================================
    # FIGURE 3D — GSEA SENSITIVITY ANALYSIS
    # ==========================================
    
    library(dplyr)
    library(ggplot2)
    
    # ---- Top positive GO ----
    go_pos <- as.data.frame(gsea_go_no_signature) %>%
      filter(NES > 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = 8) %>%
      mutate(
        Direction = "Positive enrichment"
      )
    
    # ---- Top negative GO ----
    go_neg <- as.data.frame(gsea_go_no_signature) %>%
      filter(NES < 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = 8) %>%
      mutate(
        Direction = "Negative enrichment"
      )
    
    # Combine
    go_plot_final <- bind_rows(
      go_pos,
      go_neg
    ) %>%
      mutate(
        Description = factor(
          Description,
          levels = Description[order(NES)]
        ),
        neg_log10_FDR = -log10(p.adjust)
      )
    
    # Plot
    p3D_GO <- ggplot(
      go_plot_final,
      aes(
        x = NES,
        y = Description,
        size = neg_log10_FDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        size = expression(-log[10](FDR)),
        title = "GO GSEA after removal of signature genes"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.title.x = element_text(
          face = "bold"
        )
      )
    
    p3D_GO
    # ==========================================
    # KEGG
    # ==========================================
    
    kegg_pos <- as.data.frame(gsea_kegg_no_signature) %>%
      filter(NES > 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = 8)
    
    kegg_neg <- as.data.frame(gsea_kegg_no_signature) %>%
      filter(NES < 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = 8)
    
    kegg_plot_final <- bind_rows(
      kegg_pos,
      kegg_neg
    ) %>%
      mutate(
        Description = factor(
          Description,
          levels = Description[order(NES)]
        ),
        neg_log10_FDR = -log10(p.adjust)
      )
    
    p3D_KEGG <- ggplot(
      kegg_plot_final,
      aes(
        x = NES,
        y = Description,
        size = neg_log10_FDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        size = expression(-log[10](FDR)),
        title = "KEGG GSEA after removal of signature genes"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.title.x = element_text(
          face = "bold"
        )
      )
    
    p3D_KEGG
    # ==========================================
    # SAVE GSEA SENSITIVITY PANELS
    # ==========================================
    
    ggsave(
      "Figure3D_GO_GSEA_sensitivity.pdf",
      p3D_GO,
      width = 8,
      height = 7,
      units = "in"
    )
    
    ggsave(
      "Figure3D_KEGG_GSEA_sensitivity.pdf",
      p3D_KEGG,
      width = 8,
      height = 7,
      units = "in"
    )
    
    ggsave(
      "Figure3D_GO_GSEA_sensitivity.png",
      p3D_GO,
      width = 8,
      height = 7,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      "Figure3D_KEGG_GSEA_sensitivity.png",
      p3D_KEGG,
      width = 8,
      height = 7,
      units = "in",
      dpi = 600
    )
    # ==========================================
    # SENSITIVITY ANALYSIS — KEY PATHWAYS
    # ==========================================
    
    go_sensitivity_summary <- as.data.frame(
      gsea_go_no_signature
    ) %>%
      filter(
        grepl(
          "ribosome|rRNA|translation|ribonucleoprotein|muscle|contractile|actin|sarcomere|adhesion",
          Description,
          ignore.case = TRUE
        )
      ) %>%
      select(
        Description,
        NES,
        p.adjust
      ) %>%
      arrange(p.adjust)
    
    go_sensitivity_summary
    go_sensitivity_summary <- as.data.frame(
      gsea_go_no_signature
    ) %>%
      dplyr::filter(
        grepl(
          "ribosome|rRNA|translation|ribonucleoprotein|muscle|contractile|actin|sarcomere|adhesion",
          Description,
          ignore.case = TRUE
        )
      ) %>%
      dplyr::select(
        Description,
        NES,
        p.adjust
      ) %>%
      dplyr::arrange(p.adjust)
    
    go_sensitivity_summary
    kegg_sensitivity_summary <- as.data.frame(
      gsea_kegg_no_signature
    ) %>%
      dplyr::filter(
        grepl(
          "ribosome|translation|oxidative phosphorylation|proteasome|calcium signaling|focal adhesion|cytoskeleton|cardiomyopathy|proteoglycans",
          Description,
          ignore.case = TRUE
        )
      ) %>%
      dplyr::select(
        Description,
        NES,
        p.adjust
      ) %>%
      dplyr::arrange(p.adjust)
    
    kegg_sensitivity_summary
    # ==========================================
    # COMPARE ORIGINAL vs SENSITIVITY GSEA
    # ==========================================
    
    go_original_df <- as.data.frame(gsea_go)
    
    go_sensitivity_df <- as.data.frame(
      gsea_go_no_signature
    )
    
    # Match pathways present in both analyses
    go_common <- intersect(
      rownames(go_original_df),
      rownames(go_sensitivity_df)
    )
    
    go_gsea_comparison <- data.frame(
      GO_ID = go_common,
      Original_NES = go_original_df[go_common, "NES"],
      Sensitivity_NES = go_sensitivity_df[go_common, "NES"]
    )
    
    # Correlation of pathway NES values
    cor(
      go_gsea_comparison$Original_NES,
      go_gsea_comparison$Sensitivity_NES,
      method = "spearman"
    )
    kegg_original_df <- as.data.frame(gsea_kegg)
    
    kegg_sensitivity_df <- as.data.frame(
      gsea_kegg_no_signature
    )
    
    kegg_common <- intersect(
      rownames(kegg_original_df),
      rownames(kegg_sensitivity_df)
    )
    
    kegg_gsea_comparison <- data.frame(
      KEGG_ID = kegg_common,
      Original_NES = kegg_original_df[kegg_common, "NES"],
      Sensitivity_NES = kegg_sensitivity_df[kegg_common, "NES"]
    )
    
    cor(
      kegg_gsea_comparison$Original_NES,
      kegg_gsea_comparison$Sensitivity_NES,
      method = "spearman"
    )
    library(ggplot2)
    
    # ============================================================
    # FIGURE 4A — GO GSEA sensitivity analysis
    # ============================================================
    
    go_original_df <- as.data.frame(gsea_go)
    go_sensitivity_df <- as.data.frame(gsea_go_no_signature)
    
    go_common <- intersect(
      rownames(go_original_df),
      rownames(go_sensitivity_df)
    )
    
    go_gsea_comparison <- data.frame(
      GO_ID = go_common,
      Description = go_original_df[go_common, "Description"],
      Original_NES = go_original_df[go_common, "NES"],
      Sensitivity_NES = go_sensitivity_df[go_common, "NES"]
    )
    
    # Correlation
    go_rho <- cor(
      go_gsea_comparison$Original_NES,
      go_gsea_comparison$Sensitivity_NES,
      method = "spearman"
    )
    
    go_rho
    ggplot(
      go_gsea_comparison,
      aes(
        x = Original_NES,
        y = Sensitivity_NES
      )
    ) +
      geom_point(size = 2.5, alpha = 0.65) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0(
          "Spearman \u03c1 = ",
          sprintf("%.3f", go_rho)
        ),
        hjust = -0.05,
        vjust = 1.5,
        size = 5
      ) +
      labs(
        title = "GO GSEA sensitivity analysis",
        subtitle = "Original signature vs. after removal of MYLK, PRKCB and PDGFC",
        x = "Original NES",
        y = "NES after signature-gene removal"
      ) +
      theme_classic(base_size = 14)
    ggplot(
      kegg_gsea_comparison,
      aes(
        x = Original_NES,
        y = Sensitivity_NES
      )
    ) +
      geom_point(
        size = 3,
        alpha = 0.75
      ) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0(
          "Spearman \u03c1 = ",
          sprintf("%.3f", kegg_rho)
        ),
        hjust = -0.05,
        vjust = 1.5,
        size = 5
      ) +
      labs(
        title = "KEGG GSEA sensitivity analysis",
        subtitle = "Original signature vs. after removal of MYLK, PRKCB and PDGFC",
        x = "Original NES",
        y = "NES after signature-gene removal"
      ) +
      theme_classic(base_size = 14)
    kegg_original_df <- as.data.frame(gsea_kegg)
    
    kegg_sensitivity_df <- as.data.frame(
      gsea_kegg_no_signature
    )
    
    kegg_common <- intersect(
      rownames(kegg_original_df),
      rownames(kegg_sensitivity_df)
    )
    
    kegg_gsea_comparison <- data.frame(
      KEGG_ID = kegg_common,
      Description = kegg_original_df[
        kegg_common, "Description"
      ],
      Original_NES = kegg_original_df[
        kegg_common, "NES"
      ],
      Sensitivity_NES = kegg_sensitivity_df[
        kegg_common, "NES"
      ]
    )
    
    kegg_rho <- cor(
      kegg_gsea_comparison$Original_NES,
      kegg_gsea_comparison$Sensitivity_NES,
      method = "spearman"
    )
    
    kegg_rho
    ggplot(
      kegg_gsea_comparison,
      aes(
        x = Original_NES,
        y = Sensitivity_NES
      )
    ) +
      geom_point(
        size = 3,
        alpha = 0.75
      ) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0(
          "Spearman \u03c1 = ",
          sprintf("%.3f", kegg_rho)
        ),
        hjust = -0.05,
        vjust = 1.5,
        size = 5
      ) +
      labs(
        title = "KEGG GSEA sensitivity analysis",
        subtitle = "Original signature vs. after removal of MYLK, PRKCB and PDGFC",
        x = "Original NES",
        y = "NES after signature-gene removal"
      ) +
      theme_classic(base_size = 14)
    library(ggplot2)
    library(patchwork)
    
    # ================================
    # FIGURE 4A — GO
    # ================================
    
    p_go_sensitivity <- ggplot(
      go_gsea_comparison,
      aes(
        x = Original_NES,
        y = Sensitivity_NES
      )
    ) +
      geom_point(
        size = 2.5,
        alpha = 0.65
      ) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0(
          "Spearman \u03c1 = ",
          sprintf("%.3f", go_rho)
        ),
        hjust = -0.05,
        vjust = 1.5,
        size = 5
      ) +
      labs(
        title = "A  GO GSEA",
        x = "Original NES",
        y = "NES after signature-gene removal"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_blank()
      )
    
    
    # ================================
    # FIGURE 4B — KEGG
    # ================================
    
    p_kegg_sensitivity <- ggplot(
      kegg_gsea_comparison,
      aes(
        x = Original_NES,
        y = Sensitivity_NES
      )
    ) +
      geom_point(
        size = 2.8,
        alpha = 0.7
      ) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.8
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0(
          "Spearman \u03c1 = ",
          sprintf("%.3f", kegg_rho)
        ),
        hjust = -0.05,
        vjust = 1.5,
        size = 5
      ) +
      labs(
        title = "B  KEGG GSEA",
        x = "Original NES",
        y = "NES after signature-gene removal"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_blank()
      )
    figure4 <- p_go_sensitivity +
      p_kegg_sensitivity +
      plot_layout(ncol = 2)
    
    figure4
    ggsave(
      "Figure_4_GSEA_sensitivity.pdf",
      figure4,
      width = 12,
      height = 5.5
    )
    
    ggsave(
      "Figure_4_GSEA_sensitivity.png",
      figure4,
      width = 12,
      height = 5.5,
      dpi = 600
    )
    tcga_weights
    tcga_gene_means
    tcga_gene_sds    
    data.frame(
      Gene = names(tcga_weights),
      Weight = as.numeric(tcga_weights),
      Mean = as.numeric(tcga_gene_means[names(tcga_weights)]),
      SD = as.numeric(tcga_gene_sds[names(tcga_weights)])
    )
    signature_direction_check <- data.frame(
      Gene = c("MYLK", "PRKCB", "PDGFC"),
      Spearman_rho = c(
        cor(
          tcga_signature_497,
          tcga_signature_expr[, "MYLK"],
          method = "spearman"
        ),
        cor(
          tcga_signature_497,
          tcga_signature_expr[, "PRKCB"],
          method = "spearman"
        ),
        cor(
          tcga_signature_497,
          tcga_signature_expr[, "PDGFC"],
          method = "spearman"
        )
      )
    )
    
    signature_direction_check
    tcga_signature_reconstructed <- as.numeric(
      scale(
        tcga_signature_expr[, c("MYLK", "PRKCB", "PDGFC")],
        center = tcga_gene_means[c("MYLK", "PRKCB", "PDGFC")],
        scale = tcga_gene_sds[c("MYLK", "PRKCB", "PDGFC")]
      ) %*%
        tcga_weights[c("MYLK", "PRKCB", "PDGFC")]
    )
    
    cor(
      tcga_signature_reconstructed,
      tcga_signature_497,
      method = "pearson"
    )
    # ============================================================
    # FINAL 3-GENE SIGNATURE DEFINITION
    # ============================================================
    
    signature_definition <- data.frame(
      Gene = names(tcga_weights),
      Weight = as.numeric(tcga_weights),
      Mean = as.numeric(tcga_gene_means[names(tcga_weights)]),
      SD = as.numeric(tcga_gene_sds[names(tcga_weights)])
    )
    
    signature_definition
    signature_formula <- paste0(
      "Signature = ",
      round(tcga_weights["MYLK"], 7),
      "*Z(MYLK) + ",
      round(tcga_weights["PRKCB"], 7),
      "*Z(PRKCB) + ",
      round(tcga_weights["PDGFC"], 7),
      "*Z(PDGFC)"
    )
    
    signature_formula
    levels(geo_gleason_ordinal)
    table(geo_gleason_ordinal)   
    table(geo_gleason_ordinal)    
    coef(geo_signature_ordinal_z)    
    # Check proportional-odds assumption
    library(brant)
    
    brant::brant(geo_signature_ordinal_z)    
    # ==========================================
    # PROPORTIONAL-ODDS CHECK — GEO
    # ==========================================
    
    geo_po_data <- geo_gleason
    
    # Cutoff 1: 6 vs 7-9
    geo_po_data$cut1 <- ifelse(
      geo_po_data$gleason_score >= 7,
      1,
      0
    )
    
    # Cutoff 2: 6-7 vs 8-9
    geo_po_data$cut2 <- ifelse(
      geo_po_data$gleason_score >= 8,
      1,
      0
    )
    
    # Cutoff 3: 6-8 vs 9
    geo_po_data$cut3 <- ifelse(
      geo_po_data$gleason_score >= 9,
      1,
      0
    )
    
    po_model_1 <- glm(
      cut1 ~ geo_signature_z + age_years,
      data = geo_po_data,
      family = binomial
    )
    
    po_model_2 <- glm(
      cut2 ~ geo_signature_z + age_years,
      data = geo_po_data,
      family = binomial
    )
    
    po_model_3 <- glm(
      cut3 ~ geo_signature_z + age_years,
      data = geo_po_data,
      family = binomial
    )
    
    po_results <- data.frame(
      Cutoff = c(
        "6 vs 7-9",
        "6-7 vs 8-9",
        "6-8 vs 9"
      ),
      Signature_Beta = c(
        coef(po_model_1)["geo_signature_z"],
        coef(po_model_2)["geo_signature_z"],
        coef(po_model_3)["geo_signature_z"]
      ),
      Signature_OR = c(
        exp(coef(po_model_1)["geo_signature_z"]),
        exp(coef(po_model_2)["geo_signature_z"]),
        exp(coef(po_model_3)["geo_signature_z"])
      ),
      P_value = c(
        summary(po_model_1)$coefficients[
          "geo_signature_z", "Pr(>|z|)"
        ],
        summary(po_model_2)$coefficients[
          "geo_signature_z", "Pr(>|z|)"
        ],
        summary(po_model_3)$coefficients[
          "geo_signature_z", "Pr(>|z|)"
        ]
      )
    )
    
    po_results
    # ==========================================
    # GEO GLEASON SENSITIVITY ANALYSIS
    # ==========================================
    
    geo_sens <- geo_gleason
    
    # 6 vs >=7
    geo_sens$gleason_high7 <- ifelse(
      geo_sens$gleason_score >= 7,
      1,
      0
    )
    
    # <=7 vs >=8
    geo_sens$gleason_high8 <- ifelse(
      geo_sens$gleason_score >= 8,
      1,
      0
    )
    
    # Logistic regression: 6 vs >=7
    geo_binary_7 <- glm(
      gleason_high7 ~ geo_signature_z + age_years,
      data = geo_sens,
      family = binomial
    )
    
    # Logistic regression: <=7 vs >=8
    geo_binary_8 <- glm(
      gleason_high8 ~ geo_signature_z + age_years,
      data = geo_sens,
      family = binomial
    )
    
    summary(geo_binary_7)
    summary(geo_binary_8)
    geo_binary_results <- data.frame(
      Comparison = c(
        "Gleason 6 vs >=7",
        "Gleason <=7 vs >=8"
      ),
      OR = c(
        exp(coef(geo_binary_7)["geo_signature_z"]),
        exp(coef(geo_binary_8)["geo_signature_z"])
      ),
      P_value = c(
        summary(geo_binary_7)$coefficients[
          "geo_signature_z", "Pr(>|z|)"
        ],
        summary(geo_binary_8)$coefficients[
          "geo_signature_z", "Pr(>|z|)"
        ]
      )
    )
    
    geo_binary_results
    ls()
    ls(pattern = "fig|plot|gleason|model|signature",
       ignore.case = TRUE)    
    objects <- ls()
    
    objects[
      grepl(
        "fig|plot|gleason|model|signature",
        objects,
        ignore.case = TRUE
      )
    ]
    # ============================================================
    # INVENTORY OF EXISTING FIGURE OBJECTS
    # ============================================================
    
    figure_objects <- c(
      "fig3_gleason",
      "fig3_T",
      "figure4",
      "final_GSEA_figure",
      "geo_gleason_plot",
      "geo_plot",
      "geo_T_plot",
      "gleason_plot",
      "go_fig",
      "go_plot_final",
      "gsea_down_plot",
      "gsea_up_plot",
      "kegg_fig",
      "kegg_plot_final",
      "p_fig3",
      "p_gleason_final",
      "tcga_gleason_plot",
      "tcga_plot",
      "tcga_T_plot"
    )
    
    figure_inventory <- data.frame(
      Object = figure_objects,
      Exists = sapply(
        figure_objects,
        exists
      ),
      Class = sapply(
        figure_objects,
        function(x) {
          if (exists(x)) {
            paste(class(get(x)), collapse = ", ")
          } else {
            NA_character_
          }
        }
      )
    )
    
    figure_inventory
    figure_inventory[
      grepl(
        "ggplot",
        figure_inventory$Class,
        ignore.case = TRUE
      ),
    ]
    data_inventory <- c(
      "fig3_gleason",
      "fig3_T",
      "geo_gleason",
      "gsea_go_plot",
      "gsea_kegg_plot",
      "gsea_go_no_signature_plot",
      "gsea_kegg_no_signature_plot",
      "tcga_survival_signature",
      "tcga_survival_complete",
      "tcga_t_model"
    )
    
    data.frame(
      Object = data_inventory,
      Exists = sapply(data_inventory, exists),
      Rows = sapply(
        data_inventory,
        function(x) {
          if (exists(x) && !is.null(dim(get(x)))) {
            nrow(get(x))
          } else {
            NA
          }
        }
      ),
      Columns = sapply(
        data_inventory,
        function(x) {
          if (exists(x) && !is.null(dim(get(x)))) {
            ncol(get(x))
          } else {
            NA
          }
        }
      )
    )
    p_fig3
    p_gleason_final   
    final_GSEA_figure 
    ############################################################
    # FINAL STATISTICAL AUDIT — 3-GENE SIGNATURE
    ############################################################
    
    cat("\n==============================\n")
    cat("1. SIGNATURE DEFINITION\n")
    cat("==============================\n")
    
    print(signature_definition)
    cat("\nFormula:\n")
    print(signature_formula)
    
    
    ############################################################
    # 2. GEO ORDINAL REGRESSION
    ############################################################
    
    cat("\n==============================\n")
    cat("2. GEO ORDINAL REGRESSION\n")
    cat("==============================\n")
    
    geo_ord_sum <- summary(geo_signature_ordinal_z)
    
    print(geo_ord_sum)
    
    geo_coef <- coef(geo_ord_sum)
    
    geo_beta <- geo_coef["geo_signature_z", "Value"]
    geo_se   <- geo_coef["geo_signature_z", "Std. Error"]
    geo_t    <- geo_coef["geo_signature_z", "t value"]
    
    geo_p <- 2 * pnorm(abs(geo_t), lower.tail = FALSE)
    geo_OR <- exp(geo_beta)
    
    geo_CI <- exp(
      confint(
        geo_signature_ordinal_z,
        parm = "geo_signature_z"
      )
    )
    
    cat("\nGEO signature:\n")
    cat("Beta =", geo_beta, "\n")
    cat("SE   =", geo_se, "\n")
    cat("OR   =", geo_OR, "\n")
    cat("95% CI =", geo_CI[1], "-", geo_CI[2], "\n")
    cat("P    =", geo_p, "\n")
    
    
    ############################################################
    # 3. GEO BINARY SENSITIVITY
    ############################################################
    
    cat("\n==============================\n")
    cat("3. GEO BINARY SENSITIVITY\n")
    cat("==============================\n")
    
    print(geo_binary_results)
    
    
    ############################################################
    # 4. GEO GENE-LEVEL MODELS
    ############################################################
    
    cat("\n==============================\n")
    cat("4. GEO GENE-LEVEL MODELS\n")
    cat("==============================\n")
    
    print(geo_model_comparison)
    
    
    ############################################################
    # 5. TCGA GENE-LEVEL MODELS
    ############################################################
    
    cat("\n==============================\n")
    cat("5. TCGA GENE-LEVEL MODELS\n")
    cat("==============================\n")
    
    print(model_comparison)
    
    
    ############################################################
    # 6. TCGA SIGNATURE ASSOCIATION WITH GLEASON
    ############################################################
    
    cat("\n==============================\n")
    cat("6. TCGA SIGNATURE MODEL\n")
    cat("==============================\n")
    
    if (exists("tcga_signature_ordinal")) {
      print(summary(tcga_signature_ordinal))
    }
    
    if (exists("tcga_signature_ordinal_no10")) {
      cat("\nTCGA sensitivity model excluding Gleason 10:\n")
      print(summary(tcga_signature_ordinal_no10))
    }
    
    
    ############################################################
    # 7. TCGA T-STAGE
    ############################################################
    
    cat("\n==============================\n")
    cat("7. T-STAGE\n")
    cat("==============================\n")
    
    if (exists("tcga_signature_T_model")) {
      print(summary(tcga_signature_T_model))
    }
    
    if (exists("geo_signature_T_model")) {
      cat("\nGEO T-stage model:\n")
      print(summary(geo_signature_T_model))
    }
    
    
    ############################################################
    # 8. TCGA SURVIVAL
    ############################################################
    
    cat("\n==============================\n")
    cat("8. SURVIVAL ANALYSIS\n")
    cat("==============================\n")
    
    if (exists("tcga_signature_cox")) {
      cat("\nUnadjusted Cox model:\n")
      print(summary(tcga_signature_cox))
    }
    
    if (exists("tcga_signature_cox_adj")) {
      cat("\nAdjusted Cox model:\n")
      print(summary(tcga_signature_cox_adj))
    }
    
    if (exists("tcga_signature_cox_z")) {
      cat("\nZ-scaled Cox model:\n")
      print(summary(tcga_signature_cox_z))
    }
    
    
    ############################################################
    # 9. GSEA ROBUSTNESS
    ############################################################
    
    cat("\n==============================\n")
    cat("9. GSEA SENSITIVITY\n")
    cat("==============================\n")
    
    go_original_df <- as.data.frame(gsea_go)
    go_sensitivity_df <- as.data.frame(gsea_go_no_signature)
    
    go_common <- intersect(
      rownames(go_original_df),
      rownames(go_sensitivity_df)
    )
    
    go_gsea_comparison <- data.frame(
      GO_ID = go_common,
      Original_NES =
        go_original_df[go_common, "NES"],
      Sensitivity_NES =
        go_sensitivity_df[go_common, "NES"]
    )
    
    go_rho <- cor(
      go_gsea_comparison$Original_NES,
      go_gsea_comparison$Sensitivity_NES,
      method = "spearman"
    )
    
    cat("GO Spearman rho =", go_rho, "\n")
    
    
    kegg_original_df <- as.data.frame(gsea_kegg)
    kegg_sensitivity_df <- as.data.frame(gsea_kegg_no_signature)
    
    kegg_common <- intersect(
      rownames(kegg_original_df),
      rownames(kegg_sensitivity_df)
    )
    
    kegg_gsea_comparison <- data.frame(
      KEGG_ID = kegg_common,
      Original_NES =
        kegg_original_df[kegg_common, "NES"],
      Sensitivity_NES =
        kegg_sensitivity_df[kegg_common, "NES"]
    )
    
    kegg_rho <- cor(
      kegg_gsea_comparison$Original_NES,
      kegg_gsea_comparison$Sensitivity_NES,
      method = "spearman"
    )
    
    cat("KEGG Spearman rho =", kegg_rho, "\n")
    
    
    ############################################################
    # 10. SAMPLE COUNTS
    ############################################################
    
    cat("\n==============================\n")
    cat("10. SAMPLE COUNTS\n")
    cat("==============================\n")
    
    if (exists("geo_gleason")) {
      cat("GEO N =", nrow(geo_gleason), "\n")
      print(table(geo_gleason_ordinal))
    }
    
    if (exists("tcga_t_model")) {
      cat("TCGA N =", nrow(tcga_t_model), "\n")
    }
    
    if (exists("tcga_survival_signature")) {
      cat(
        "TCGA survival signature N =",
        nrow(tcga_survival_signature),
        "\n"
      )
    }
    
    if (exists("tcga_survival_complete")) {
      cat(
        "TCGA complete survival N =",
        nrow(tcga_survival_complete),
        "\n"
      )
    }
    
    
    ############################################################
    # 11. PCA RELATIONSHIP
    ############################################################
    
    cat("\n==============================\n")
    cat("11. PCA VALIDATION\n")
    cat("==============================\n")
    
    cat(
      "Pearson signature vs PC1 =",
      cor(
        tcga_signature_pca,
        tcga_pca$x[, 1],
        method = "pearson"
      ),
      "\n"
    )
    
    cat(
      "Spearman signature vs PC1 =",
      cor(
        tcga_signature_pca,
        tcga_pca$x[, 1],
        method = "spearman"
      ),
      "\n"
    )
    
    
    ############################################################
    # 12. SIGNATURE RECONSTRUCTION
    ############################################################
    
    cat("\n==============================\n")
    cat("12. SIGNATURE RECONSTRUCTION\n")
    cat("==============================\n")
    
    cat(
      "Correlation reconstructed vs original =",
      cor(
        tcga_signature_reconstructed,
        tcga_signature_497,
        method = "pearson"
      ),
      "\n"
    )
    
    
    ############################################################
    # 13. FINAL SUMMARY TABLE
    ############################################################
    
    cat("\n==============================\n")
    cat("13. KEY RESULTS\n")
    cat("==============================\n")
    
    key_results <- data.frame(
      Analysis = c(
        "GEO ordinal regression",
        "GEO Gleason 6 vs >=7",
        "GEO Gleason <=7 vs >=8",
        "GO GSEA sensitivity",
        "KEGG GSEA sensitivity",
        "Signature vs TCGA PC1"
      ),
      Effect = c(
        geo_OR,
        geo_binary_results$OR[1],
        geo_binary_results$OR[2],
        go_rho,
        kegg_rho,
        cor(
          tcga_signature_pca,
          tcga_pca$x[, 1],
          method = "spearman"
        )
      ),
      P_value = c(
        geo_p,
        geo_binary_results$P_value[1],
        geo_binary_results$P_value[2],
        NA,
        NA,
        NA
      )
    )
    
    print(key_results)    
    ############################################################
    # PROPORTIONAL ODDS / ORDINAL MODEL SENSITIVITY
    ############################################################
    
    # Compare the ordinal model with separate binary models
    # at each Gleason threshold.
    
    po_results
    tcga_coef <- coef(summary(tcga_signature_ordinal))
    
    tcga_beta <- tcga_coef["signature", "Value"]
    tcga_se   <- tcga_coef["signature", "Std. Error"]
    tcga_t    <- tcga_coef["signature", "t value"]
    
    tcga_p <- 2 * pnorm(
      abs(tcga_t),
      lower.tail = FALSE
    )
    
    tcga_OR <- exp(tcga_beta)
    
    tcga_CI_wald <- exp(
      c(
        tcga_beta - 1.96 * tcga_se,
        tcga_beta + 1.96 * tcga_se
      )
    )
    
    data.frame(
      Beta = tcga_beta,
      SE = tcga_se,
      OR = tcga_OR,
      CI_lower = tcga_CI_wald[1],
      CI_upper = tcga_CI_wald[2],
      P_value = tcga_p
    )
    ############################################################
    # TABLE 1 — COHORT CHARACTERISTICS
    ############################################################
    
    cat("========== TCGA ==========\n")
    
    cat("N =", nrow(tcga_t_model), "\n")
    
    summary(tcga_t_model[, c(
      "age_years",
      "gleason_score"
    )])
    
    cat("\nGleason distribution:\n")
    print(table(tcga_t_model$gleason_score))
    
    cat("\nT-stage distribution:\n")
    print(table(tcga_t_model$ajcc_pathologic_t, useNA = "ifany"))
    
    cat("\nN-stage distribution:\n")
    print(table(tcga_t_model$ajcc_pathologic_n, useNA = "ifany"))
    
    
    cat("\n========== GEO ==========\n")
    
    cat("N =", nrow(geo_gleason), "\n")
    
    summary(geo_gleason[, c(
      "age_years",
      "gleason_score"
    )])
    
    cat("\nGleason distribution:\n")
    print(table(geo_gleason$gleason_score))
    
    if ("ajcc_pathologic_t" %in% colnames(geo_gleason)) {
      cat("\nT-stage distribution:\n")
      print(table(
        geo_gleason$ajcc_pathologic_t,
        useNA = "ifany"
      ))
    }
    ############################################################
    # TABLE 2 — CLINICAL ASSOCIATION
    ############################################################
    
    table2 <- data.frame(
      Cohort = c(
        "GEO",
        "TCGA"
      ),
      N = c(
        nrow(geo_gleason),
        nrow(tcga_t_model)
      ),
      Beta = c(
        coef(geo_signature_ordinal_z)["geo_signature_z"],
        coef(tcga_signature_ordinal)["signature"]
      ),
      OR = c(
        exp(coef(geo_signature_ordinal_z)["geo_signature_z"]),
        exp(coef(tcga_signature_ordinal)["signature"])
      ),
      P_value = c(
        geo_z_p,
        tcga_p
      )
    )
    
    table2
    table2$CI_lower <- c(
      geo_z_CI[1],
      tcga_CI_wald[1]
    )
    
    table2$CI_upper <- c(
      geo_z_CI[2],
      tcga_CI_wald[2]
    )
    
    table2 <- table2[, c(
      "Cohort",
      "N",
      "Beta",
      "OR",
      "CI_lower",
      "CI_upper",
      "P_value"
    )]
    
    table2
    ############################################################
    # TABLE 3 — SENSITIVITY / ROBUSTNESS
    ############################################################
    
    table3 <- data.frame(
      Analysis = c(
        "GEO: Gleason 6 vs >=7",
        "GEO: Gleason <=7 vs >=8",
        "GEO: Gleason <=8 vs 9",
        "GO GSEA: signature genes removed",
        "KEGG GSEA: signature genes removed"
      ),
      Effect = c(
        geo_binary_results$OR[1],
        geo_binary_results$OR[2],
        geo_binary_results$OR[3],
        go_rho,
        kegg_rho
      ),
      P_value = c(
        geo_binary_results$P_value[1],
        geo_binary_results$P_value[2],
        geo_binary_results$P_value[3],
        NA,
        NA
      )
    )
    
    table3
    ############################################################
    # TABLE 4 — INDIVIDUAL GENES VS SIGNATURE
    ############################################################
    
    model_comparison
    geo_model_comparison
    ############################################################
    # FIGURE AUDIT — EXISTING MANUSCRIPT FIGURES
    ############################################################
    
    figure_objects <- c(
      "p_fig3",
      "p_gleason_final",
      "figure4",
      "final_GSEA_figure",
      "fig3_gleason",
      "fig3_T",
      "geo_gleason_plot",
      "geo_plot",
      "geo_T_plot",
      "tcga_gleason_plot",
      "tcga_plot",
      "tcga_T_plot",
      "go_plot_final",
      "kegg_plot_final"
    )
    
    figure_audit <- data.frame(
      Object = figure_objects,
      Exists = sapply(figure_objects, exists),
      Class = sapply(
        figure_objects,
        function(x) {
          if (exists(x)) {
            paste(class(get(x)), collapse = ", ")
          } else {
            NA_character_
          }
        }
      ),
      Dimensions = sapply(
        figure_objects,
        function(x) {
          if (exists(x) && !is.null(dim(get(x)))) {
            paste(dim(get(x)), collapse = " x ")
          } else {
            NA_character_
          }
        }
      )
    )
    
    figure_audit
    ############################################################
    # INSPECT CLINICAL FIGURE DATA
    ############################################################
    
    cat("========== fig3_gleason ==========\n")
    print(names(fig3_gleason))
    print(head(fig3_gleason))
    print(summary(fig3_gleason))
    
    cat("\n========== fig3_T ==========\n")
    print(names(fig3_T))
    print(head(fig3_T))
    print(summary(fig3_T))
    
    cat("\n========== geo_gleason_plot ==========\n")
    print(names(geo_gleason_plot))
    print(head(geo_gleason_plot))
    
    cat("\n========== tcga_gleason_plot ==========\n")
    print(names(tcga_gleason_plot))
    print(head(tcga_gleason_plot))
    ############################################################
    # INSPECT EXISTING COMPOSITE FIGURES
    ############################################################
    
    print(p_fig3)
    print(p_gleason_final)
    print(figure4)
    print(final_GSEA_figure)
    ############################################################
    # FIGURE 2 — INSPECT EXISTING GLEASON PLOTS
    ############################################################
    
    print(geo_plot)
    print(tcga_plot)
    
    cat("\n========== GEO PLOT DATA ==========\n")
    print(names(geo_plot))
    print(head(geo_plot))
    print(summary(geo_plot))
    
    cat("\n========== TCGA PLOT DATA ==========\n")
    print(names(tcga_plot))
    print(head(tcga_plot))
    print(summary(tcga_plot))
    ############################################################
    # INSPECT PLOT LAYERS
    ############################################################
    
    cat("\n========== p_gleason_final ==========\n")
    print(p_gleason_final$layers)
    
    cat("\n========== p_fig3 ==========\n")
    print(p_fig3$layers)
    ############################################################
    # CHECK SCALE / AXIS DEFINITIONS
    ############################################################
    
    cat("\n========== p_gleason_final SCALES ==========\n")
    print(p_gleason_final$scales)
    
    cat("\n========== p_fig3 SCALES ==========\n")
    print(p_fig3$scales)
    ############################################################
    # FINAL FIGURE 2 — SIGNATURE VS GLEASON
    ############################################################
    
    library(ggplot2)
    library(patchwork)
    
    # GEO
    geo_fig2 <- ggplot(
      geo_plot,
      aes(
        x = signature_z,
        y = Gleason
      )
    ) +
      geom_jitter(
        width = 0.06,
        height = 0.08,
        alpha = 0.55,
        size = 1.8
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0(
          "Ordinal logistic regression\n",
          "OR = ", round(geo_z_OR, 2),
          " (95% CI ",
          round(geo_z_CI[1], 2), "–",
          round(geo_z_CI[2], 2), ")\n",
          "P = ", format.pval(geo_z_p, digits = 3, eps = 0.001)
        ),
        hjust = 1.05,
        vjust = 1.2,
        size = 3.5
      ) +
      scale_y_continuous(
        breaks = sort(unique(geo_plot$Gleason))
      ) +
      labs(
        title = "GEO validation cohort",
        x = "3-gene signature (standardized)",
        y = "Gleason score"
      ) +
      theme_classic(base_size = 12)
    
    
    # TCGA
    tcga_fig2 <- ggplot(
      tcga_plot,
      aes(
        x = signature_z,
        y = Gleason
      )
    ) +
      geom_jitter(
        width = 0.06,
        height = 0.08,
        alpha = 0.40,
        size = 1.6
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0(
          "Ordinal logistic regression\n",
          "OR = ", round(tcga_OR, 2),
          " (95% CI ",
          round(tcga_CI_wald[1], 2), "–",
          round(tcga_CI_wald[2], 2), ")\n",
          "P = ", format.pval(tcga_p, digits = 3, eps = 0.001)
        ),
        hjust = 1.05,
        vjust = 1.2,
        size = 3.5
      ) +
      scale_y_continuous(
        breaks = sort(unique(tcga_plot$Gleason))
      ) +
      labs(
        title = "TCGA validation cohort",
        x = "3-gene signature (standardized)",
        y = "Gleason score"
      ) +
      theme_classic(base_size = 12)
    
    
    ############################################################
    # COMBINE
    ############################################################
    
    p_gleason_final_v2 <-
      geo_fig2 + tcga_fig2 +
      plot_annotation(
        title = "Association of the 3-gene signature with Gleason score"
      )
    
    print(p_gleason_final_v2)
    ############################################################
    # FIGURE 3 — T STAGE ASSOCIATION
    ############################################################
    
    cat("========== GEO T-STAGE MODEL ==========\n")
    summary(t_stage_model)
    
    cat("\n========== TCGA T-STAGE MODEL ==========\n")
    summary(tcga_t_model)
    
    cat("\n========== GEO T-STAGE DATA ==========\n")
    print(table(geo_T_plot$T_stage))
    
    cat("\n========== TCGA T-STAGE DATA ==========\n")
    print(table(tcga_T_plot$T_stage))
    ############################################################
    # EXTRACT SIGNATURE EFFECTS
    ############################################################
    
    cat("\n========== GEO ==========\n")
    
    geo_t_coef <- coef(summary(t_stage_model))
    
    print(geo_t_coef)
    
    
    cat("\n========== TCGA ==========\n")
    
    tcga_t_coef <- coef(summary(tcga_t_model))
    
    print(tcga_t_coef)
    class(tcga_signature_T_model)
    summary(tcga_signature_T_model)
    tcga_t_coef <- coef(summary(tcga_signature_T_model))
    tcga_t_coef
    ############################################################
    # FIGURE 3 — IDENTIFY THE CORRECT T-STAGE MODELS
    ############################################################
    
    cat("========== tcga_signature_T_model ==========\n")
    class(tcga_signature_T_model)
    print(formula(tcga_signature_T_model))
    summary(tcga_signature_T_model)
    
    cat("\n========== geo_signature_T_model ==========\n")
    class(geo_signature_T_model)
    print(formula(geo_signature_T_model))
    summary(geo_signature_T_model)
    ############################################################
    # FINAL T-STAGE SIGNATURE RESULTS
    ############################################################
    
    extract_t_model <- function(model, cohort) {
      
      co <- summary(model)$coefficients
      
      beta <- co["signature", "Estimate"]
      se   <- co["signature", "Std. Error"]
      p    <- co["signature", "Pr(>|z|)"]
      
      data.frame(
        Cohort = cohort,
        N = nobs(model),
        Beta = beta,
        SE = se,
        OR = exp(beta),
        CI_lower = exp(beta - 1.96 * se),
        CI_upper = exp(beta + 1.96 * se),
        P_value = p
      )
    }
    
    t_stage_results <- rbind(
      extract_t_model(
        geo_signature_T_model,
        "GEO"
      ),
      extract_t_model(
        tcga_signature_T_model,
        "TCGA"
      )
    )
    
    t_stage_results
    ############################################################
    # AUDIT THE T-STAGE FIGURE
    ############################################################
    
    cat("========== FIGURE 3 T DATA ==========\n")
    print(names(fig3_T))
    
    cat("\n========== EXISTING T-STAGE MODELS ==========\n")
    
    cat("\nGEO signature T model:\n")
    print(formula(geo_signature_T_model))
    
    cat("\nTCGA signature T model:\n")
    print(formula(tcga_signature_T_model))
    
    cat("\nExisting T-stage objects:\n")
    objects[grepl("T_model|t_model|T_plot|t_plot",
                  objects,
                  ignore.case = TRUE)]
    ############################################################
    # AUDIT EXACT FIGURE 3 STATISTICS
    ############################################################
    
    cat("========== FIG3 T OBJECT ==========\n")
    str(fig3_T)
    
    cat("\n========== GEO T PLOT ==========\n")
    str(geo_T_plot)
    
    cat("\n========== TCGA T PLOT ==========\n")
    str(tcga_T_plot)
    
    cat("\n========== FIG3 GLEASON ==========\n")
    str(fig3_gleason)
    
    cat("\n========== EXISTING FIGURE OBJECTS ==========\n")
    print(p_fig3)
    ############################################################
    # FIGURE 3 — IDENTIFY EXACT T-STAGE STATISTICAL TEST
    ############################################################
    
    cat("========== p_fig3 ==========\n")
    
    print(class(p_fig3))
    print(p_fig3$labels)
    
    cat("\n========== p_fig3 layers ==========\n")
    print(p_fig3$layers)
    
    cat("\n========== FIG3 T SUMMARY ==========\n")
    
    for (cohort in c("GEO", "TCGA")) {
      
      d <- subset(fig3_T, Cohort == cohort)
      
      cat("\n------------------------------\n")
      cat(cohort, "\n")
      cat("------------------------------\n")
      
      print(table(d$T_stage))
      
      cat("\nWilcoxon test:\n")
      print(wilcox.test(signature ~ T_stage, data = d))
      
      cat("\nLinear model:\n")
      lm_tmp <- lm(signature ~ T_stage, data = d)
      print(summary(lm_tmp))
      
      cat("\nSignature SD:\n")
      print(sd(d$signature, na.rm = TRUE))
    }
    ############################################################
    # FINAL STATISTICAL SUMMARY — MASTER TABLE
    ############################################################
    
    ############################################################
    # 1. GLEASON — ORDINAL LOGISTIC REGRESSION
    ############################################################
    
    gleason_summary <- data.frame(
      Analysis = "Gleason score — ordinal logistic regression",
      
      GEO_OR = exp(
        coef(geo_signature_ordinal_z)["geo_signature_z"]
      ),
      GEO_CI_lower = geo_z_CI[1],
      GEO_CI_upper = geo_z_CI[2],
      GEO_P = geo_z_p,
      
      TCGA_OR = exp(
        coef(tcga_signature_ordinal)["signature"]
      ),
      TCGA_CI_lower = tcga_CI_wald[1],
      TCGA_CI_upper = tcga_CI_wald[2],
      TCGA_P = tcga_p
    )
    
    gleason_summary
    ############################################################
    # 2. T-STAGE — UNADJUSTED WILCOXON
    ############################################################
    
    geo_T_wilcox <- wilcox.test(
      signature ~ T_stage,
      data = geo_T_plot
    )
    
    tcga_T_wilcox <- wilcox.test(
      signature ~ T_stage,
      data = tcga_T_plot
    )
    
    T_unadjusted_summary <- data.frame(
      Analysis = "T stage — unadjusted Wilcoxon",
      
      GEO_T2_N = sum(geo_T_plot$T_stage == "T2"),
      GEO_T3T4_N = sum(geo_T_plot$T_stage == "T3/T4"),
      GEO_P = geo_T_wilcox$p.value,
      
      TCGA_T2_N = sum(tcga_T_plot$T_stage == "T2"),
      TCGA_T3T4_N = sum(tcga_T_plot$T_stage == "T3/T4"),
      TCGA_P = tcga_T_wilcox$p.value
    )
    
    T_unadjusted_summary
    ############################################################
    # 3. T-STAGE — MULTIVARIABLE LOGISTIC REGRESSION
    ############################################################
    
    extract_binary_model <- function(model, cohort) {
      
      co <- summary(model)$coefficients
      
      beta <- co["signature", "Estimate"]
      se   <- co["signature", "Std. Error"]
      p    <- co["signature", "Pr(>|z|)"]
      
      data.frame(
        Cohort = cohort,
        N = nobs(model),
        Beta = beta,
        OR = exp(beta),
        CI_lower = exp(beta - 1.96 * se),
        CI_upper = exp(beta + 1.96 * se),
        P_value = p
      )
    }
    
    T_adjusted_summary <- rbind(
      extract_binary_model(
        geo_signature_T_model,
        "GEO"
      ),
      extract_binary_model(
        tcga_signature_T_model,
        "TCGA"
      )
    )
    
    T_adjusted_summary
    ############################################################
    # 4. GSEA SENSITIVITY ANALYSIS
    ############################################################
    
    GSEA_sensitivity_summary <- data.frame(
      Analysis = c(
        "GO GSEA — signature genes removed",
        "KEGG GSEA — signature genes removed"
      ),
      
      Spearman_rho = c(
        go_rho,
        kegg_rho
      )
    )
    
    GSEA_sensitivity_summary
    ############################################################
    # 5. INDIVIDUAL GENES VS 3-GENE SIGNATURE
    ############################################################
    
    gene_comparison_final <- rbind(
      
      data.frame(
        Cohort = "GEO",
        Model = rownames(geo_model_comparison),
        Beta = geo_model_comparison$Beta,
        P_value = geo_model_comparison$P_value,
        Adjusted_R2 = geo_model_comparison$Adjusted_R2,
        AIC = geo_model_comparison$AIC
      ),
      
      data.frame(
        Cohort = "TCGA",
        Model = rownames(model_comparison),
        Beta = model_comparison$Beta,
        P_value = model_comparison$P_value,
        Adjusted_R2 = model_comparison$Adjusted_R2,
        AIC = model_comparison$AIC
      )
    )
    
    rownames(gene_comparison_final) <- NULL
    
    gene_comparison_final
    ############################################################
    # 6. MASTER RESULTS TABLE
    ############################################################
    
    master_results <- data.frame(
      
      Analysis = c(
        "Gleason — ordinal regression",
        "T stage — unadjusted",
        "T stage — adjusted",
        "GO GSEA sensitivity",
        "KEGG GSEA sensitivity"
      ),
      
      GEO_Result = c(
        sprintf(
          "OR %.3f (95%% CI %.3f–%.3f), P=%.4g",
          gleason_summary$GEO_OR,
          gleason_summary$GEO_CI_lower,
          gleason_summary$GEO_CI_upper,
          gleason_summary$GEO_P
        ),
        
        sprintf(
          "Wilcoxon P=%.4f",
          geo_T_wilcox$p.value
        ),
        
        sprintf(
          "OR %.3f (95%% CI %.3f–%.3f), P=%.4f",
          T_adjusted_summary$OR[
            T_adjusted_summary$Cohort == "GEO"
          ],
          T_adjusted_summary$CI_lower[
            T_adjusted_summary$Cohort == "GEO"
          ],
          T_adjusted_summary$CI_upper[
            T_adjusted_summary$Cohort == "GEO"
          ],
          T_adjusted_summary$P_value[
            T_adjusted_summary$Cohort == "GEO"
          ]
        ),
        
        sprintf(
          "Spearman rho=%.3f",
          go_rho
        ),
        
        sprintf(
          "Spearman rho=%.3f",
          kegg_rho
        )
      ),
      
      TCGA_Result = c(
        sprintf(
          "OR %.3f (95%% CI %.3f–%.3f), P=%.4g",
          gleason_summary$TCGA_OR,
          gleason_summary$TCGA_CI_lower,
          gleason_summary$TCGA_CI_upper,
          gleason_summary$TCGA_P
        ),
        
        sprintf(
          "Wilcoxon P=%.4f",
          tcga_T_wilcox$p.value
        ),
        
        sprintf(
          "OR %.3f (95%% CI %.3f–%.3f), P=%.4f",
          T_adjusted_summary$OR[
            T_adjusted_summary$Cohort == "TCGA"
          ],
          T_adjusted_summary$CI_lower[
            T_adjusted_summary$Cohort == "TCGA"
          ],
          T_adjusted_summary$CI_upper[
            T_adjusted_summary$Cohort == "TCGA"
          ],
          T_adjusted_summary$P_value[
            T_adjusted_summary$Cohort == "TCGA"
          ]
        ),
        
        NA,
        NA
      )
    )
    
    master_results
    ############################################################
    # FINAL FIGURE AUDIT
    ############################################################
    
    cat("========== FIGURE 4 ==========\n")
    print(class(figure4))
    print(figure4)
    
    cat("\n========== FINAL GSEA FIGURE ==========\n")
    print(class(final_GSEA_figure))
    print(final_GSEA_figure)
    
    cat("\n========== SURVIVAL OBJECTS ==========\n")
    
    for (x in c(
      "tcga_signature_cox",
      "tcga_signature_cox_adj",
      "tcga_signature_cox_z",
      "tcga_survival_signature",
      "tcga_survival_complete"
    )) {
      
      cat("\n------------------------------\n")
      cat(x, "\n")
      cat("------------------------------\n")
      
      print(class(get(x)))
      
      if (inherits(get(x), "coxph")) {
        print(formula(get(x)))
        print(summary(get(x)))
      } else if (is.data.frame(get(x))) {
        print(dim(get(x)))
        print(names(get(x)))
        print(head(get(x)))
      }
    }
    ############################################################
    # FIGURE OBJECT INVENTORY
    ############################################################
    
    final_figures <- c(
      "p_fig3",
      "p_gleason_final",
      "figure4",
      "final_GSEA_figure"
    )
    
    data.frame(
      Figure = final_figures,
      Exists = sapply(final_figures, exists),
      Class = sapply(
        final_figures,
        function(x) paste(class(get(x)), collapse = ", ")
      )
    )
    ############################################################
    # FINAL SURVIVAL RESULTS FOR MANUSCRIPT
    ############################################################
    
    cat("========== UNADJUSTED COX ==========\n")
    print(summary(tcga_signature_cox))
    
    cat("\n========== STANDARDIZED SIGNATURE COX ==========\n")
    print(summary(tcga_signature_cox_z))
    
    cat("\n========== ADJUSTED COX ==========\n")
    print(summary(tcga_signature_cox_adj))
    
    cat("\n========== ADJUSTED COX RESULTS TABLE ==========\n")
    print(tcga_survival_adj_results)
    
    cat("\n========== SURVIVAL SAMPLE ==========\n")
    cat("N =", nrow(tcga_survival_complete), "\n")
    cat("Events =", sum(tcga_survival_complete$event), "\n")
    cat("Censored =", sum(tcga_survival_complete$event == 0), "\n")
    
    cat("\n========== PROPORTIONAL HAZARDS TEST ==========\n")
    print(tcga_cox_ph_test)
    ############################################################
    # FIND THE ORIGIN OF THE SIGNATURE WEIGHTS
    ############################################################
    
    grep(
      "0\\.5896508|0\\.5686800|0\\.5735111",
      readLines("TCGA_PRAD_analysis.R"),
      value = TRUE
    )
    ############################################################
    # VERIFY PCA-DERIVED SIGNATURE
    ############################################################
    
    print(tcga_signature_pca)
    
    cat("\n========== PCA LOADINGS ==========\n")
    print(tcga_signature_pca$rotation)
    
    cat("\n========== PCA VARIANCE ==========\n")
    print(summary(tcga_signature_pca))
    
    cat("\n========== STORED WEIGHTS ==========\n")
    print(tcga_weights)
    
    cat("\n========== CORRELATION WITH PC1 ==========\n")
    
    pc1 <- tcga_signature_pca$x[, 1]
    
    print(cor(
      pc1,
      tcga_signature_497,
      method = "pearson"
    ))
    
    print(cor(
      -pc1,
      tcga_signature_497,
      method = "pearson"
    ))
    class(tcga_signature_pca)
    typeof(tcga_signature_pca)
    length(tcga_signature_pca)
    print(tcga_signature_pca)
    cat("\n========== SIGNATURE WEIGHTS ==========\n")
    print(tcga_weights)
    
    cat("\n========== SIGNATURE DEFINITION ==========\n")
    print(signature_definition)
    
    cat("\n========== SIGNATURE FORMULA ==========\n")
    print(signature_formula)
    grep(
      "tcga_signature_pca\\s*<-",
      readLines("TCGA_PRAD_analysis.R"),
      value = TRUE
    )
    sessionInfo()
    cat("========== FINAL GO ==========\n")
    print(
      go_final[, c(
        "Description",
        "NES",
        "p.adjust"
      )]
    )
    
    cat("\n========== FINAL KEGG ==========\n")
    print(
      kegg_final[, c(
        "Description",
        "NES",
        "p.adjust"
      )]
    )    
    cat("\n===== ORIGINAL 11-GENE TCGA RESULTS =====\n")
    
    print(
      cor_results[
        order(cor_results$p_adjusted),
        c("Gene", "Spearman_rho", "p_value", "p_adjusted")
      ],
      row.names = FALSE
    )
    
    cat("\n===== AGE-ADJUSTED RESULTS FOR THE 5-GENE SET =====\n")
    
    print(
      adjusted_results[
        order(adjusted_results$p_adjusted),
      ],
      row.names = FALSE
    )
    
    cat("\n===== TCGA SAMPLE SIZE =====\n")
    cat("N =", nrow(gleason_analysis), "\n")
    
    cat("\n===== GLEASON DISTRIBUTION =====\n")
    print(table(gleason_analysis$gleason_score, useNA = "ifany"))
    cat("\n===== FIVE-GENE CORRELATION MATRIX =====\n")
    
    five_genes <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    print(
      round(
        cor(
          tcga_t_model[, five_genes],
          method = "spearman",
          use = "pairwise.complete.obs"
        ),
        3
      )
    )
    
    
    cat("\n===== FIVE-GENE VIF =====\n")
    
    print(vif_results)
    
    
    cat("\n===== FIVE-GENE AGE-ADJUSTED MODELS =====\n")
    
    print(adjusted_results)
    
    
    cat("\n===== FIVE-GENE MULTIVARIABLE MODEL =====\n")
    
    if (exists("five_gene_model")) {
      print(summary(five_gene_model))
    } else {
      cat("five_gene_model does not exist\n")
    }
    
    
    cat("\n===== OBJECTS RELATED TO GENE SELECTION =====\n")
    
    print(
      objects[
        grepl(
          "five|gene|vif|network|pca|candidate|selection",
          objects,
          ignore.case = TRUE
        )
      ]
    )
    cat("\n===== SEARCH SCRIPT FOR 5-GENE SELECTION =====\n")
    
    script_lines <- readLines("TCGA_PRAD_analysis.R", warn = FALSE)
    
    hits <- grep(
      "MYLK4|VEGFD|MYLK|PRKCB|PDGFC|five_genes|five_gene|vif|VIF|PCA|network",
      script_lines,
      ignore.case = TRUE
    )
    
    for (i in hits) {
      cat("\n--- LINE", i, "---\n")
      start <- max(1, i - 5)
      end   <- min(length(script_lines), i + 8)
      cat(
        paste(
          sprintf("%4d: %s", start:end, script_lines[start:end]),
          collapse = "\n"
        )
      )
      cat("\n")
    }
    cat("\n===== FIVE-GENE VIF VALUES =====\n")
    
    print(
      vif_manual(
        tcga_t_model,
        c("MYLK", "PRKCB", "PDGFC", "VEGFD", "MYLK4")
      )
    )
    
    cat("\n===== FIVE-GENE MULTIVARIABLE MODEL =====\n")
    
    print(summary(
      lm(
        gleason_score ~ MYLK + PRKCB + PDGFC + VEGFD + MYLK4 + age_years,
        data = tcga_t_model
      )
    ))
    cat("\n========== TCGA TABLE 1 ==========\n")
    
    cat("\nAge:\n")
    print(summary(tcga_t_model$age_years))
    
    cat("\nGleason:\n")
    print(table(tcga_t_model$gleason_score, useNA = "ifany"))
    
    cat("\nT stage:\n")
    print(table(tcga_t_model$ajcc_pathologic_t, useNA = "ifany"))
    
    cat("\nN stage:\n")
    print(table(tcga_t_model$ajcc_pathologic_n, useNA = "ifany"))
    
    
    cat("\n========== GEO TABLE 1 ==========\n")
    
    cat("\nAge:\n")
    print(summary(geo_ordinal$age_years))
    
    cat("\nGleason:\n")
    print(table(geo_ordinal$gleason_score, useNA = "ifany"))
    
    cat("\nT stage:\n")
    print(table(geo_ordinal$T_stage, useNA = "ifany"))
    cat("\n========== FIGURE OBJECTS ==========\n")
    
    cat("\n--- p_fig3 ---\n")
    print(p_fig3)
    
    cat("\n--- figure4 ---\n")
    print(figure4)
    
    cat("\n--- final_GSEA_figure ---\n")
    print(final_GSEA_figure)
    
    cat("\n========== LABELS ==========\n")
    
    cat("\np_fig3 labels:\n")
    print(p_fig3$labels)
    
    cat("\nfigure4 labels:\n")
    print(figure4$labels)
    
    cat("\nfinal_GSEA_figure labels:\n")
    print(final_GSEA_figure$labels)
    ############################################################
    # CLEAN SIGNATURE COHORT + TRUE SIGNATURE-ASSOCIATED GSEA
    # Do NOT overwrite the original analysis objects
    ############################################################
    
    library(dplyr)
    library(tidyr)
    library(MASS)
    library(limma)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(survival)
    
    cat("\n====================================================\n")
    cat("1. CLEAN TCGA SIGNATURE-DEVELOPMENT COHORT\n")
    cat("====================================================\n")
    
    # ---------------------------------------------------------
    # Build the clinical annotation independently of T stage.
    # Signature development should not require T-stage variables.
    # ---------------------------------------------------------
    
    tcga_clinical_clean <- clinical_prad %>%
      transmute(
        patient_id = submitter_id,
        age_years = age_at_diagnosis / 365.25,
        gleason_score = gleason_score
      )
    
    # Match to the 497 tumor-patient expression dataset
    tcga_signature_clean <- gleason_analysis %>%
      select(
        patient_id,
        sample_id,
        gleason_score,
        MYLK,
        PRKCB,
        PDGFC
      ) %>%
      left_join(
        tcga_clinical_clean %>%
          select(patient_id, age_years),
        by = "patient_id"
      )
    
    cat("\nInitial rows:\n")
    print(nrow(tcga_signature_clean))
    
    cat("\nMissingness:\n")
    print(
      colSums(
        is.na(
          tcga_signature_clean[
            ,
            c(
              "gleason_score",
              "age_years",
              "MYLK",
              "PRKCB",
              "PDGFC"
            )
          ]
        )
      )
    )
    
    # Complete-case cohort for signature development
    tcga_signature_clean <- tcga_signature_clean %>%
      filter(
        complete.cases(
          gleason_score,
          age_years,
          MYLK,
          PRKCB,
          PDGFC
        )
      )
    
    cat("\nFINAL CLEAN SIGNATURE COHORT N =", nrow(tcga_signature_clean), "\n")
    
    cat("\nGleason distribution:\n")
    print(table(tcga_signature_clean$gleason_score))
    
    cat("\nAge summary:\n")
    print(summary(tcga_signature_clean$age_years))
    
    
    ############################################################
    # 2. REBUILD THE THREE-GENE PCA SIGNATURE IN CLEAN TCGA
    ############################################################
    
    cat("\n====================================================\n")
    cat("2. CLEAN PCA SIGNATURE\n")
    cat("====================================================\n")
    
    clean_network <- tcga_signature_clean %>%
      select(MYLK, PRKCB, PDGFC)
    
    clean_pca <- prcomp(
      clean_network,
      scale. = TRUE
    )
    
    cat("\nPCA variance explained:\n")
    print(summary(clean_pca))
    
    cat("\nPCA rotation:\n")
    print(clean_pca$rotation)
    
    # Fix PCA orientation deterministically:
    # make the MYLK loading positive
    clean_weights <- clean_pca$rotation[, 1]
    
    if (clean_weights["MYLK"] < 0) {
      clean_weights <- -clean_weights
    }
    
    cat("\nFINAL CLEAN PCA WEIGHTS:\n")
    print(clean_weights)
    
    # TCGA means / SDs used for transport
    clean_means <- colMeans(clean_network)
    clean_sds <- apply(clean_network, 2, sd)
    
    cat("\nTCGA CLEAN MEANS:\n")
    print(clean_means)
    
    cat("\nTCGA CLEAN SDs:\n")
    print(clean_sds)
    
    # Construct fixed signature
    tcga_signature_clean$signature_clean <- as.numeric(
      as.matrix(
        scale(
          clean_network,
          center = clean_means,
          scale = clean_sds
        )
      ) %*% clean_weights
    )
    
    cat("\nSignature summary:\n")
    print(summary(tcga_signature_clean$signature_clean))
    
    # Verify that our weighted score reproduces the oriented PC1
    pc1_clean <- as.numeric(clean_pca$x[, 1])
    
    if (cor(pc1_clean, tcga_signature_clean$signature_clean) < 0) {
      pc1_clean <- -pc1_clean
    }
    
    cat(
      "\nCorrelation signature vs oriented PC1 = ",
      cor(
        tcga_signature_clean$signature_clean,
        pc1_clean,
        method = "pearson"
      ),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 3. CLEAN TCGA GLEASON MODEL
    ############################################################
    
    cat("\n====================================================\n")
    cat("3. CLEAN TCGA GLEASON MODEL\n")
    cat("====================================================\n")
    
    tcga_signature_clean$gleason_ordinal <- factor(
      tcga_signature_clean$gleason_score,
      levels = sort(unique(tcga_signature_clean$gleason_score)),
      ordered = TRUE
    )
    
    tcga_clean_ordinal <- MASS::polr(
      gleason_ordinal ~ signature_clean + age_years,
      data = tcga_signature_clean,
      Hess = TRUE
    )
    
    print(summary(tcga_clean_ordinal))
    
    co <- coef(summary(tcga_clean_ordinal))
    beta <- co["signature_clean", "Value"]
    se <- co["signature_clean", "Std. Error"]
    
    tcga_clean_gleason_result <- data.frame(
      N = nrow(tcga_signature_clean),
      OR = exp(beta),
      CI_lower = exp(beta - 1.96 * se),
      CI_upper = exp(beta + 1.96 * se),
      P_value = 2 * pnorm(abs(co["signature_clean", "t value"]),
                          lower.tail = FALSE)
    )
    
    cat("\nCLEAN TCGA GLEASON RESULT:\n")
    print(tcga_clean_gleason_result)
    
    
    ############################################################
    # 4. APPLY THE FIXED TCGA SIGNATURE TO GEO
    ############################################################
    
    cat("\n====================================================\n")
    cat("4. CLEAN GEO VALIDATION SIGNATURE\n")
    cat("====================================================\n")
    
    # geo_gleason already contains the three expression variables
    geo_signature_clean <- as.matrix(
      scale(
        geo_gleason[, c("MYLK", "PRKCB", "PDGFC")],
        center = clean_means,
        scale = clean_sds
      )
    ) %*% clean_weights
    
    geo_signature_clean <- as.numeric(geo_signature_clean)
    
    geo_clean <- geo_gleason
    
    geo_clean$signature_clean <- geo_signature_clean
    
    cat("\nGEO N =", nrow(geo_clean), "\n")
    cat("\nGleason distribution:\n")
    print(table(geo_clean$gleason_score))
    
    # Ordinal validation
    geo_clean$gleason_ordinal <- factor(
      geo_clean$gleason_score,
      levels = sort(unique(geo_clean$gleason_score)),
      ordered = TRUE
    )
    
    geo_clean_ordinal <- MASS::polr(
      gleason_ordinal ~ signature_clean + age_years,
      data = geo_clean,
      Hess = TRUE
    )
    
    print(summary(geo_clean_ordinal))
    
    co_geo <- coef(summary(geo_clean_ordinal))
    beta_geo <- co_geo["signature_clean", "Value"]
    se_geo <- co_geo["signature_clean", "Std. Error"]
    
    geo_clean_gleason_result <- data.frame(
      N = nrow(geo_clean),
      OR = exp(beta_geo),
      CI_lower = exp(beta_geo - 1.96 * se_geo),
      CI_upper = exp(beta_geo + 1.96 * se_geo),
      P_value = 2 * pnorm(abs(co_geo["signature_clean", "t value"]),
                          lower.tail = FALSE)
    )
    
    cat("\nCLEAN GEO GLEASON RESULT:\n")
    print(geo_clean_gleason_result)
    
    
    ############################################################
    # 5. CLEAN T-STAGE ANALYSIS
    # Endpoint-specific complete cases only
    ############################################################
    
    cat("\n====================================================\n")
    cat("5. CLEAN T-STAGE ANALYSIS\n")
    cat("====================================================\n")
    
    # Map T stage from clinical data
    tcga_t_clean <- gleason_analysis %>%
      select(
        patient_id,
        sample_id,
        gleason_score,
        MYLK,
        PRKCB,
        PDGFC
      ) %>%
      left_join(
        tcga_clinical_clean %>%
          left_join(
            clinical_prad %>%
              select(
                submitter_id,
                ajcc_pathologic_t
              ),
            by = c("patient_id" = "submitter_id")
          ) %>%
          select(patient_id, age_years, ajcc_pathologic_t),
        by = "patient_id"
      ) %>%
      left_join(
        tcga_signature_clean %>%
          select(patient_id, signature_clean),
        by = "patient_id"
      ) %>%
      mutate(
        T_advanced_clean = case_when(
          grepl("^T2", ajcc_pathologic_t) ~ 0,
          grepl("^T3|^T4", ajcc_pathologic_t) ~ 1,
          TRUE ~ NA_real_
        )
      )
    
    # Endpoint-specific complete cases
    tcga_T_clean <- tcga_t_clean %>%
      filter(
        complete.cases(
          T_advanced_clean,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat("\nCLEAN TCGA T-STAGE N =", nrow(tcga_T_clean), "\n")
    print(table(tcga_T_clean$T_advanced_clean))
    
    # Unadjusted Wilcoxon
    tcga_T_wilcox_clean <- wilcox.test(
      signature_clean ~ T_advanced_clean,
      data = tcga_T_clean
    )
    
    cat("\nTCGA T-stage Wilcoxon:\n")
    print(tcga_T_wilcox_clean)
    
    # Adjusted logistic model
    tcga_T_model_clean <- glm(
      T_advanced_clean ~ signature_clean + gleason_score + age_years,
      data = tcga_T_clean,
      family = binomial
    )
    
    cat("\nTCGA CLEAN T-STAGE MODEL:\n")
    print(summary(tcga_T_model_clean))
    
    
    # GEO T stage
    geo_T_clean <- geo_ordinal %>%
      select(
        geo_accession,
        T_stage_simple,
        gleason_score,
        age_years
      ) %>%
      left_join(
        geo_clean %>%
          select(geo_accession, signature_clean),
        by = "geo_accession"
      ) %>%
      mutate(
        T_advanced_clean = ifelse(
          T_stage_simple == "T2",
          0,
          ifelse(
            T_stage_simple == "T3/T4",
            1,
            NA
          )
        )
      ) %>%
      filter(
        complete.cases(
          T_advanced_clean,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat("\nCLEAN GEO T-STAGE N =", nrow(geo_T_clean), "\n")
    print(table(geo_T_clean$T_advanced_clean))
    
    geo_T_wilcox_clean <- wilcox.test(
      signature_clean ~ T_advanced_clean,
      data = geo_T_clean
    )
    
    cat("\nGEO T-stage Wilcoxon:\n")
    print(geo_T_wilcox_clean)
    
    geo_T_model_clean <- glm(
      T_advanced_clean ~ signature_clean + gleason_score + age_years,
      data = geo_T_clean,
      family = binomial
    )
    
    cat("\nGEO CLEAN T-STAGE MODEL:\n")
    print(summary(geo_T_model_clean))
    
    
    ############################################################
    # 6. TRUE SIGNATURE-ASSOCIATED GSEA
    ############################################################
    
    cat("\n====================================================\n")
    cat("6. TRUE SIGNATURE-ASSOCIATED GSEA\n")
    cat("====================================================\n")
    
    # We use all eligible TCGA primary-tumor samples in the
    # clean signature cohort, not the old 448-patient T-stage set.
    
    # Make expression matrix correspond to clean patient order
    clean_patient_ids <- tcga_signature_clean$patient_id
    
    clean_sample_ids <- tcga_signature_clean$sample_id
    
    expr_clean <- vsd_tumor[
      ,
      match(clean_sample_ids, colnames(vsd_tumor)),
      drop = FALSE
    ]
    
    cat("\nExpression matrix dimensions:\n")
    print(dim(expr_clean))
    
    # Sanity check
    stopifnot(
      identical(
        colnames(expr_clean),
        tcga_signature_clean$sample_id
      )
    )
    
    # Build covariate matrix
    design_clean <- model.matrix(
      ~ signature_clean + age_years,
      data = tcga_signature_clean
    )
    
    # Fit limma model gene-by-gene
    fit_clean <- limma::lmFit(
      expr_clean,
      design_clean
    )
    
    fit_clean <- limma::eBayes(fit_clean)
    
    signature_coef_col <- which(
      colnames(design_clean) == "signature_clean"
    )
    
    signature_t <- fit_clean$t[, signature_coef_col]
    
    # Annotation: Ensembl -> Entrez
    ensembl_clean <- sub(
      "\\..*$",
      "",
      rownames(expr_clean)
    )
    
    entrez_clean <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_clean,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    # Remove genes with missing Entrez IDs
    keep <- !is.na(entrez_clean) & !is.na(signature_t)
    
    signature_t <- signature_t[keep]
    entrez_clean <- entrez_clean[keep]
    
    # Remove duplicated Entrez IDs
    signature_gsea_df <- data.frame(
      ENTREZID = entrez_clean,
      t_stat = signature_t,
      stringsAsFactors = FALSE
    )
    
    signature_gsea_df <- signature_gsea_df[
      !duplicated(signature_gsea_df$ENTREZID),
    ]
    
    # Rank high-to-low
    gene_list_signature <- signature_gsea_df$t_stat
    names(gene_list_signature) <- signature_gsea_df$ENTREZID
    
    gene_list_signature <- sort(
      gene_list_signature,
      decreasing = TRUE
    )
    
    cat("\nNumber of ranked genes =", length(gene_list_signature), "\n")
    cat("\nTop positive genes:\n")
    print(head(gene_list_signature, 10))
    
    cat("\nTop negative genes:\n")
    print(tail(gene_list_signature, 10))
    
    
    ############################################################
    # 7. PRIMARY TRUE SIGNATURE GSEA
    ############################################################
    
    signature_gsea_go <- clusterProfiler::gseGO(
      geneList = gene_list_signature,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_gsea_kegg <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_go_df <- as.data.frame(signature_gsea_go)
    signature_kegg_df <- as.data.frame(signature_gsea_kegg)
    
    cat("\n===== TRUE SIGNATURE GO GSEA =====\n")
    print(
      signature_go_df[
        order(signature_go_df$p.adjust),
        c(
          "ID",
          "Description",
          "NES",
          "pvalue",
          "p.adjust"
        )
      ][1:min(20, nrow(signature_go_df)), ]
    )
    
    cat("\n===== TRUE SIGNATURE KEGG GSEA =====\n")
    print(
      signature_kegg_df[
        order(signature_kegg_df$p.adjust),
        c(
          "ID",
          "Description",
          "NES",
          "pvalue",
          "p.adjust"
        )
      ][1:min(20, nrow(signature_kegg_df)), ]
    )
    
    
    ############################################################
    # 8. TRUE SIGNATURE GSEA SENSITIVITY:
    #    REMOVE MYLK, PRKCB, PDGFC
    ############################################################
    
    cat("\n====================================================\n")
    cat("8. TRUE SIGNATURE GSEA SENSITIVITY\n")
    cat("====================================================\n")
    
    signature_entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = c("MYLK", "PRKCB", "PDGFC"),
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    signature_entrez <- unname(
      signature_entrez[!is.na(signature_entrez)]
    )
    
    gene_list_signature_no3 <- gene_list_signature[
      !names(gene_list_signature) %in% signature_entrez
    ]
    
    signature_gsea_go_no3 <- clusterProfiler::gseGO(
      geneList = gene_list_signature_no3,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_gsea_kegg_no3 <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature_no3,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    go_no3_df <- as.data.frame(signature_gsea_go_no3)
    kegg_no3_df <- as.data.frame(signature_gsea_kegg_no3)
    
    # Compare common pathways
    go_common_clean <- intersect(
      signature_go_df$ID,
      go_no3_df$ID
    )
    
    kegg_common_clean <- intersect(
      signature_kegg_df$ID,
      kegg_no3_df$ID
    )
    
    go_gsea_clean_rho <- cor(
      signature_go_df$NES[
        match(go_common_clean, signature_go_df$ID)
      ],
      go_no3_df$NES[
        match(go_common_clean, go_no3_df$ID)
      ],
      method = "spearman"
    )
    
    kegg_gsea_clean_rho <- cor(
      signature_kegg_df$NES[
        match(kegg_common_clean, signature_kegg_df$ID)
      ],
      kegg_no3_df$NES[
        match(kegg_common_clean, kegg_no3_df$ID)
      ],
      method = "spearman"
    )
    
    cat("\nGO sensitivity Spearman rho = ",
        go_gsea_clean_rho, "\n", sep = "")
    
    cat(
      "KEGG sensitivity Spearman rho = ",
      kegg_gsea_clean_rho,
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 9. CLEAN TCGA SURVIVAL ANALYSIS
    ############################################################
    
    cat("\n====================================================\n")
    cat("9. CLEAN TCGA SURVIVAL\n")
    cat("====================================================\n")
    
    tcga_survival_clean <- survival_clinical %>%
      select(
        submitter_id,
        survival_time,
        event
      ) %>%
      left_join(
        tcga_signature_clean %>%
          select(
            patient_id,
            signature_clean,
            gleason_score,
            age_years
          ),
        by = c("submitter_id" = "patient_id")
      )
    
    cat("\nSurvival signature cohort N =",
        nrow(tcga_survival_clean), "\n")
    
    cat("Events =",
        sum(tcga_survival_clean$event, na.rm = TRUE), "\n")
    
    # Unadjusted
    tcga_clean_cox <- coxph(
      Surv(survival_time, event) ~ signature_clean,
      data = tcga_survival_clean
    )
    
    cat("\nUNADJUSTED CLEAN COX:\n")
    print(summary(tcga_clean_cox))
    
    # Adjusted complete-case model
    tcga_survival_clean_complete <- tcga_survival_clean %>%
      filter(
        complete.cases(
          survival_time,
          event,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat(
      "\nAdjusted survival cohort N = ",
      nrow(tcga_survival_clean_complete),
      "\n",
      sep = ""
    )
    
    cat(
      "Adjusted survival events = ",
      sum(tcga_survival_clean_complete$event),
      "\n",
      sep = ""
    )
    
    tcga_clean_cox_adj <- coxph(
      Surv(survival_time, event) ~
        scale(signature_clean) +
        gleason_score +
        age_years,
      data = tcga_survival_clean_complete
    )
    
    cat("\nADJUSTED CLEAN COX:\n")
    print(summary(tcga_clean_cox_adj))
    
    cat("\nPROPORTIONAL HAZARDS TEST:\n")
    tcga_clean_ph <- cox.zph(tcga_clean_cox_adj)
    print(tcga_clean_ph)
    
    
    ############################################################
    # 10. MASTER COMPARISON SUMMARY
    ############################################################
    
    cat("\n====================================================\n")
    cat("10. FINAL CLEAN ANALYSIS SUMMARY\n")
    cat("====================================================\n")
    
    cat("\nOLD TCGA SIGNATURE COHORT: 448\n")
    cat("CLEAN TCGA SIGNATURE COHORT: ",
        nrow(tcga_signature_clean), "\n", sep = "")
    
    cat("\nOLD SIGNATURE WEIGHTS:\n")
    print(c(
      MYLK = 0.5896508,
      PRKCB = 0.5686800,
      PDGFC = 0.5735111
    ))
    
    cat("\nNEW CLEAN SIGNATURE WEIGHTS:\n")
    print(clean_weights)
    
    cat("\nOLD TCGA OR:\n")
    print(gleason_summary$TCGA_OR)
    
    cat("\nNEW TCGA OR:\n")
    print(tcga_clean_gleason_result)
    
    cat("\nOLD GEO OR:\n")
    print(gleason_summary$GEO_OR)
    
    cat("\nNEW GEO OR:\n")
    print(geo_clean_gleason_result)
    
    cat("\nOLD TCGA SURVIVAL N:\n")
    print(nrow(tcga_survival_signature))
    
    cat("\nNEW TCGA SURVIVAL N:\n")
    print(nrow(tcga_survival_clean))
    
    cat("\nNEW CLEAN T-STAGE N TCGA:\n")
    print(nrow(tcga_T_clean))
    
    cat("\nNEW CLEAN T-STAGE N GEO:\n")
    print(nrow(geo_T_clean))
    
    cat("\nNEW TRUE GSEA GO RHO:\n")
    print(go_gsea_clean_rho)
    
    cat("\nNEW TRUE GSEA KEGG RHO:\n")
    print(kegg_gsea_clean_rho)
    
    
    ############################################################
    # 11. SAVE CLEAN RESULTS
    ############################################################
    
    write.csv(
      tcga_clean_gleason_result,
      "clean_TCGA_Gleason_result.csv",
      row.names = FALSE
    )
    
    write.csv(
      geo_clean_gleason_result,
      "clean_GEO_Gleason_result.csv",
      row.names = FALSE
    )
    
    write.csv(
      data.frame(
        Gene = names(clean_weights),
        Weight = as.numeric(clean_weights),
        Mean = as.numeric(clean_means[names(clean_weights)]),
        SD = as.numeric(clean_sds[names(clean_weights)])
      ),
      "clean_signature_definition.csv",
      row.names = FALSE
    )
    
    write.csv(
      signature_go_df,
      "true_signature_GSEA_GO.csv",
      row.names = FALSE
    )
    
    write.csv(
      signature_kegg_df,
      "true_signature_GSEA_KEGG.csv",
      row.names = FALSE
    )
    
    write.csv(
      go_no3_df,
      "true_signature_GSEA_GO_no3.csv",
      row.names = FALSE
    )
    
    write.csv(
      kegg_no3_df,
      "true_signature_GSEA_KEGG_no3.csv",
      row.names = FALSE
    )
    
    save(
      tcga_signature_clean,
      clean_pca,
      clean_weights,
      clean_means,
      clean_sds,
      geo_clean,
      tcga_T_clean,
      geo_T_clean,
      signature_gsea_go,
      signature_gsea_kegg,
      signature_gsea_go_no3,
      signature_gsea_kegg_no3,
      tcga_clean_cox,
      tcga_clean_cox_adj,
      tcga_clean_ph,
      file = "clean_signature_GSEA_analysis.RData"
    )
    
    cat("\n\n===== DONE =====\n")
    cat("Clean analysis objects and CSV outputs have been created.\n")
    dplyr::select
    ############################################################
    # CLEAN SIGNATURE COHORT + TRUE SIGNATURE-ASSOCIATED GSEA
    # Do NOT overwrite the original analysis objects
    ############################################################
    
    library(dplyr)
    
    # Explicitly use dplyr functions throughout this analysis
    select <- dplyr::select
    filter <- dplyr::filter
    mutate <- dplyr::mutate
    transmute <- dplyr::transmute
    summarise <- dplyr::summarise
    library(tidyr)
    library(MASS)
    library(limma)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(survival)
    
    cat("\n====================================================\n")
    cat("1. CLEAN TCGA SIGNATURE-DEVELOPMENT COHORT\n")
    cat("====================================================\n")
    
    # ---------------------------------------------------------
    # Build the clinical annotation independently of T stage.
    # Signature development should not require T-stage variables.
    # ---------------------------------------------------------
    
    tcga_clinical_clean <- clinical_prad %>%
      transmute(
        patient_id = submitter_id,
        age_years = age_at_diagnosis / 365.25,
        gleason_score = gleason_score
      )
    
    # Match to the 497 tumor-patient expression dataset
    tcga_signature_clean <- gleason_analysis %>%
      select(
        patient_id,
        sample_id,
        gleason_score,
        MYLK,
        PRKCB,
        PDGFC
      ) %>%
      left_join(
        tcga_clinical_clean %>%
          select(patient_id, age_years),
        by = "patient_id"
      )
    
    cat("\nInitial rows:\n")
    print(nrow(tcga_signature_clean))
    
    cat("\nMissingness:\n")
    print(
      colSums(
        is.na(
          tcga_signature_clean[
            ,
            c(
              "gleason_score",
              "age_years",
              "MYLK",
              "PRKCB",
              "PDGFC"
            )
          ]
        )
      )
    )
    
    # Complete-case cohort for signature development
    tcga_signature_clean <- tcga_signature_clean %>%
      filter(
        complete.cases(
          gleason_score,
          age_years,
          MYLK,
          PRKCB,
          PDGFC
        )
      )
    
    cat("\nFINAL CLEAN SIGNATURE COHORT N =", nrow(tcga_signature_clean), "\n")
    
    cat("\nGleason distribution:\n")
    print(table(tcga_signature_clean$gleason_score))
    
    cat("\nAge summary:\n")
    print(summary(tcga_signature_clean$age_years))
    
    
    ############################################################
    # 2. REBUILD THE THREE-GENE PCA SIGNATURE IN CLEAN TCGA
    ############################################################
    
    cat("\n====================================================\n")
    cat("2. CLEAN PCA SIGNATURE\n")
    cat("====================================================\n")
    
    clean_network <- tcga_signature_clean %>%
      select(MYLK, PRKCB, PDGFC)
    
    clean_pca <- prcomp(
      clean_network,
      scale. = TRUE
    )
    
    cat("\nPCA variance explained:\n")
    print(summary(clean_pca))
    
    cat("\nPCA rotation:\n")
    print(clean_pca$rotation)
    
    # Fix PCA orientation deterministically:
    # make the MYLK loading positive
    clean_weights <- clean_pca$rotation[, 1]
    
    if (clean_weights["MYLK"] < 0) {
      clean_weights <- -clean_weights
    }
    
    cat("\nFINAL CLEAN PCA WEIGHTS:\n")
    print(clean_weights)
    
    # TCGA means / SDs used for transport
    clean_means <- colMeans(clean_network)
    clean_sds <- apply(clean_network, 2, sd)
    
    cat("\nTCGA CLEAN MEANS:\n")
    print(clean_means)
    
    cat("\nTCGA CLEAN SDs:\n")
    print(clean_sds)
    
    # Construct fixed signature
    tcga_signature_clean$signature_clean <- as.numeric(
      as.matrix(
        scale(
          clean_network,
          center = clean_means,
          scale = clean_sds
        )
      ) %*% clean_weights
    )
    
    cat("\nSignature summary:\n")
    print(summary(tcga_signature_clean$signature_clean))
    
    # Verify that our weighted score reproduces the oriented PC1
    pc1_clean <- as.numeric(clean_pca$x[, 1])
    
    if (cor(pc1_clean, tcga_signature_clean$signature_clean) < 0) {
      pc1_clean <- -pc1_clean
    }
    
    cat(
      "\nCorrelation signature vs oriented PC1 = ",
      cor(
        tcga_signature_clean$signature_clean,
        pc1_clean,
        method = "pearson"
      ),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 3. CLEAN TCGA GLEASON MODEL
    ############################################################
    
    cat("\n====================================================\n")
    cat("3. CLEAN TCGA GLEASON MODEL\n")
    cat("====================================================\n")
    
    tcga_signature_clean$gleason_ordinal <- factor(
      tcga_signature_clean$gleason_score,
      levels = sort(unique(tcga_signature_clean$gleason_score)),
      ordered = TRUE
    )
    
    tcga_clean_ordinal <- MASS::polr(
      gleason_ordinal ~ signature_clean + age_years,
      data = tcga_signature_clean,
      Hess = TRUE
    )
    
    print(summary(tcga_clean_ordinal))
    
    co <- coef(summary(tcga_clean_ordinal))
    beta <- co["signature_clean", "Value"]
    se <- co["signature_clean", "Std. Error"]
    
    tcga_clean_gleason_result <- data.frame(
      N = nrow(tcga_signature_clean),
      OR = exp(beta),
      CI_lower = exp(beta - 1.96 * se),
      CI_upper = exp(beta + 1.96 * se),
      P_value = 2 * pnorm(abs(co["signature_clean", "t value"]),
                          lower.tail = FALSE)
    )
    
    cat("\nCLEAN TCGA GLEASON RESULT:\n")
    print(tcga_clean_gleason_result)
    
    
    ############################################################
    # 4. APPLY THE FIXED TCGA SIGNATURE TO GEO
    ############################################################
    
    cat("\n====================================================\n")
    cat("4. CLEAN GEO VALIDATION SIGNATURE\n")
    cat("====================================================\n")
    
    # geo_gleason already contains the three expression variables
    geo_signature_clean <- as.matrix(
      scale(
        geo_gleason[, c("MYLK", "PRKCB", "PDGFC")],
        center = clean_means,
        scale = clean_sds
      )
    ) %*% clean_weights
    
    geo_signature_clean <- as.numeric(geo_signature_clean)
    
    geo_clean <- geo_gleason
    
    geo_clean$signature_clean <- geo_signature_clean
    
    cat("\nGEO N =", nrow(geo_clean), "\n")
    cat("\nGleason distribution:\n")
    print(table(geo_clean$gleason_score))
    
    # Ordinal validation
    geo_clean$gleason_ordinal <- factor(
      geo_clean$gleason_score,
      levels = sort(unique(geo_clean$gleason_score)),
      ordered = TRUE
    )
    
    geo_clean_ordinal <- MASS::polr(
      gleason_ordinal ~ signature_clean + age_years,
      data = geo_clean,
      Hess = TRUE
    )
    
    print(summary(geo_clean_ordinal))
    
    co_geo <- coef(summary(geo_clean_ordinal))
    beta_geo <- co_geo["signature_clean", "Value"]
    se_geo <- co_geo["signature_clean", "Std. Error"]
    
    geo_clean_gleason_result <- data.frame(
      N = nrow(geo_clean),
      OR = exp(beta_geo),
      CI_lower = exp(beta_geo - 1.96 * se_geo),
      CI_upper = exp(beta_geo + 1.96 * se_geo),
      P_value = 2 * pnorm(abs(co_geo["signature_clean", "t value"]),
                          lower.tail = FALSE)
    )
    
    cat("\nCLEAN GEO GLEASON RESULT:\n")
    print(geo_clean_gleason_result)
    
    
    ############################################################
    # 5. CLEAN T-STAGE ANALYSIS
    # Endpoint-specific complete cases only
    ############################################################
    
    cat("\n====================================================\n")
    cat("5. CLEAN T-STAGE ANALYSIS\n")
    cat("====================================================\n")
    
    # Map T stage from clinical data
    tcga_t_clean <- gleason_analysis %>%
      select(
        patient_id,
        sample_id,
        gleason_score,
        MYLK,
        PRKCB,
        PDGFC
      ) %>%
      left_join(
        tcga_clinical_clean %>%
          left_join(
            clinical_prad %>%
              select(
                submitter_id,
                ajcc_pathologic_t
              ),
            by = c("patient_id" = "submitter_id")
          ) %>%
          select(patient_id, age_years, ajcc_pathologic_t),
        by = "patient_id"
      ) %>%
      left_join(
        tcga_signature_clean %>%
          select(patient_id, signature_clean),
        by = "patient_id"
      ) %>%
      mutate(
        T_advanced_clean = case_when(
          grepl("^T2", ajcc_pathologic_t) ~ 0,
          grepl("^T3|^T4", ajcc_pathologic_t) ~ 1,
          TRUE ~ NA_real_
        )
      )
    
    # Endpoint-specific complete cases
    tcga_T_clean <- tcga_t_clean %>%
      filter(
        complete.cases(
          T_advanced_clean,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat("\nCLEAN TCGA T-STAGE N =", nrow(tcga_T_clean), "\n")
    print(table(tcga_T_clean$T_advanced_clean))
    
    # Unadjusted Wilcoxon
    tcga_T_wilcox_clean <- wilcox.test(
      signature_clean ~ T_advanced_clean,
      data = tcga_T_clean
    )
    
    cat("\nTCGA T-stage Wilcoxon:\n")
    print(tcga_T_wilcox_clean)
    
    # Adjusted logistic model
    tcga_T_model_clean <- glm(
      T_advanced_clean ~ signature_clean + gleason_score + age_years,
      data = tcga_T_clean,
      family = binomial
    )
    
    cat("\nTCGA CLEAN T-STAGE MODEL:\n")
    print(summary(tcga_T_model_clean))
    
    
    # GEO T stage
    geo_T_clean <- geo_ordinal %>%
      select(
        geo_accession,
        T_stage_simple,
        gleason_score,
        age_years
      ) %>%
      left_join(
        geo_clean %>%
          select(geo_accession, signature_clean),
        by = "geo_accession"
      ) %>%
      mutate(
        T_advanced_clean = ifelse(
          T_stage_simple == "T2",
          0,
          ifelse(
            T_stage_simple == "T3/T4",
            1,
            NA
          )
        )
      ) %>%
      filter(
        complete.cases(
          T_advanced_clean,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat("\nCLEAN GEO T-STAGE N =", nrow(geo_T_clean), "\n")
    print(table(geo_T_clean$T_advanced_clean))
    
    geo_T_wilcox_clean <- wilcox.test(
      signature_clean ~ T_advanced_clean,
      data = geo_T_clean
    )
    
    cat("\nGEO T-stage Wilcoxon:\n")
    print(geo_T_wilcox_clean)
    
    geo_T_model_clean <- glm(
      T_advanced_clean ~ signature_clean + gleason_score + age_years,
      data = geo_T_clean,
      family = binomial
    )
    
    cat("\nGEO CLEAN T-STAGE MODEL:\n")
    print(summary(geo_T_model_clean))
    
    
    ############################################################
    # 6. TRUE SIGNATURE-ASSOCIATED GSEA
    ############################################################
    
    cat("\n====================================================\n")
    cat("6. TRUE SIGNATURE-ASSOCIATED GSEA\n")
    cat("====================================================\n")
    
    # We use all eligible TCGA primary-tumor samples in the
    # clean signature cohort, not the old 448-patient T-stage set.
    
    # Make expression matrix correspond to clean patient order
    clean_patient_ids <- tcga_signature_clean$patient_id
    
    clean_sample_ids <- tcga_signature_clean$sample_id
    
    expr_clean <- vsd_tumor[
      ,
      match(clean_sample_ids, colnames(vsd_tumor)),
      drop = FALSE
    ]
    
    cat("\nExpression matrix dimensions:\n")
    print(dim(expr_clean))
    
    # Sanity check
    stopifnot(
      identical(
        colnames(expr_clean),
        tcga_signature_clean$sample_id
      )
    )
    
    # Build covariate matrix
    design_clean <- model.matrix(
      ~ signature_clean + age_years,
      data = tcga_signature_clean
    )
    
    # Fit limma model gene-by-gene
    fit_clean <- limma::lmFit(
      expr_clean,
      design_clean
    )
    
    fit_clean <- limma::eBayes(fit_clean)
    
    signature_coef_col <- which(
      colnames(design_clean) == "signature_clean"
    )
    
    signature_t <- fit_clean$t[, signature_coef_col]
    
    # Annotation: Ensembl -> Entrez
    ensembl_clean <- sub(
      "\\..*$",
      "",
      rownames(expr_clean)
    )
    
    entrez_clean <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_clean,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    # Remove genes with missing Entrez IDs
    keep <- !is.na(entrez_clean) & !is.na(signature_t)
    
    signature_t <- signature_t[keep]
    entrez_clean <- entrez_clean[keep]
    
    # Remove duplicated Entrez IDs
    signature_gsea_df <- data.frame(
      ENTREZID = entrez_clean,
      t_stat = signature_t,
      stringsAsFactors = FALSE
    )
    
    signature_gsea_df <- signature_gsea_df[
      !duplicated(signature_gsea_df$ENTREZID),
    ]
    
    # Rank high-to-low
    gene_list_signature <- signature_gsea_df$t_stat
    names(gene_list_signature) <- signature_gsea_df$ENTREZID
    
    gene_list_signature <- sort(
      gene_list_signature,
      decreasing = TRUE
    )
    
    cat("\nNumber of ranked genes =", length(gene_list_signature), "\n")
    cat("\nTop positive genes:\n")
    print(head(gene_list_signature, 10))
    
    cat("\nTop negative genes:\n")
    print(tail(gene_list_signature, 10))
    
    
    ############################################################
    # 7. PRIMARY TRUE SIGNATURE GSEA
    ############################################################
    
    signature_gsea_go <- clusterProfiler::gseGO(
      geneList = gene_list_signature,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_gsea_kegg <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_go_df <- as.data.frame(signature_gsea_go)
    signature_kegg_df <- as.data.frame(signature_gsea_kegg)
    
    cat("\n===== TRUE SIGNATURE GO GSEA =====\n")
    print(
      signature_go_df[
        order(signature_go_df$p.adjust),
        c(
          "ID",
          "Description",
          "NES",
          "pvalue",
          "p.adjust"
        )
      ][1:min(20, nrow(signature_go_df)), ]
    )
    
    cat("\n===== TRUE SIGNATURE KEGG GSEA =====\n")
    print(
      signature_kegg_df[
        order(signature_kegg_df$p.adjust),
        c(
          "ID",
          "Description",
          "NES",
          "pvalue",
          "p.adjust"
        )
      ][1:min(20, nrow(signature_kegg_df)), ]
    )
    
    
    ############################################################
    # 8. TRUE SIGNATURE GSEA SENSITIVITY:
    #    REMOVE MYLK, PRKCB, PDGFC
    ############################################################
    
    cat("\n====================================================\n")
    cat("8. TRUE SIGNATURE GSEA SENSITIVITY\n")
    cat("====================================================\n")
    
    signature_entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = c("MYLK", "PRKCB", "PDGFC"),
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    signature_entrez <- unname(
      signature_entrez[!is.na(signature_entrez)]
    )
    
    gene_list_signature_no3 <- gene_list_signature[
      !names(gene_list_signature) %in% signature_entrez
    ]
    
    signature_gsea_go_no3 <- clusterProfiler::gseGO(
      geneList = gene_list_signature_no3,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_gsea_kegg_no3 <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature_no3,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    go_no3_df <- as.data.frame(signature_gsea_go_no3)
    kegg_no3_df <- as.data.frame(signature_gsea_kegg_no3)
    
    # Compare common pathways
    go_common_clean <- intersect(
      signature_go_df$ID,
      go_no3_df$ID
    )
    
    kegg_common_clean <- intersect(
      signature_kegg_df$ID,
      kegg_no3_df$ID
    )
    
    go_gsea_clean_rho <- cor(
      signature_go_df$NES[
        match(go_common_clean, signature_go_df$ID)
      ],
      go_no3_df$NES[
        match(go_common_clean, go_no3_df$ID)
      ],
      method = "spearman"
    )
    
    kegg_gsea_clean_rho <- cor(
      signature_kegg_df$NES[
        match(kegg_common_clean, signature_kegg_df$ID)
      ],
      kegg_no3_df$NES[
        match(kegg_common_clean, kegg_no3_df$ID)
      ],
      method = "spearman"
    )
    
    cat("\nGO sensitivity Spearman rho = ",
        go_gsea_clean_rho, "\n", sep = "")
    
    cat(
      "KEGG sensitivity Spearman rho = ",
      kegg_gsea_clean_rho,
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 9. CLEAN TCGA SURVIVAL ANALYSIS
    ############################################################
    
    cat("\n====================================================\n")
    cat("9. CLEAN TCGA SURVIVAL\n")
    cat("====================================================\n")
    
    tcga_survival_clean <- survival_clinical %>%
      select(
        submitter_id,
        survival_time,
        event
      ) %>%
      left_join(
        tcga_signature_clean %>%
          select(
            patient_id,
            signature_clean,
            gleason_score,
            age_years
          ),
        by = c("submitter_id" = "patient_id")
      )
    
    cat("\nSurvival signature cohort N =",
        nrow(tcga_survival_clean), "\n")
    
    cat("Events =",
        sum(tcga_survival_clean$event, na.rm = TRUE), "\n")
    
    # Unadjusted
    tcga_clean_cox <- coxph(
      Surv(survival_time, event) ~ signature_clean,
      data = tcga_survival_clean
    )
    
    cat("\nUNADJUSTED CLEAN COX:\n")
    print(summary(tcga_clean_cox))
    
    # Adjusted complete-case model
    tcga_survival_clean_complete <- tcga_survival_clean %>%
      filter(
        complete.cases(
          survival_time,
          event,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat(
      "\nAdjusted survival cohort N = ",
      nrow(tcga_survival_clean_complete),
      "\n",
      sep = ""
    )
    
    cat(
      "Adjusted survival events = ",
      sum(tcga_survival_clean_complete$event),
      "\n",
      sep = ""
    )
    
    tcga_clean_cox_adj <- coxph(
      Surv(survival_time, event) ~
        scale(signature_clean) +
        gleason_score +
        age_years,
      data = tcga_survival_clean_complete
    )
    
    cat("\nADJUSTED CLEAN COX:\n")
    print(summary(tcga_clean_cox_adj))
    
    cat("\nPROPORTIONAL HAZARDS TEST:\n")
    tcga_clean_ph <- cox.zph(tcga_clean_cox_adj)
    print(tcga_clean_ph)
    
    
    ############################################################
    # 10. MASTER COMPARISON SUMMARY
    ############################################################
    
    cat("\n====================================================\n")
    cat("10. FINAL CLEAN ANALYSIS SUMMARY\n")
    cat("====================================================\n")
    
    cat("\nOLD TCGA SIGNATURE COHORT: 448\n")
    cat("CLEAN TCGA SIGNATURE COHORT: ",
        nrow(tcga_signature_clean), "\n", sep = "")
    
    cat("\nOLD SIGNATURE WEIGHTS:\n")
    print(c(
      MYLK = 0.5896508,
      PRKCB = 0.5686800,
      PDGFC = 0.5735111
    ))
    
    cat("\nNEW CLEAN SIGNATURE WEIGHTS:\n")
    print(clean_weights)
    
    cat("\nOLD TCGA OR:\n")
    print(gleason_summary$TCGA_OR)
    
    cat("\nNEW TCGA OR:\n")
    print(tcga_clean_gleason_result)
    
    cat("\nOLD GEO OR:\n")
    print(gleason_summary$GEO_OR)
    
    cat("\nNEW GEO OR:\n")
    print(geo_clean_gleason_result)
    
    cat("\nOLD TCGA SURVIVAL N:\n")
    print(nrow(tcga_survival_signature))
    
    cat("\nNEW TCGA SURVIVAL N:\n")
    print(nrow(tcga_survival_clean))
    
    cat("\nNEW CLEAN T-STAGE N TCGA:\n")
    print(nrow(tcga_T_clean))
    
    cat("\nNEW CLEAN T-STAGE N GEO:\n")
    print(nrow(geo_T_clean))
    
    cat("\nNEW TRUE GSEA GO RHO:\n")
    print(go_gsea_clean_rho)
    
    cat("\nNEW TRUE GSEA KEGG RHO:\n")
    print(kegg_gsea_clean_rho)
    
    
    ############################################################
    # 11. SAVE CLEAN RESULTS
    ############################################################
    
    write.csv(
      tcga_clean_gleason_result,
      "clean_TCGA_Gleason_result.csv",
      row.names = FALSE
    )
    
    write.csv(
      geo_clean_gleason_result,
      "clean_GEO_Gleason_result.csv",
      row.names = FALSE
    )
    
    write.csv(
      data.frame(
        Gene = names(clean_weights),
        Weight = as.numeric(clean_weights),
        Mean = as.numeric(clean_means[names(clean_weights)]),
        SD = as.numeric(clean_sds[names(clean_weights)])
      ),
      "clean_signature_definition.csv",
      row.names = FALSE
    )
    
    write.csv(
      signature_go_df,
      "true_signature_GSEA_GO.csv",
      row.names = FALSE
    )
    
    write.csv(
      signature_kegg_df,
      "true_signature_GSEA_KEGG.csv",
      row.names = FALSE
    )
    
    write.csv(
      go_no3_df,
      "true_signature_GSEA_GO_no3.csv",
      row.names = FALSE
    )
    
    write.csv(
      kegg_no3_df,
      "true_signature_GSEA_KEGG_no3.csv",
      row.names = FALSE
    )
    
    save(
      tcga_signature_clean,
      clean_pca,
      clean_weights,
      clean_means,
      clean_sds,
      geo_clean,
      tcga_T_clean,
      geo_T_clean,
      signature_gsea_go,
      signature_gsea_kegg,
      signature_gsea_go_no3,
      signature_gsea_kegg_no3,
      tcga_clean_cox,
      tcga_clean_cox_adj,
      tcga_clean_ph,
      file = "clean_signature_GSEA_analysis.RData"
    )
    
    cat("\n\n===== DONE =====\n")
    cat("Clean analysis objects and CSV outputs have been created.\n")    
    cat("\n===== GEO CLEAN T-STAGE DEBUG =====\n")
    
    cat("\nT_stage_simple:\n")
    print(table(geo_T_clean$T_stage_simple, useNA = "ifany"))
    
    cat("\nT_advanced_clean:\n")
    print(table(geo_T_clean$T_advanced_clean, useNA = "ifany"))
    
    cat("\nUnique values:\n")
    print(unique(geo_T_clean$T_advanced_clean))
    
    cat("\nStructure:\n")
    print(str(geo_T_clean$T_advanced_clean))
    ############################################################
    # FIX GEO T-STAGE MAPPING
    ############################################################
    
    cat("\n===== SOURCE GEO PATHOLOGICAL STAGE =====\n")
    
    print(
      table(
        geo_ordinal$T_stage,
        useNA = "ifany"
      )
    )
    
    cat("\n===== SOURCE GEO T_STAGE_SIMPLE =====\n")
    
    print(
      table(
        geo_ordinal$T_stage_simple,
        useNA = "ifany"
      )
    )
    
    # Reconstruct the binary T-stage directly from the
    # underlying pathological-stage variable.
    
    geo_T_clean <- geo_ordinal %>%
      dplyr::select(
        geo_accession,
        gleason_score,
        age_years,
        T_stage
      ) %>%
      dplyr::left_join(
        geo_clean %>%
          dplyr::select(
            geo_accession,
            signature_clean
          ),
        by = "geo_accession"
      ) %>%
      dplyr::mutate(
        T_advanced_clean = dplyr::case_when(
          grepl("^T2", T_stage, ignore.case = TRUE) ~ 0,
          grepl("^T3|^T4", T_stage, ignore.case = TRUE) ~ 1,
          TRUE ~ NA_real_
        )
      ) %>%
      dplyr::filter(
        complete.cases(
          T_advanced_clean,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat("\n===== REBUILT GEO T-STAGE =====\n")
    
    print(
      table(
        geo_T_clean$T_advanced_clean,
        useNA = "ifany"
      )
    )
    
    cat("\nExpected:\n")
    cat("T2 =", sum(geo_T_clean$T_advanced_clean == 0), "\n")
    cat("T3/T4 =", sum(geo_T_clean$T_advanced_clean == 1), "\n")
    cat("\n===== CLEAN GEO T-STAGE ANALYSIS =====\n")
    
    geo_T_wilcox_clean <- wilcox.test(
      signature_clean ~ T_advanced_clean,
      data = geo_T_clean
    )
    
    print(geo_T_wilcox_clean)
    
    cat("\n===== CLEAN GEO ADJUSTED T-STAGE MODEL =====\n")
    
    geo_T_model_clean <- glm(
      T_advanced_clean ~ signature_clean + gleason_score + age_years,
      data = geo_T_clean,
      family = binomial
    )
    
    print(summary(geo_T_model_clean))
    print(geo_T_wilcox_clean)
    cat("\n===== TRUE SIGNATURE-ASSOCIATED GSEA =====\n")
    
    cat("\nGO:\n")
    print(
      signature_go_df[
        order(signature_go_df$p.adjust),
        c("Description", "NES", "p.adjust")
      ][1:min(20, nrow(signature_go_df)), ]
    )
    
    cat("\nKEGG:\n")
    print(
      signature_kegg_df[
        order(signature_kegg_df$p.adjust),
        c("Description", "NES", "p.adjust")
      ][1:min(20, nrow(signature_kegg_df)), ]
    )
    
    cat("\n===== GSEA SENSITIVITY =====\n")
    cat("GO Spearman rho =", go_gsea_clean_rho, "\n")
    cat("KEGG Spearman rho =", kegg_gsea_clean_rho, "\n")
    cat("\n===== CLEAN GSEA INPUT CHECK =====\n")
    
    cat("TCGA clean N =", nrow(tcga_signature_clean), "\n")
    
    cat("Expression matrix:\n")
    print(dim(expr_clean))
    
    cat("Signature summary:\n")
    print(summary(tcga_signature_clean$signature_clean))
    cat("\n===== CHECK CLEAN TCGA OBJECTS =====\n")
    
    cat("tcga_signature_clean exists: ",
        exists("tcga_signature_clean"), "\n")
    
    cat("vsd_tumor exists: ",
        exists("vsd_tumor"), "\n")
    
    cat("clean_means exists: ",
        exists("clean_means"), "\n")
    
    cat("clean_sds exists: ",
        exists("clean_sds"), "\n")
    
    cat("\nTCGA clean sample IDs:\n")
    print(head(tcga_signature_clean$sample_id))
    cat("\n===== BUILD CLEAN EXPRESSION MATRIX =====\n")
    
    expr_clean <- vsd_tumor[
      ,
      match(
        tcga_signature_clean$sample_id,
        colnames(vsd_tumor)
      ),
      drop = FALSE
    ]
    
    cat("Expression dimensions:\n")
    print(dim(expr_clean))
    
    cat("Matched samples:\n")
    print(sum(!is.na(
      match(
        tcga_signature_clean$sample_id,
        colnames(vsd_tumor)
      )
    )))
    
    cat("Missing samples:\n")
    print(sum(is.na(
      match(
        tcga_signature_clean$sample_id,
        colnames(vsd_tumor)
      )
    )))
    cat("\n===== TRUE SIGNATURE-ASSOCIATED GSEA =====\n")
    
    # Gene-wise association with the patient-level signature,
    # adjusted for age
    design_clean <- model.matrix(
      ~ signature_clean + age_years,
      data = tcga_signature_clean
    )
    
    fit_clean <- limma::lmFit(
      expr_clean,
      design_clean
    )
    
    fit_clean <- limma::eBayes(fit_clean)
    
    signature_coef_col <- which(
      colnames(design_clean) == "signature_clean"
    )
    
    signature_t <- fit_clean$t[, signature_coef_col]
    
    # Ensembl -> Entrez
    ensembl_clean <- sub(
      "\\..*$",
      "",
      rownames(expr_clean)
    )
    
    entrez_clean <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_clean,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    keep <- !is.na(entrez_clean) & !is.na(signature_t)
    
    signature_gsea_df_input <- data.frame(
      ENTREZID = entrez_clean[keep],
      t_stat = signature_t[keep],
      stringsAsFactors = FALSE
    )
    
    signature_gsea_df_input <- signature_gsea_df_input[
      !duplicated(signature_gsea_df_input$ENTREZID),
    ]
    
    gene_list_signature <- signature_gsea_df_input$t_stat
    names(gene_list_signature) <- signature_gsea_df_input$ENTREZID
    
    gene_list_signature <- sort(
      gene_list_signature,
      decreasing = TRUE
    )
    
    cat("\nNumber of ranked genes:",
        length(gene_list_signature), "\n")
    
    cat("\nTop 10 genes:\n")
    print(head(gene_list_signature, 10))
    
    cat("\nBottom 10 genes:\n")
    print(tail(gene_list_signature, 10))
    
    
    cat("\n===== GO GSEA =====\n")
    
    signature_gsea_go <- clusterProfiler::gseGO(
      geneList = gene_list_signature,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_go_df <- as.data.frame(signature_gsea_go)
    
    print(
      signature_go_df[
        order(signature_go_df$p.adjust),
        c("Description", "NES", "p.adjust")
      ][1:min(20, nrow(signature_go_df)), ]
    )
    
    
    cat("\n===== KEGG GSEA =====\n")
    
    signature_gsea_kegg <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_kegg_df <- as.data.frame(signature_gsea_kegg)
    
    print(
      signature_kegg_df[
        order(signature_kegg_df$p.adjust),
        c("Description", "NES", "p.adjust")
      ][1:min(20, nrow(signature_kegg_df)), ]
    )
    cat("\n===== GSEA WITHOUT MYLK / PRKCB / PDGFC =====\n")
    
    signature_entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = c("MYLK", "PRKCB", "PDGFC"),
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    signature_entrez <- unname(
      signature_entrez[!is.na(signature_entrez)]
    )
    
    gene_list_signature_no3 <- gene_list_signature[
      !names(gene_list_signature) %in% signature_entrez
    ]
    
    signature_gsea_go_no3 <- clusterProfiler::gseGO(
      geneList = gene_list_signature_no3,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    signature_gsea_kegg_no3 <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature_no3,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    go_no3_df <- as.data.frame(signature_gsea_go_no3)
    kegg_no3_df <- as.data.frame(signature_gsea_kegg_no3)
    
    go_common <- intersect(
      signature_go_df$ID,
      go_no3_df$ID
    )
    
    kegg_common <- intersect(
      signature_kegg_df$ID,
      kegg_no3_df$ID
    )
    
    go_gsea_clean_rho <- cor(
      signature_go_df$NES[
        match(go_common, signature_go_df$ID)
      ],
      go_no3_df$NES[
        match(go_common, go_no3_df$ID)
      ],
      method = "spearman"
    )
    
    kegg_gsea_clean_rho <- cor(
      signature_kegg_df$NES[
        match(kegg_common, signature_kegg_df$ID)
      ],
      kegg_no3_df$NES[
        match(kegg_common, kegg_no3_df$ID)
      ],
      method = "spearman"
    )
    
    cat("\nGO sensitivity rho =", go_gsea_clean_rho, "\n")
    cat("KEGG sensitivity rho =", kegg_gsea_clean_rho, "\n")
    
    signature_go_df <- as.data.frame(signature_gsea_go)
    signature_kegg_df <- as.data.frame(signature_gsea_kegg)
    go_no3_df <- as.data.frame(signature_gsea_go_no3)
    kegg_no3_df <- as.data.frame(signature_gsea_kegg_no3)
    
    cat("\n===== FINAL TRUE GO GSEA =====\n")
    print(
      signature_go_df[
        order(signature_go_df$p.adjust),
        c("Description", "NES", "pvalue", "p.adjust")
      ][1:min(20, nrow(signature_go_df)), ]
    )
    
    cat("\n===== FINAL TRUE KEGG GSEA =====\n")
    print(
      signature_kegg_df[
        order(signature_kegg_df$p.adjust),
        c("Description", "NES", "pvalue", "p.adjust")
      ][1:min(20, nrow(signature_kegg_df)), ]
    )
    
    go_common <- intersect(signature_go_df$ID, go_no3_df$ID)
    kegg_common <- intersect(signature_kegg_df$ID, kegg_no3_df$ID)
    
    go_gsea_clean_rho <- cor(
      signature_go_df$NES[match(go_common, signature_go_df$ID)],
      go_no3_df$NES[match(go_common, go_no3_df$ID)],
      method = "spearman"
    )
    
    kegg_gsea_clean_rho <- cor(
      signature_kegg_df$NES[match(kegg_common, signature_kegg_df$ID)],
      kegg_no3_df$NES[match(kegg_common, kegg_no3_df$ID)],
      method = "spearman"
    )
    
    cat("\nGO sensitivity rho =", go_gsea_clean_rho, "\n")
    cat("KEGG sensitivity rho =", kegg_gsea_clean_rho, "\n")
    cat("\n===== CLEAN SIGNATURE DEFINITION =====\n")
    print(data.frame(
      Gene = names(clean_weights),
      Weight = as.numeric(clean_weights),
      Mean = as.numeric(clean_means[names(clean_weights)]),
      SD = as.numeric(clean_sds[names(clean_weights)])
    ))
    
    cat("\n===== CLEAN TCGA GLEASON =====\n")
    print(tcga_clean_gleason_result)
    
    cat("\n===== CLEAN GEO GLEASON =====\n")
    print(geo_clean_gleason_result)
    
    cat("\n===== CLEAN TCGA N =====\n")
    cat(nrow(tcga_signature_clean), "\n")
    
    cat("\n===== CLEAN GEO N =====\n")
    cat(nrow(geo_clean), "\n")
    set.seed(20260821)
    
    signature_gsea_go <- clusterProfiler::gseGO(
      geneList = gene_list_signature,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      eps = 0,
      nPermSimple = 100000,
      verbose = FALSE
    )
    
    set.seed(20260821)
    
    signature_gsea_kegg <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0,
      nPermSimple = 100000,
      verbose = FALSE
    )
    
    set.seed(20260821)
    
    signature_gsea_go_no3 <- clusterProfiler::gseGO(
      geneList = gene_list_signature_no3,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      eps = 0,
      nPermSimple = 100000,
      verbose = FALSE
    )
    
    set.seed(20260821)
    
    signature_gsea_kegg_no3 <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature_no3,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0,
      nPermSimple = 100000,
      verbose = FALSE
    )
    signature_go_final <- as.data.frame(signature_gsea_go)
    signature_kegg_final <- as.data.frame(signature_gsea_kegg)
    
    print(
      signature_go_final[
        order(signature_go_final$p.adjust),
        c("Description", "NES", "pvalue", "p.adjust")
      ][1:10, ]
    )
    
    print(
      signature_kegg_final[
        order(signature_kegg_final$p.adjust),
        c("Description", "NES", "pvalue", "p.adjust")
      ][1:10, ]
    )
    cat("\n===== CLEAN TCGA T-STAGE =====\n")
    print(summary(tcga_T_model_clean))
    print(tcga_T_wilcox_clean)
    
    cat("\n===== CLEAN TCGA SURVIVAL =====\n")
    print(summary(tcga_clean_cox))
    print(summary(tcga_clean_cox_adj))
    
    cat("\n===== CLEAN PH TEST =====\n")
    print(tcga_clean_ph)
    cat("\n===== CLEAN SURVIVAL INPUTS =====\n")
    
    cat("survival_clinical exists: ",
        exists("survival_clinical"), "\n")
    
    cat("tcga_signature_clean exists: ",
        exists("tcga_signature_clean"), "\n")
    
    cat("survival_clinical columns:\n")
    print(
      intersect(
        c("submitter_id", "survival_time", "event"),
        colnames(survival_clinical)
      )
    )
    cat("\n===== BUILD CLEAN SURVIVAL COHORT =====\n")
    
    tcga_survival_clean <- survival_clinical %>%
      dplyr::select(
        submitter_id,
        survival_time,
        event
      ) %>%
      dplyr::left_join(
        tcga_signature_clean %>%
          dplyr::select(
            patient_id,
            signature_clean,
            gleason_score,
            age_years
          ),
        by = c("submitter_id" = "patient_id")
      )
    
    cat("Survival cohort N = ",
        nrow(tcga_survival_clean),
        "\n",
        sep = "")
    
    cat("Events = ",
        sum(tcga_survival_clean$event, na.rm = TRUE),
        "\n",
        sep = "")
    
    cat("\n===== UNADJUSTED CLEAN COX =====\n")
    
    tcga_clean_cox <- survival::coxph(
      survival::Surv(survival_time, event) ~ signature_clean,
      data = tcga_survival_clean
    )
    
    print(summary(tcga_clean_cox))
    
    
    cat("\n===== ADJUSTED CLEAN COX =====\n")
    
    tcga_survival_clean_complete <- tcga_survival_clean %>%
      dplyr::filter(
        complete.cases(
          survival_time,
          event,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cat(
      "Complete survival N = ",
      nrow(tcga_survival_clean_complete),
      "\n",
      sep = ""
    )
    
    cat(
      "Events = ",
      sum(tcga_survival_clean_complete$event),
      "\n",
      sep = ""
    )
    
    tcga_clean_cox_adj <- survival::coxph(
      survival::Surv(survival_time, event) ~
        scale(signature_clean) +
        gleason_score +
        age_years,
      data = tcga_survival_clean_complete
    )
    
    print(summary(tcga_clean_cox_adj))
    
    
    cat("\n===== PROPORTIONAL HAZARDS =====\n")
    
    tcga_clean_ph <- survival::cox.zph(
      tcga_clean_cox_adj
    )
    
    print(tcga_clean_ph)
    cat("\n===== FINAL CLEAN TCGA TABLE 1 =====\n")
    
    cat("\nN:\n")
    cat("Clean signature cohort N =", nrow(tcga_signature_clean), "\n")
    
    cat("\nAGE:\n")
    print(summary(tcga_signature_clean$age_years))
    
    cat("\nGLEASON:\n")
    print(table(tcga_signature_clean$gleason_score, useNA = "ifany"))
    
    cat("\nPATHOLOGICAL T STAGE IN CLEAN COHORT:\n")
    
    tcga_table1_clean <- tcga_signature_clean %>%
      dplyr::select(patient_id) %>%
      dplyr::left_join(
        clinical_prad %>%
          dplyr::select(
            submitter_id,
            ajcc_pathologic_t,
            ajcc_pathologic_n
          ),
        by = c("patient_id" = "submitter_id")
      )
    
    print(
      table(
        tcga_table1_clean$ajcc_pathologic_t,
        useNA = "ifany"
      )
    )
    
    cat("\nN STAGE:\n")
    
    print(
      table(
        tcga_table1_clean$ajcc_pathologic_n,
        useNA = "ifany"
      )
    )
    
    cat("\n===== FINAL CLEAN GEO TABLE 1 =====\n")
    
    cat("\nN:\n")
    cat("GEO N =", nrow(geo_clean), "\n")
    
    cat("\nAGE:\n")
    print(summary(geo_clean$age_years))
    
    cat("\nGLEASON:\n")
    print(table(geo_clean$gleason_score, useNA = "ifany"))
    
    cat("\nPATHOLOGICAL T STAGE:\n")
    print(
      table(
        geo_ordinal$T_stage,
        useNA = "ifany"
      )
    )
    ############################################################
    # OBJECTIVE 3-GENE SIGNATURE SELECTION
    # 10 possible combinations of the 5 candidate genes
    # Repeated stratified 5-fold cross-validation
    #
    # IMPORTANT:
    # PCA scaling + PCA weights are fitted WITHIN each training fold.
    # GEO is NOT used for selection.
    ############################################################
    
    library(MASS)
    
    set.seed(20260822)
    
    candidate_genes <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    all_combos <- combn(
      candidate_genes,
      3,
      simplify = FALSE
    )
    
    cat("\n===== CANDIDATE COMBINATIONS =====\n")
    
    for (i in seq_along(all_combos)) {
      cat(
        i, ": ",
        paste(all_combos[[i]], collapse = " + "),
        "\n",
        sep = ""
      )
    }
    
    
    ############################################################
    # PREPARE DATA
    ############################################################
    
    cv_data <- tcga_signature_clean %>%
      dplyr::select(
        patient_id,
        gleason_score,
        age_years,
        dplyr::all_of(candidate_genes)
      ) %>%
      dplyr::filter(
        complete.cases(.)
      )
    
    cv_data$gleason_ordinal <- factor(
      cv_data$gleason_score,
      levels = c(6, 7, 8, 9, 10),
      ordered = TRUE
    )
    
    cat("\n===== CV DATA =====\n")
    cat("N =", nrow(cv_data), "\n")
    print(table(cv_data$gleason_ordinal))
    
    
    ############################################################
    # STRATIFIED REPEATED FOLD GENERATOR
    ############################################################
    
    make_stratified_folds <- function(y, K = 5, seed = NULL) {
      
      if (!is.null(seed)) {
        set.seed(seed)
      }
      
      folds <- integer(length(y))
      
      for (lev in levels(y)) {
        
        idx <- which(y == lev)
        
        idx <- sample(idx)
        
        fold_sequence <- rep(
          seq_len(K),
          length.out = length(idx)
        )
        
        folds[idx] <- fold_sequence
      }
      
      folds
    }
    
    
    ############################################################
    # FUNCTION TO FIT ONE PCA SIGNATURE
    ############################################################
    
    fit_fold_signature <- function(train_df, test_df, genes) {
      
      train_x <- as.matrix(
        train_df[, genes, drop = FALSE]
      )
      
      test_x <- as.matrix(
        test_df[, genes, drop = FALSE]
      )
      
      # Standardize using TRAINING data only
      train_scaled <- scale(
        train_x,
        center = TRUE,
        scale = TRUE
      )
      
      train_center <- attr(
        train_scaled,
        "scaled:center"
      )
      
      train_scale <- attr(
        train_scaled,
        "scaled:scale"
      )
      
      # PCA using training data only
      pca <- prcomp(
        train_scaled,
        center = FALSE,
        scale. = FALSE
      )
      
      weights <- pca$rotation[, 1]
      
      # Deterministic orientation:
      # largest absolute loading is positive
      anchor <- which.max(abs(weights))
      
      if (weights[anchor] < 0) {
        weights <- -weights
      }
      
      # Training signature
      train_signature <- as.numeric(
        train_scaled %*% weights
      )
      
      # Test data standardized using TRAINING parameters
      test_scaled <- scale(
        test_x,
        center = train_center,
        scale = train_scale
      )
      
      test_signature <- as.numeric(
        test_scaled %*% weights
      )
      
      list(
        train_signature = train_signature,
        test_signature = test_signature,
        weights = weights,
        center = train_center,
        scale = train_scale
      )
    }
    
    
    ############################################################
    # PERFORMANCE METRICS
    ############################################################
    
    multiclass_logloss <- function(prob, truth) {
      
      truth_chr <- as.character(truth)
      
      true_prob <- prob[
        cbind(
          seq_along(truth_chr),
          match(truth_chr, colnames(prob))
        )
      ]
      
      true_prob <- pmax(
        true_prob,
        1e-15
      )
      
      -mean(log(true_prob))
    }
    
    
    ordinal_mae <- function(prob, truth) {
      
      class_values <- as.numeric(
        colnames(prob)
      )
      
      expected_value <- as.numeric(
        prob %*% class_values
      )
      
      mean(
        abs(
          expected_value -
            as.numeric(as.character(truth))
        )
      )
    }
    
    
    ############################################################
    # REPEATED CROSS-VALIDATION
    ############################################################
    
    K <- 5
    N_REPEATS <- 5
    
    cv_results <- list()
    
    counter <- 1
    
    for (combo_id in seq_along(all_combos)) {
      
      genes <- all_combos[[combo_id]]
      
      combo_name <- paste(
        genes,
        collapse = " + "
      )
      
      cat(
        "\nRunning combination:",
        combo_id,
        "/",
        length(all_combos),
        combo_name,
        "\n"
      )
      
      for (rep_id in seq_len(N_REPEATS)) {
        
        folds <- make_stratified_folds(
          cv_data$gleason_ordinal,
          K = K,
          seed = 20260822 + rep_id
        )
        
        for (fold_id in seq_len(K)) {
          
          test_idx <- which(
            folds == fold_id
          )
          
          train_idx <- which(
            folds != fold_id
          )
          
          train_df <- cv_data[train_idx, , drop = FALSE]
          test_df <- cv_data[test_idx, , drop = FALSE]
          
          # Fit PCA-derived signature inside training fold
          sig_fit <- fit_fold_signature(
            train_df,
            test_df,
            genes
          )
          
          train_df$cv_signature <- sig_fit$train_signature
          test_df$cv_signature <- sig_fit$test_signature
          
          # Ordinal regression
          model <- MASS::polr(
            gleason_ordinal ~ cv_signature + age_years,
            data = train_df,
            Hess = FALSE
          )
          
          # Predict probabilities on held-out patients
          prob <- predict(
            model,
            newdata = test_df,
            type = "probs"
          )
          
          logloss <- multiclass_logloss(
            prob,
            test_df$gleason_ordinal
          )
          
          mae <- ordinal_mae(
            prob,
            test_df$gleason_ordinal
          )
          
          cv_results[[counter]] <- data.frame(
            Combination = combo_name,
            Combo_ID = combo_id,
            Repeat = rep_id,
            Fold = fold_id,
            LogLoss = logloss,
            MAE = mae,
            stringsAsFactors = FALSE
          )
          
          counter <- counter + 1
        }
      }
    }
    
    
    ############################################################
    # COMBINE RESULTS
    ############################################################
    
    cv_results <- dplyr::bind_rows(
      cv_results
    )
    
    cat("\n===== RAW CV RESULTS =====\n")
    print(head(cv_results))
    
    
    ############################################################
    # RANK ALL COMBINATIONS
    ############################################################
    
    cv_summary <- cv_results %>%
      dplyr::group_by(
        Combo_ID,
        Combination
      ) %>%
      dplyr::summarise(
        Mean_LogLoss = mean(LogLoss),
        SD_LogLoss = sd(LogLoss),
        Mean_MAE = mean(MAE),
        SD_MAE = sd(MAE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        Mean_LogLoss
      ) %>%
      dplyr::mutate(
        LogLoss_Rank = row_number()
      )
    
    cat("\n====================================================\n")
    cat("FINAL 3-GENE CROSS-VALIDATION RANKING\n")
    cat("====================================================\n")
    
    print(
      cv_summary,
      row.names = FALSE
    )
    
    
    ############################################################
    # HOW OFTEN EACH COMBINATION WON EACH FOLD
    ############################################################
    
    fold_winners <- cv_results %>%
      dplyr::group_by(
        Repeat,
        Fold
      ) %>%
      dplyr::slice_min(
        LogLoss,
        n = 1,
        with_ties = FALSE
      ) %>%
      dplyr::ungroup()
    
    winner_stability <- fold_winners %>%
      dplyr::count(
        Combination,
        name = "Wins"
      ) %>%
      dplyr::arrange(
        dplyr::desc(Wins)
      )
    
    cat("\n====================================================\n")
    cat("CROSS-VALIDATION WIN FREQUENCY\n")
    cat("====================================================\n")
    
    print(
      winner_stability,
      row.names = FALSE
    )
    
    
    ############################################################
    # CURRENT MYLK + PRKCB + PDGFC COMBINATION
    ############################################################
    
    current_combo_name <- paste(
      c("MYLK", "PRKCB", "PDGFC"),
      collapse = " + "
    )
    
    cat("\n====================================================\n")
    cat("CURRENT SIGNATURE CHECK\n")
    cat("====================================================\n")
    
    print(
      cv_summary[
        cv_summary$Combination == current_combo_name,
        ,
        drop = FALSE
      ]
    )
    
    
    ############################################################
    # BEST COMBINATION
    ############################################################
    
    best_combo <- cv_summary$Combination[1]
    
    cat("\n====================================================\n")
    cat("BEST THREE-GENE COMBINATION\n")
    cat("====================================================\n")
    
    cat(best_combo, "\n")
    
    cat("\nFull ranking:\n")
    print(cv_summary)
    
    
    ############################################################
    # SAVE RESULTS
    ############################################################
    
    write.csv(
      cv_results,
      "three_gene_CV_raw_results.csv",
      row.names = FALSE
    )
    
    write.csv(
      cv_summary,
      "three_gene_CV_summary.csv",
      row.names = FALSE
    )
    
    write.csv(
      winner_stability,
      "three_gene_CV_winner_stability.csv",
      row.names = FALSE
    )
    
    cat("\n===== DONE =====\n")
    ############################################################
    # BUILD THE FIVE-GENE DISCOVERY DATASET
    # Used ONLY for comparing all 10 possible 3-gene signatures
    ############################################################
    
    candidate_genes <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    cat("\n===== SOURCE COLUMN CHECK =====\n")
    
    print(
      intersect(
        c(
          "patient_id",
          "sample_id",
          "gleason_score",
          "age_years",
          candidate_genes
        ),
        colnames(gleason_analysis)
      )
    )
    
    # Start from the original TCGA tumor analysis object,
    # which contains the five candidate genes.
    cv_source <- gleason_analysis %>%
      dplyr::select(
        patient_id,
        sample_id,
        gleason_score,
        dplyr::all_of(candidate_genes)
      ) %>%
      dplyr::left_join(
        clinical_prad %>%
          dplyr::select(
            submitter_id,
            age_at_diagnosis
          ) %>%
          dplyr::mutate(
            age_years = age_at_diagnosis / 365.25
          ) %>%
          dplyr::select(
            submitter_id,
            age_years
          ),
        by = c("patient_id" = "submitter_id")
      )
    
    cat("\n===== FIVE-GENE SOURCE =====\n")
    
    cat("Initial N =", nrow(cv_source), "\n")
    
    cat("\nMissingness:\n")
    
    print(
      colSums(
        is.na(
          cv_source[
            ,
            c(
              "gleason_score",
              "age_years",
              candidate_genes
            )
          ]
        )
      )
    )
    
    # Complete cases required for model selection
    cv_data <- cv_source %>%
      dplyr::filter(
        complete.cases(
          gleason_score,
          age_years,
          dplyr::across(
            dplyr::all_of(candidate_genes)
          )
        )
      )
    
    cat("\n===== FINAL FIVE-GENE CV DATASET =====\n")
    
    cat("N =", nrow(cv_data), "\n")
    
    cat("\nGleason distribution:\n")
    print(
      table(
        cv_data$gleason_score,
        useNA = "ifany"
      )
    )
    
    cat("\nAge summary:\n")
    print(
      summary(cv_data$age_years)
    )
    
    cat("\nGene columns:\n")
    print(
      colnames(
        cv_data[
          ,
          candidate_genes,
          drop = FALSE
        ]
      )
    )
    ############################################################
    # OBJECTIVE 3-GENE SIGNATURE SELECTION
    # ALL 10 COMBINATIONS OF THE 5 CANDIDATE GENES
    #
    # Discovery cohort:
    #   455 TCGA patients
    #
    # Method:
    #   Repeated stratified 5-fold cross-validation
    #   PCA fitted INSIDE each training fold
    #   Ordinal logistic regression
    #   Age-adjusted
    #
    # GEO is NOT used for selection.
    ############################################################
    
    library(MASS)
    library(dplyr)
    
    set.seed(20260822)
    
    ############################################################
    # 1. DEFINE ALL 10 THREE-GENE COMBINATIONS
    ############################################################
    
    candidate_genes <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    all_combos <- combn(
      candidate_genes,
      3,
      simplify = FALSE
    )
    
    cat("\n===== CANDIDATE COMBINATIONS =====\n")
    
    for (i in seq_along(all_combos)) {
      
      cat(
        i,
        ": ",
        paste(all_combos[[i]], collapse = " + "),
        "\n",
        sep = ""
      )
      
    }
    
    
    ############################################################
    # 2. CHECK CV DATA
    ############################################################
    
    cat("\n===== CV DATA CHECK =====\n")
    
    cat("N =", nrow(cv_data), "\n")
    
    print(
      table(
        cv_data$gleason_score,
        useNA = "ifany"
      )
    )
    
    print(
      colSums(
        is.na(
          cv_data[
            ,
            c(
              "gleason_score",
              "age_years",
              candidate_genes
            )
          ]
        )
      )
    )
    
    # Ordered Gleason outcome
    cv_data$gleason_ordinal <- factor(
      cv_data$gleason_score,
      levels = c(6, 7, 8, 9, 10),
      ordered = TRUE
    )
    
    
    ############################################################
    # 3. STRATIFIED FOLD FUNCTION
    ############################################################
    
    make_stratified_folds <- function(
    y,
    K = 5,
    seed = NULL
    ) {
      
      if (!is.null(seed)) {
        set.seed(seed)
      }
      
      folds <- integer(length(y))
      
      for (lev in levels(y)) {
        
        idx <- which(y == lev)
        
        idx <- sample(idx)
        
        fold_assignments <- rep(
          seq_len(K),
          length.out = length(idx)
        )
        
        folds[idx] <- fold_assignments
      }
      
      folds
    }
    
    
    ############################################################
    # 4. FIT PCA SIGNATURE WITHIN TRAINING FOLD
    ############################################################
    
    fit_fold_signature <- function(
    train_df,
    test_df,
    genes
    ) {
      
      train_x <- as.matrix(
        train_df[, genes, drop = FALSE]
      )
      
      test_x <- as.matrix(
        test_df[, genes, drop = FALSE]
      )
      
      # Standardize using TRAINING data only
      train_scaled <- scale(
        train_x,
        center = TRUE,
        scale = TRUE
      )
      
      train_center <- attr(
        train_scaled,
        "scaled:center"
      )
      
      train_scale <- attr(
        train_scaled,
        "scaled:scale"
      )
      
      # PCA using TRAINING data only
      pca <- prcomp(
        train_scaled,
        center = FALSE,
        scale. = FALSE
      )
      
      weights <- pca$rotation[, 1]
      
      # Fix arbitrary PCA sign
      anchor <- which.max(abs(weights))
      
      if (weights[anchor] < 0) {
        weights <- -weights
      }
      
      # Training signature
      train_signature <- as.numeric(
        train_scaled %*% weights
      )
      
      # Test standardized using TRAINING parameters
      test_scaled <- scale(
        test_x,
        center = train_center,
        scale = train_scale
      )
      
      test_signature <- as.numeric(
        test_scaled %*% weights
      )
      
      list(
        train_signature = train_signature,
        test_signature = test_signature,
        weights = weights
      )
    }
    
    
    ############################################################
    # 5. PERFORMANCE METRICS
    ############################################################
    
    multiclass_logloss <- function(
    prob,
    truth
    ) {
      
      truth_chr <- as.character(truth)
      
      truth_columns <- match(
        truth_chr,
        colnames(prob)
      )
      
      true_probability <- prob[
        cbind(
          seq_along(truth_chr),
          truth_columns
        )
      ]
      
      true_probability <- pmax(
        true_probability,
        1e-15
      )
      
      -mean(
        log(true_probability)
      )
    }
    
    
    ordinal_mae <- function(
    prob,
    truth
    ) {
      
      class_values <- as.numeric(
        colnames(prob)
      )
      
      expected_value <- as.numeric(
        prob %*% class_values
      )
      
      mean(
        abs(
          expected_value -
            as.numeric(
              as.character(truth)
            )
        )
      )
    }
    
    
    ############################################################
    # 6. CROSS-VALIDATION SETTINGS
    ############################################################
    
    K <- 5
    
    N_REPEATS <- 5
    
    cv_results <- list()
    
    counter <- 1
    
    
    ############################################################
    # 7. RUN ALL 10 COMBINATIONS
    ############################################################
    
    for (combo_id in seq_along(all_combos)) {
      
      genes <- all_combos[[combo_id]]
      
      combo_name <- paste(
        genes,
        collapse = " + "
      )
      
      cat(
        "\n============================================\n"
      )
      
      cat(
        "Combination ",
        combo_id,
        "/",
        length(all_combos),
        ": ",
        combo_name,
        "\n",
        sep = ""
      )
      
      cat(
        "============================================\n"
      )
      
      for (rep_id in seq_len(N_REPEATS)) {
        
        folds <- make_stratified_folds(
          cv_data$gleason_ordinal,
          K = K,
          seed = 20260822 + rep_id
        )
        
        for (fold_id in seq_len(K)) {
          
          test_idx <- which(
            folds == fold_id
          )
          
          train_idx <- which(
            folds != fold_id
          )
          
          train_df <- cv_data[
            train_idx,
            ,
            drop = FALSE
          ]
          
          test_df <- cv_data[
            test_idx,
            ,
            drop = FALSE
          ]
          
          ######################################################
          # FIT PCA USING TRAINING SET ONLY
          ######################################################
          
          sig_fit <- fit_fold_signature(
            train_df = train_df,
            test_df = test_df,
            genes = genes
          )
          
          train_df$cv_signature <-
            sig_fit$train_signature
          
          test_df$cv_signature <-
            sig_fit$test_signature
          
          
          ######################################################
          # ORDINAL LOGISTIC REGRESSION
          ######################################################
          
          model <- MASS::polr(
            gleason_ordinal ~
              cv_signature +
              age_years,
            data = train_df,
            Hess = FALSE
          )
          
          
          ######################################################
          # PREDICT HELD-OUT PATIENTS
          ######################################################
          
          prob <- predict(
            model,
            newdata = test_df,
            type = "probs"
          )
          
          
          ######################################################
          # PERFORMANCE
          ######################################################
          
          logloss <- multiclass_logloss(
            prob = prob,
            truth = test_df$gleason_ordinal
          )
          
          mae <- ordinal_mae(
            prob = prob,
            truth = test_df$gleason_ordinal
          )
          
          
          ######################################################
          # STORE RESULT
          ######################################################
          
          cv_results[[counter]] <- data.frame(
            Combination = combo_name,
            Combo_ID = combo_id,
            Repeat = rep_id,
            Fold = fold_id,
            LogLoss = logloss,
            MAE = mae,
            stringsAsFactors = FALSE
          )
          
          counter <- counter + 1
          
        }
      }
    }
    
    
    ############################################################
    # 8. COMBINE RESULTS
    ############################################################
    
    cv_results <- dplyr::bind_rows(
      cv_results
    )
    
    cat(
      "\n===== RAW CV RESULTS =====\n"
    )
    
    print(
      head(cv_results, 20)
    )
    
    
    ############################################################
    # 9. SUMMARIZE ALL 10 COMBINATIONS
    ############################################################
    
    cv_summary <- cv_results %>%
      dplyr::group_by(
        Combo_ID,
        Combination
      ) %>%
      dplyr::summarise(
        Mean_LogLoss =
          mean(LogLoss),
        
        SD_LogLoss =
          sd(LogLoss),
        
        Mean_MAE =
          mean(MAE),
        
        SD_MAE =
          sd(MAE),
        
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        Mean_LogLoss
      ) %>%
      dplyr::mutate(
        LogLoss_Rank =
          dplyr::row_number()
      )
    
    
    ############################################################
    # 10. FINAL RANKING
    ############################################################
    
    cat(
      "\n====================================================\n"
    )
    
    cat(
      "FINAL 3-GENE CROSS-VALIDATION RANKING\n"
    )
    
    cat(
      "====================================================\n"
    )
    
    print(
      cv_summary,
      row.names = FALSE
    )
    
    
    ############################################################
    # 11. WIN FREQUENCY
    ############################################################
    
    fold_winners <- cv_results %>%
      dplyr::group_by(
        Repeat,
        Fold
      ) %>%
      dplyr::slice_min(
        order_by = LogLoss,
        n = 1,
        with_ties = FALSE
      ) %>%
      dplyr::ungroup()
    
    
    winner_stability <- fold_winners %>%
      dplyr::count(
        Combination,
        name = "Wins"
      ) %>%
      dplyr::arrange(
        dplyr::desc(Wins)
      )
    
    
    cat(
      "\n====================================================\n"
    )
    
    cat(
      "CROSS-VALIDATION WIN FREQUENCY\n"
    )
    
    cat(
      "====================================================\n"
    )
    
    print(
      winner_stability,
      row.names = FALSE
    )
    
    
    ############################################################
    # 12. CHECK THE CURRENT SIGNATURE
    ############################################################
    
    current_combo_name <- paste(
      c(
        "MYLK",
        "PRKCB",
        "PDGFC"
      ),
      collapse = " + "
    )
    
    cat(
      "\n====================================================\n"
    )
    
    cat(
      "CURRENT MYLK + PRKCB + PDGFC SIGNATURE\n"
    )
    
    cat(
      "====================================================\n"
    )
    
    print(
      cv_summary[
        cv_summary$Combination ==
          current_combo_name,
        ,
        drop = FALSE
      ],
      row.names = FALSE
    )
    
    
    ############################################################
    # 13. BEST COMBINATION
    ############################################################
    
    best_combo <- cv_summary$Combination[1]
    
    cat(
      "\n====================================================\n"
    )
    
    cat(
      "BEST THREE-GENE COMBINATION\n"
    )
    
    cat(
      "====================================================\n"
    )
    
    cat(
      best_combo,
      "\n"
    )
    
    
    ############################################################
    # 14. SAVE RESULTS
    ############################################################
    
    write.csv(
      cv_results,
      "three_gene_CV_raw_results.csv",
      row.names = FALSE
    )
    
    write.csv(
      cv_summary,
      "three_gene_CV_summary.csv",
      row.names = FALSE
    )
    
    write.csv(
      winner_stability,
      "three_gene_CV_winner_stability.csv",
      row.names = FALSE
    )
    
    
    ############################################################
    # DONE
    ############################################################
    
    cat(
      "\n===== DONE =====\n"
    )
    ############################################################
    # COMPARE THE CV DIFFERENCES
    ############################################################
    
    cat("\n===== FULL-PRECISION CV SUMMARY =====\n")
    
    print(
      cv_summary %>%
        dplyr::arrange(Mean_LogLoss) %>%
        dplyr::mutate(
          Delta_vs_Best =
            Mean_LogLoss - min(Mean_LogLoss)
        ),
      digits = 8,
      row.names = FALSE
    )
    
    
    ############################################################
    # CURRENT VS BEST
    ############################################################
    
    best_name <- cv_summary$Combination[1]
    
    current_name <- "MYLK + PRKCB + PDGFC"
    
    best_scores <- cv_results %>%
      dplyr::filter(
        Combination == best_name
      ) %>%
      dplyr::arrange(Repeat, Fold)
    
    current_scores <- cv_results %>%
      dplyr::filter(
        Combination == current_name
      ) %>%
      dplyr::arrange(Repeat, Fold)
    
    paired_comparison <- data.frame(
      Repeat = best_scores$Repeat,
      Fold = best_scores$Fold,
      Best_LogLoss = best_scores$LogLoss,
      Current_LogLoss = current_scores$LogLoss
    )
    
    paired_comparison$Difference <-
      paired_comparison$Current_LogLoss -
      paired_comparison$Best_LogLoss
    
    cat("\n===== BEST VS CURRENT SIGNATURE =====\n")
    
    print(
      summary(paired_comparison$Difference)
    )
    
    cat("\nMean difference:\n")
    print(
      mean(paired_comparison$Difference)
    )
    
    cat("\nSD difference:\n")
    print(
      sd(paired_comparison$Difference)
    )
    
    cat("\nPaired test:\n")
    print(
      t.test(
        paired_comparison$Difference,
        mu = 0
      )
    )
    
    cat("\nWilcoxon paired test:\n")
    print(
      wilcox.test(
        paired_comparison$Best_LogLoss,
        paired_comparison$Current_LogLoss,
        paired = TRUE,
        exact = FALSE
      )
    )
    
    
    ############################################################
    # AGE-ONLY BASELINE
    ############################################################
    
    cat("\n===== AGE-ONLY BASELINE =====\n")
    
    age_only_results <- list()
    counter <- 1
    
    for (rep_id in seq_len(N_REPEATS)) {
      
      folds <- make_stratified_folds(
        cv_data$gleason_ordinal,
        K = K,
        seed = 20260822 + rep_id
      )
      
      for (fold_id in seq_len(K)) {
        
        test_idx <- which(folds == fold_id)
        train_idx <- which(folds != fold_id)
        
        train_df <- cv_data[train_idx, , drop = FALSE]
        test_df <- cv_data[test_idx, , drop = FALSE]
        
        age_model <- MASS::polr(
          gleason_ordinal ~ age_years,
          data = train_df,
          Hess = FALSE
        )
        
        prob <- predict(
          age_model,
          newdata = test_df,
          type = "probs"
        )
        
        age_only_results[[counter]] <- data.frame(
          Repeat = rep_id,
          Fold = fold_id,
          LogLoss = multiclass_logloss(
            prob,
            test_df$gleason_ordinal
          ),
          MAE = ordinal_mae(
            prob,
            test_df$gleason_ordinal
          )
        )
        
        counter <- counter + 1
      }
    }
    
    age_only_results <- dplyr::bind_rows(
      age_only_results
    )
    
    cat("\nAge-only mean log-loss:\n")
    print(
      mean(age_only_results$LogLoss)
    )
    
    cat("\nAge-only mean MAE:\n")
    print(
      mean(age_only_results$MAE)
    )
    
    cat("\nCurrent signature mean log-loss:\n")
    print(
      mean(
        cv_results$LogLoss[
          cv_results$Combination == current_name
        ]
      )
    )
    
    cat("\nCurrent signature mean MAE:\n")
    print(
      mean(
        cv_results$MAE[
          cv_results$Combination == current_name
        ]
      )
    )
    ############################################################
    # FINAL SIGNATURE SELECTION CHECK
    #
    # Compare:
    #   1. Age only
    #   2. Current 3-gene PCA + age
    #   3. Five-gene PCA + age
    #
    # Same 5-fold x 5-repeat CV splits for every model
    # PCA/scaling are fitted inside each training fold
    ############################################################
    
    set.seed(20260822)
    
    library(MASS)
    library(dplyr)
    
    ############################################################
    # 1. DEFINE CURRENT AND FIVE-GENE SIGNATURES
    ############################################################
    
    genes_3 <- c(
      "MYLK",
      "PRKCB",
      "PDGFC"
    )
    
    genes_5 <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    ############################################################
    # 2. MAKE ONE FIXED SET OF CV SPLITS
    ############################################################
    
    K <- 5
    N_REPEATS <- 5
    
    fold_list <- vector(
      "list",
      N_REPEATS
    )
    
    for (rep_id in seq_len(N_REPEATS)) {
      
      fold_list[[rep_id]] <- make_stratified_folds(
        cv_data$gleason_ordinal,
        K = K,
        seed = 20260822 + rep_id
      )
    }
    
    ############################################################
    # 3. FUNCTION TO RUN A PCA SIGNATURE MODEL
    ############################################################
    
    run_pca_cv_model <- function(
    genes,
    model_name
    ) {
      
      results <- list()
      counter <- 1
      
      for (rep_id in seq_len(N_REPEATS)) {
        
        folds <- fold_list[[rep_id]]
        
        for (fold_id in seq_len(K)) {
          
          test_idx <- which(
            folds == fold_id
          )
          
          train_idx <- which(
            folds != fold_id
          )
          
          train_df <- cv_data[
            train_idx,
            ,
            drop = FALSE
          ]
          
          test_df <- cv_data[
            test_idx,
            ,
            drop = FALSE
          ]
          
          # PCA/scaling ONLY on training data
          sig_fit <- fit_fold_signature(
            train_df = train_df,
            test_df = test_df,
            genes = genes
          )
          
          train_df$cv_signature <-
            sig_fit$train_signature
          
          test_df$cv_signature <-
            sig_fit$test_signature
          
          # Ordinal model
          model <- MASS::polr(
            gleason_ordinal ~
              cv_signature +
              age_years,
            data = train_df,
            Hess = FALSE
          )
          
          # Held-out probabilities
          prob <- predict(
            model,
            newdata = test_df,
            type = "probs"
          )
          
          results[[counter]] <- data.frame(
            Model = model_name,
            Repeat = rep_id,
            Fold = fold_id,
            LogLoss = multiclass_logloss(
              prob,
              test_df$gleason_ordinal
            ),
            MAE = ordinal_mae(
              prob,
              test_df$gleason_ordinal
            )
          )
          
          counter <- counter + 1
        }
      }
      
      dplyr::bind_rows(results)
    }
    
    ############################################################
    # 4. RUN 3-GENE MODEL
    ############################################################
    
    cat(
      "\n===== RUNNING 3-GENE MODEL =====\n"
    )
    
    cv_3gene <- run_pca_cv_model(
      genes = genes_3,
      model_name = "3-gene PCA + age"
    )
    
    ############################################################
    # 5. RUN 5-GENE MODEL
    ############################################################
    
    cat(
      "\n===== RUNNING 5-GENE MODEL =====\n"
    )
    
    cv_5gene <- run_pca_cv_model(
      genes = genes_5,
      model_name = "5-gene PCA + age"
    )
    
    ############################################################
    # 6. RUN AGE-ONLY MODEL ON SAME FOLDS
    ############################################################
    
    cat(
      "\n===== RUNNING AGE-ONLY MODEL =====\n"
    )
    
    age_results <- list()
    counter <- 1
    
    for (rep_id in seq_len(N_REPEATS)) {
      
      folds <- fold_list[[rep_id]]
      
      for (fold_id in seq_len(K)) {
        
        test_idx <- which(
          folds == fold_id
        )
        
        train_idx <- which(
          folds != fold_id
        )
        
        train_df <- cv_data[
          train_idx,
          ,
          drop = FALSE
        ]
        
        test_df <- cv_data[
          test_idx,
          ,
          drop = FALSE
        ]
        
        model <- MASS::polr(
          gleason_ordinal ~ age_years,
          data = train_df,
          Hess = FALSE
        )
        
        prob <- predict(
          model,
          newdata = test_df,
          type = "probs"
        )
        
        age_results[[counter]] <- data.frame(
          Model = "Age only",
          Repeat = rep_id,
          Fold = fold_id,
          LogLoss = multiclass_logloss(
            prob,
            test_df$gleason_ordinal
          ),
          MAE = ordinal_mae(
            prob,
            test_df$gleason_ordinal
          )
        )
        
        counter <- counter + 1
      }
    }
    
    cv_age <- dplyr::bind_rows(
      age_results
    )
    
    ############################################################
    # 7. COMBINE ALL MODELS
    ############################################################
    
    final_cv <- dplyr::bind_rows(
      cv_age,
      cv_3gene,
      cv_5gene
    )
    
    ############################################################
    # 8. SUMMARY
    ############################################################
    
    final_cv_summary <- final_cv %>%
      dplyr::group_by(Model) %>%
      dplyr::summarise(
        Mean_LogLoss = mean(LogLoss),
        SD_LogLoss = sd(LogLoss),
        Mean_MAE = mean(MAE),
        SD_MAE = sd(MAE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        Mean_LogLoss
      )
    
    cat(
      "\n====================================================\n"
    )
    
    cat(
      "FINAL MODEL COMPARISON\n"
    )
    
    cat(
      "====================================================\n"
    )
    
    print(
      final_cv_summary,
      row.names = FALSE
    )
    
    ############################################################
    # 9. PAIRWISE COMPARISONS USING THE SAME FOLDS
    ############################################################
    
    wide_logloss <- final_cv %>%
      dplyr::select(
        Model,
        Repeat,
        Fold,
        LogLoss
      ) %>%
      tidyr::pivot_wider(
        names_from = Model,
        values_from = LogLoss
      )
    
    wide_mae <- final_cv %>%
      dplyr::select(
        Model,
        Repeat,
        Fold,
        MAE
      ) %>%
      tidyr::pivot_wider(
        names_from = Model,
        values_from = MAE
      )
    
    ############################################################
    # 10. 3-GENE VS 5-GENE
    ############################################################
    
    delta_3_vs_5 <- wide_logloss$`3-gene PCA + age` -
      wide_logloss$`5-gene PCA + age`
    
    cat(
      "\n===== 3-GENE VS 5-GENE LOG-LOSS =====\n"
    )
    
    cat(
      "Mean difference (3 - 5) = ",
      mean(delta_3_vs_5),
      "\n",
      sep = ""
    )
    
    print(
      t.test(
        delta_3_vs_5,
        mu = 0
      )
    )
    
    print(
      wilcox.test(
        wide_logloss$`3-gene PCA + age`,
        wide_logloss$`5-gene PCA + age`,
        paired = TRUE,
        exact = FALSE
      )
    )
    
    ############################################################
    # 11. 3-GENE VS AGE ONLY
    ############################################################
    
    delta_3_vs_age <- wide_logloss$`3-gene PCA + age` -
      wide_logloss$`Age only`
    
    cat(
      "\n===== 3-GENE VS AGE-ONLY LOG-LOSS =====\n"
    )
    
    cat(
      "Mean difference (3 - age) = ",
      mean(delta_3_vs_age),
      "\n",
      sep = ""
    )
    
    print(
      t.test(
        delta_3_vs_age,
        mu = 0
      )
    )
    
    print(
      wilcox.test(
        wide_logloss$`3-gene PCA + age`,
        wide_logloss$`Age only`,
        paired = TRUE,
        exact = FALSE
      )
    )
    
    ############################################################
    # 12. 5-GENE VS AGE ONLY
    ############################################################
    
    delta_5_vs_age <- wide_logloss$`5-gene PCA + age` -
      wide_logloss$`Age only`
    
    cat(
      "\n===== 5-GENE VS AGE-ONLY LOG-LOSS =====\n"
    )
    
    cat(
      "Mean difference (5 - age) = ",
      mean(delta_5_vs_age),
      "\n",
      sep = ""
    )
    
    print(
      t.test(
        delta_5_vs_age,
        mu = 0
      )
    )
    
    print(
      wilcox.test(
        wide_logloss$`5-gene PCA + age`,
        wide_logloss$`Age only`,
        paired = TRUE,
        exact = FALSE
      )
    )
    
    ############################################################
    # 13. MAE COMPARISON
    ############################################################
    
    cat(
      "\n===== MAE SUMMARY =====\n"
    )
    
    print(
      final_cv_summary %>%
        dplyr::select(
          Model,
          Mean_MAE,
          SD_MAE
        ),
      row.names = FALSE
    )
    
    ############################################################
    # 14. SAVE EVERYTHING
    ############################################################
    
    write.csv(
      final_cv,
      "final_3_vs_5_vs_age_CV_raw.csv",
      row.names = FALSE
    )
    
    write.csv(
      final_cv_summary,
      "final_3_vs_5_vs_age_CV_summary.csv",
      row.names = FALSE
    )
    
    cat(
      "\n===== FINAL COMPARISON COMPLETE =====\n"
    )
    ############################################################
    # EXTERNAL GEO VALIDATION:
    # 3-GENE VS 5-GENE FIXED TCGA-DERIVED SIGNATURES
    ############################################################
    
    # Five-gene signature was NOT selected using GEO.
    # We are only evaluating external performance.
    
    genes_5 <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    # Fit five-gene PCA on the entire TCGA discovery cohort
    tcga_5x <- as.matrix(
      cv_data[, genes_5]
    )
    
    tcga_5x_scaled <- scale(tcga_5x)
    
    tcga_5_pca <- prcomp(
      tcga_5x_scaled,
      center = FALSE,
      scale. = FALSE
    )
    
    weights_5 <- tcga_5_pca$rotation[, 1]
    
    # Fix arbitrary PCA sign
    anchor <- which.max(abs(weights_5))
    
    if (weights_5[anchor] < 0) {
      weights_5 <- -weights_5
    }
    
    tcga_5_means <- attr(
      tcga_5x_scaled,
      "scaled:center"
    )
    
    tcga_5_sds <- attr(
      tcga_5x_scaled,
      "scaled:scale"
    )
    
    cat("\n===== FIVE-GENE PCA WEIGHTS =====\n")
    
    print(
      data.frame(
        Gene = names(weights_5),
        Weight = as.numeric(weights_5),
        Mean = as.numeric(tcga_5_means),
        SD = as.numeric(tcga_5_sds)
      )
    )
    
    ############################################################
    # APPLY FIXED FIVE-GENE PARAMETERS TO GEO
    ############################################################
    
    geo_5_matrix <- as.matrix(
      geo_clean[, genes_5]
    )
    
    geo_5_scaled <- scale(
      geo_5_matrix,
      center = tcga_5_means,
      scale = tcga_5_sds
    )
    
    geo_5_signature <- as.numeric(
      geo_5_scaled %*% weights_5
    )
    
    geo_5 <- geo_clean
    
    geo_5$signature_5gene <- geo_5_signature
    
    ############################################################
    # GEO 5-GENE ORDINAL MODEL
    ############################################################
    
    geo_5$gleason_ordinal <- factor(
      geo_5$gleason_score,
      levels = c(6, 7, 8, 9),
      ordered = TRUE
    )
    
    geo_5_model <- MASS::polr(
      gleason_ordinal ~
        signature_5gene +
        age_years,
      data = geo_5,
      Hess = TRUE
    )
    
    cat("\n===== GEO 5-GENE ORDINAL MODEL =====\n")
    
    print(
      summary(geo_5_model)
    )
    
    co5 <- coef(summary(geo_5_model))
    
    beta5 <- co5[
      "signature_5gene",
      "Value"
    ]
    
    se5 <- co5[
      "signature_5gene",
      "Std. Error"
    ]
    
    p5 <- 2 * pnorm(
      abs(
        co5[
          "signature_5gene",
          "t value"
        ]
      ),
      lower.tail = FALSE
    )
    
    geo_5_result <- data.frame(
      N = nrow(geo_5),
      OR = exp(beta5),
      CI_lower = exp(beta5 - 1.96 * se5),
      CI_upper = exp(beta5 + 1.96 * se5),
      P_value = p5
    )
    
    print(geo_5_result)
    
    ############################################################
    # CURRENT 3-GENE GEO RESULT FOR COMPARISON
    ############################################################
    
    cat("\n===== GEO CURRENT 3-GENE RESULT =====\n")
    
    print(
      geo_clean_gleason_result
    )
    print(p_gleason_final_v2)
    ############################################################
    # FINAL FIGURE 2 — SIGNATURE VS GLEASON SCORE
    # Publication-ready version
    ############################################################
    
    library(ggplot2)
    library(patchwork)
    
    ############################################################
    # GEO
    ############################################################
    
    geo_fig_final <- ggplot(
      geo_clean,
      aes(
        x = factor(gleason_score),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        alpha = 0.35,
        size = 1.5
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0(
          "Ordinal logistic regression\n",
          "OR = 0.226 (95% CI 0.095–0.533)\n",
          "P = 0.000693"
        ),
        hjust = 1.05,
        vjust = 1.20,
        size = 3.5
      ) +
      labs(
        title = "A  GSE70768 validation cohort",
        x = "Gleason score",
        y = "Three-gene signature"
      ) +
      theme_classic(
        base_size = 12
      )
    
    ############################################################
    # TCGA
    ############################################################
    
    tcga_fig_final <- ggplot(
      tcga_signature_clean,
      aes(
        x = factor(gleason_score),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        alpha = 0.25,
        size = 1.3
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0(
          "Ordinal logistic regression\n",
          "OR = 0.796 (95% CI 0.713–0.890)\n",
          "P = 5.83 × 10⁻⁵"
        ),
        hjust = 1.05,
        vjust = 1.20,
        size = 3.5
      ) +
      labs(
        title = "B  TCGA-PRAD discovery cohort",
        x = "Gleason score",
        y = "Three-gene signature"
      ) +
      theme_classic(
        base_size = 12
      )
    
    ############################################################
    # COMBINE
    ############################################################
    
    figure2_final <-
      geo_fig_final +
      tcga_fig_final +
      plot_annotation(
        title =
          "Association of the MYLK–PRKCB–PDGFC signature with Gleason score"
      )
    
    print(figure2_final)
    
    ############################################################
    # SAVE HIGH-RESOLUTION FIGURE
    ############################################################
    
    ggsave(
      "Figure_2_Gleason_signature.png",
      figure2_final,
      width = 11,
      height = 5.8,
      units = "in",
      dpi = 600
    )    
    ############################################################
    # FINAL FIGURE 2 — COSMETIC POLISH
    ############################################################
    
    library(ggplot2)
    library(patchwork)
    
    ############################################################
    # A — GEO VALIDATION
    ############################################################
    
    geo_fig_final <- ggplot(
      geo_clean,
      aes(
        x = factor(gleason_score),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        alpha = 0.35,
        size = 1.5
      ) +
      annotate(
        "text",
        x = 4.15,
        y = max(geo_clean$signature_clean) + 0.15,
        label = paste0(
          "Ordinal logistic regression\n",
          "OR = 0.226 (95% CI 0.095–0.533)\n",
          "P = 0.000693"
        ),
        hjust = 1,
        vjust = 1,
        size = 3.2
      ) +
      labs(
        title = "A  GSE70768 validation cohort",
        x = "Gleason score",
        y = "MYLK–PRKCB–PDGFC signature score"
      ) +
      theme_classic(
        base_size = 12
      )
    
    ############################################################
    # B — TCGA DISCOVERY
    ############################################################
    
    tcga_fig_final <- ggplot(
      tcga_signature_clean,
      aes(
        x = factor(gleason_score),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        alpha = 0.25,
        size = 1.3
      ) +
      annotate(
        "text",
        x = 5.15,
        y = max(tcga_signature_clean$signature_clean) + 0.15,
        label = paste0(
          "Ordinal logistic regression\n",
          "OR = 0.796 (95% CI 0.713–0.890)\n",
          "P = 5.83 × 10⁻⁵"
        ),
        hjust = 1,
        vjust = 1,
        size = 3.2
      ) +
      labs(
        title = "B  TCGA-PRAD discovery cohort",
        x = "Gleason score",
        y = "MYLK–PRKCB–PDGFC signature score"
      ) +
      theme_classic(
        base_size = 12
      )
    
    ############################################################
    # COMBINE
    ############################################################
    
    figure2_final <-
      geo_fig_final +
      tcga_fig_final +
      plot_annotation(
        title =
          "Association of the MYLK–PRKCB–PDGFC signature with Gleason score"
      )
    
    print(figure2_final)
    
    ############################################################
    # SAVE PUBLICATION-QUALITY FILE
    ############################################################
    
    ggsave(
      "Figure_2_Gleason_signature_final.png",
      figure2_final,
      width = 11,
      height = 5.8,
      units = "in",
      dpi = 600
    )
    
    cat("\n===== FIGURE 2 SAVED =====\n")
    cat("Figure_2_Gleason_signature_final.png\n")
    ############################################################
    # FINAL FIGURE 3 — T-STAGE ASSOCIATION
    ############################################################
    
    library(ggplot2)
    library(patchwork)
    
    ############################################################
    # GEO
    ############################################################
    
    geo_T_plot_final <- ggplot(
      geo_T_clean,
      aes(
        x = factor(
          T_advanced_clean,
          levels = c(0, 1),
          labels = c("T2", "T3/T4")
        ),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        alpha = 0.35,
        size = 1.5
      ) +
      annotate(
        "text",
        x = 2.35,
        y = max(geo_T_clean$signature_clean) + 0.15,
        label = paste0(
          "Wilcoxon P = 0.2956\n",
          "Adjusted OR = 2.28 (95% CI 0.98–5.29)\n",
          "P = 0.0546"
        ),
        hjust = 1,
        vjust = 1,
        size = 3.2
      ) +
      labs(
        title = "A  GSE70768 validation cohort",
        x = "Pathological T stage",
        y = "MYLK–PRKCB–PDGFC signature score"
      ) +
      theme_classic(
        base_size = 12
      )
    
    ############################################################
    # TCGA
    ############################################################
    
    tcga_T_plot_final <- ggplot(
      tcga_T_clean,
      aes(
        x = factor(
          T_advanced_clean,
          levels = c(0, 1),
          labels = c("T2", "T3/T4")
        ),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.10,
        height = 0,
        alpha = 0.20,
        size = 1.3
      ) +
      annotate(
        "text",
        x = 2.35,
        y = max(tcga_T_clean$signature_clean) + 0.15,
        label = paste0(
          "Wilcoxon P = 0.0521\n",
          "Adjusted OR = 0.960 (95% CI 0.833–1.107)\n",
          "P = 0.577"
        ),
        hjust = 1,
        vjust = 1,
        size = 3.2
      ) +
      labs(
        title = "B  TCGA-PRAD discovery cohort",
        x = "Pathological T stage",
        y = "MYLK–PRKCB–PDGFC signature score"
      ) +
      theme_classic(
        base_size = 12
      )
    
    ############################################################
    # COMBINE
    ############################################################
    
    figure3_final <-
      geo_T_plot_final +
      tcga_T_plot_final +
      plot_annotation(
        title =
          "Association of the MYLK–PRKCB–PDGFC signature with pathological T stage"
      )
    
    print(figure3_final)
    
    ############################################################
    # SAVE HIGH-RESOLUTION FIGURE
    ############################################################
    
    ggsave(
      "Figure_3_T_stage_signature_final.png",
      figure3_final,
      width = 11,
      height = 5.8,
      units = "in",
      dpi = 600
    )
    
    cat("\n===== FIGURE 3 SAVED =====\n")
    cat("Figure_3_T_stage_signature_final.png\n")
    ############################################################
    # FINAL FIGURE 3 POLISH
    ############################################################
    
    geo_T_plot_final <- ggplot(
      geo_T_clean,
      aes(
        x = factor(
          T_advanced_clean,
          levels = c(0, 1),
          labels = c("T2", "T3/T4")
        ),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        alpha = 0.35,
        size = 1.5
      ) +
      annotate(
        "text",
        x = 2.35,
        y = max(geo_T_clean$signature_clean) + 0.15,
        label = paste0(
          "Wilcoxon P = 0.2956\n",
          "Adjusted OR = 2.28 (95% CI 0.98–5.29)\n",
          "P = 0.0546"
        ),
        hjust = 1,
        vjust = 1,
        size = 3.0
      ) +
      labs(
        title = "A  GSE70768 validation cohort",
        x = "Pathological T stage",
        y = "MYLK–PRKCB–PDGFC signature"
      ) +
      theme_classic(base_size = 12)
    
    tcga_T_plot_final <- ggplot(
      tcga_T_clean,
      aes(
        x = factor(
          T_advanced_clean,
          levels = c(0, 1),
          labels = c("T2", "T3/T4")
        ),
        y = signature_clean
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.60,
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.10,
        height = 0,
        alpha = 0.20,
        size = 1.3
      ) +
      annotate(
        "text",
        x = 2.35,
        y = max(tcga_T_clean$signature_clean) + 0.15,
        label = paste0(
          "Wilcoxon P = 0.0521\n",
          "Adjusted OR = 0.960 (95% CI 0.833–1.107)\n",
          "P = 0.577"
        ),
        hjust = 1,
        vjust = 1,
        size = 3.0
      ) +
      labs(
        title = "B  TCGA-PRAD discovery cohort",
        x = "Pathological T stage",
        y = "MYLK–PRKCB–PDGFC signature"
      ) +
      theme_classic(base_size = 12)
    
    figure3_final <-
      geo_T_plot_final +
      tcga_T_plot_final +
      plot_annotation(
        title =
          "Association of the MYLK–PRKCB–PDGFC signature with pathological T stage"
      )
    
    print(figure3_final)
    
    ggsave(
      "Figure_3_T_stage_signature_final.png",
      figure3_final,
      width = 11,
      height = 5.8,
      units = "in",
      dpi = 600
    )
    cat("\n===== WORKSPACE CHECK =====\n")
    
    objects_to_check <- c(
      "tcga_signature_clean",
      "geo_clean",
      "signature_go_final",
      "signature_kegg_final",
      "signature_gsea_go",
      "signature_gsea_kegg",
      "gene_list_signature",
      "gene_list_signature_no3",
      "clean_weights",
      "clean_means",
      "clean_sds",
      "tcga_clean_gleason_result",
      "geo_clean_gleason_result",
      "tcga_T_clean",
      "geo_T_clean",
      "tcga_clean_cox",
      "tcga_clean_cox_adj",
      "tcga_clean_ph",
      "cv_summary",
      "final_cv_summary"
    )
    
    print(
      data.frame(
        Object = objects_to_check,
        Exists = sapply(objects_to_check, exists)
      )
    )
    cat("\n===== GSEA RESULT CHECK =====\n")
    
    if (exists("signature_go_final")) {
      print(
        signature_go_final[
          order(signature_go_final$p.adjust),
          c("Description", "NES", "pvalue", "p.adjust")
        ][1:10, ]
      )
    }
    
    if (exists("signature_kegg_final")) {
      print(
        signature_kegg_final[
          order(signature_kegg_final$p.adjust),
          c("Description", "NES", "pvalue", "p.adjust")
        ][1:10, ]
      )
    }
    load("TCGA_PRAD_GSEA_progress.RData")
    cat("\n===== GSEA PROGRESS WORKSPACE =====\n")
    
    objects_to_check <- c(
      "tcga_signature_clean",
      "geo_clean",
      "signature_go_final",
      "signature_kegg_final",
      "signature_gsea_go",
      "signature_gsea_kegg",
      "gene_list_signature",
      "gene_list_signature_no3",
      "clean_weights",
      "clean_means",
      "clean_sds",
      "tcga_clean_gleason_result",
      "geo_clean_gleason_result",
      "tcga_T_clean",
      "geo_T_clean",
      "tcga_clean_cox",
      "tcga_clean_cox_adj",
      "tcga_clean_ph"
    )
    
    print(
      data.frame(
        Object = objects_to_check,
        Exists = sapply(objects_to_check, exists)
      )
    )    
    cat("\n===== SAVED CV RESULTS =====\n")
    
    print(
      read.csv(
        "final_3_vs_5_vs_age_CV_summary.csv"
      )
    )
    
    cat("\n===== SAVED 3-GENE CV RESULTS =====\n")
    
    print(
      read.csv(
        "three_gene_CV_summary.csv"
      )
    )
    cat("\n===== TCGA_PRAD_workspace.RData OBJECTS =====\n")
    
    tmp1 <- new.env()
    load(
      "TCGA_PRAD_workspace.RData",
      envir = tmp1
    )
    
    print(
      sort(ls(tmp1))
    )
    
    cat("\n===== TCGA_PRAD_GSEA_progress.RData OBJECTS =====\n")
    
    tmp2 <- new.env()
    load(
      "TCGA_PRAD_GSEA_progress.RData",
      envir = tmp2
    )
    
    print(
      sort(ls(tmp2))
    )
    cat("\n===== RESULTS.RDS =====\n")
    
    results_obj <- readRDS("results.rds")
    
    cat("Class:\n")
    print(class(results_obj))
    
    cat("\nStructure:\n")
    str(
      results_obj,
      max.level = 2
    )
    
    cat("\nNames if available:\n")
    if (!is.null(names(results_obj))) {
      print(names(results_obj))
    }
    
    
    cat("\n===== DF.RDS =====\n")
    
    df_obj <- readRDS("df.rds")
    
    cat("Class:\n")
    print(class(df_obj))
    
    cat("\nStructure:\n")
    str(
      df_obj,
      max.level = 2
    )
    
    cat("\nNames if available:\n")
    if (!is.null(names(df_obj))) {
      print(names(df_obj))
    }
    load("TCGA_PRAD_GSEA_progress.RData")
    cat("\n===== RECOVERY CHECK =====\n")
    
    print(
      data.frame(
        Object = c(
          "dds",
          "counts",
          "tumor_counts",
          "tumor_samples_unique",
          "sample_info",
          "res",
          "col_data"
        ),
        Exists = sapply(
          c(
            "dds",
            "counts",
            "tumor_counts",
            "tumor_samples_unique",
            "sample_info",
            "res",
            "col_data"
          ),
          exists
        )
      )
    )
    
    cat("\nDimensions where available:\n")
    
    for (x in c("counts", "tumor_counts", "sample_info", "col_data")) {
      if (exists(x)) {
        cat("\n", x, ":\n", sep = "")
        print(dim(get(x)))
      }
    }    
    cat("\n===== RECOVERY OBJECT STRUCTURE =====\n")
    
    cat("\n--- tumor_samples_unique ---\n")
    print(class(tumor_samples_unique))
    print(length(tumor_samples_unique))
    print(head(tumor_samples_unique))
    
    cat("\n--- sample_info ---\n")
    print(class(sample_info))
    print(colnames(sample_info))
    print(head(sample_info))
    
    cat("\n--- col_data ---\n")
    print(class(col_data))
    print(colnames(col_data))
    print(head(col_data))
    
    cat("\n--- counts column names ---\n")
    print(head(colnames(counts), 10))
    
    cat("\n--- tumor_counts ---\n")
    print(class(tumor_counts))
    print(length(tumor_counts))
    print(head(tumor_counts))
    
    cat("\n--- dds ---\n")
    print(class(dds))
    print(colnames(SummarizedExperiment::colData(dds)))
    print(head(SummarizedExperiment::colData(dds)))
    
    cat("\n--- res ---\n")
    print(class(res))
    print(head(as.data.frame(res)))
    ############################################################
    # RECOVER CLEAN TCGA ANALYSIS FROM SAVED RNA-SEQ OBJECTS
    ############################################################
    
    library(DESeq2)
    library(dplyr)
    library(MASS)
    
    cat("\n===== 1. BUILD TUMOR EXPRESSION MATRIX =====\n")
    
    # tumor_samples_unique contains the one-primary-tumor-per-patient mapping
    tumor_map <- tumor_samples_unique
    
    cat("Tumor patients =", nrow(tumor_map), "\n")
    
    # Extract the corresponding VST expression from the saved DESeq2 object
    vsd_recovered <- DESeq2::vst(
      dds,
      blind = FALSE
    )
    
    # Keep exactly the one primary-tumor sample per patient
    tumor_sample_ids <- tumor_map$sample_id
    
    vsd_tumor_recovered <- assay(vsd_recovered)[
      ,
      match(
        tumor_sample_ids,
        colnames(vsd_recovered)
      ),
      drop = FALSE
    ]
    
    cat(
      "Recovered VST dimensions = ",
      paste(dim(vsd_tumor_recovered), collapse = " x "),
      "\n",
      sep = ""
    )
    
    cat(
      "Matched samples = ",
      sum(
        tumor_sample_ids %in%
          colnames(vsd_tumor_recovered)
      ),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 2. CONVERT GENE IDS AND EXTRACT THE FIVE CANDIDATES
    ############################################################
    
    candidate_genes <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    # Extract clean Ensembl IDs
    recovered_ensembl <- sub(
      "\\..*$",
      "",
      rownames(vsd_tumor_recovered)
    )
    
    # Use existing annotation from res where available
    res_df <- as.data.frame(res)
    
    res_df$ensembl_clean <- sub(
      "\\..*$",
      "",
      rownames(res_df)
    )
    
    candidate_map <- res_df %>%
      dplyr::filter(
        gene_symbol %in% candidate_genes
      ) %>%
      dplyr::select(
        gene_symbol,
        ensembl_clean
      ) %>%
      dplyr::distinct()
    
    cat("\n===== CANDIDATE GENE MAPPING =====\n")
    print(candidate_map)
    
    
    ############################################################
    # 3. EXTRACT CANDIDATE EXPRESSION
    ############################################################
    
    candidate_expr_recovered <- matrix(
      NA_real_,
      nrow = nrow(tumor_map),
      ncol = length(candidate_genes)
    )
    
    colnames(candidate_expr_recovered) <- candidate_genes
    
    for (g in candidate_genes) {
      
      gene_ids <- candidate_map$ensembl_clean[
        candidate_map$gene_symbol == g
      ]
      
      gene_idx <- which(
        recovered_ensembl %in% gene_ids
      )
      
      if (length(gene_idx) == 0) {
        stop(
          paste(
            "Could not find expression for",
            g
          )
        )
      }
      
      # In case multiple Ensembl rows map to the gene,
      # use the first mapped row, matching the existing workflow.
      candidate_expr_recovered[, g] <-
        as.numeric(
          vsd_tumor_recovered[
            gene_idx[1],
            ,
            drop = TRUE
          ]
        )
    }
    
    candidate_expr_recovered <- as.data.frame(
      candidate_expr_recovered
    )
    
    candidate_expr_recovered$patient_id <-
      tumor_map$patient_id
    
    candidate_expr_recovered$sample_id <-
      tumor_map$sample_id
    
    cat("\n===== CANDIDATE EXPRESSION CHECK =====\n")
    print(
      head(
        candidate_expr_recovered[
          ,
          c(
            "patient_id",
            "sample_id",
            candidate_genes
          )
        ]
      )
    )
    
    
    ############################################################
    # 4. BUILD TCGA CLINICAL DATA
    ############################################################
    
    # Recover clinical information from the original metadata
    # contained in sample_info / df.rds if needed.
    
    # df.rds contains clinical/sample metadata, but only a few rows.
    # Therefore use the clinical object that existed in the original
    # analysis if available; otherwise reconstruct from the DESeq2
    # / downloaded clinical sources.
    
    cat("\n===== CLINICAL OBJECT CHECK =====\n")
    
    print(
      data.frame(
        Object = c(
          "clinical_prad",
          "gleason_analysis"
        ),
        Exists = c(
          exists("clinical_prad"),
          exists("gleason_analysis")
        )
      )
    )
    ############################################################
    # RECOVER TCGA CLINICAL DATA + REBUILD BASE ANALYSIS OBJECTS
    ############################################################
    
    library(TCGAbiolinks)
    library(dplyr)
    library(DESeq2)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    
    cat("\n===== 1. DOWNLOAD TCGA-PRAD CLINICAL DATA =====\n")
    
    clinical <- TCGAbiolinks::GDCquery_clinic(
      project = "TCGA-PRAD",
      type = "clinical"
    )
    
    cat("Clinical rows =", nrow(clinical), "\n")
    
    cat("\nClinical columns needed:\n")
    print(
      intersect(
        c(
          "submitter_id",
          "age_at_diagnosis",
          "gleason_score",
          "ajcc_pathologic_t",
          "ajcc_pathologic_n",
          "vital_status",
          "days_to_last_follow_up"
        ),
        colnames(clinical)
      )
    )
    
    
    ############################################################
    # 2. REBUILD clinical_prad EXACTLY AS ORIGINAL SCRIPT
    ############################################################
    
    clinical_prad <- clinical[
      clinical$submitter_id %in%
        tumor_samples_unique$patient_id,
      c(
        "submitter_id",
        "age_at_diagnosis",
        "gleason_score",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "vital_status",
        "days_to_last_follow_up"
      )
    ]
    
    cat(
      "\nclinical_prad N = ",
      nrow(clinical_prad),
      "\n",
      sep = ""
    )
    
    cat(
      "Patients matched to tumor_samples_unique = ",
      sum(
        tumor_samples_unique$patient_id %in%
          clinical_prad$submitter_id
      ),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 3. REBUILD VST TUMOR MATRIX
    ############################################################
    
    vsd <- DESeq2::vst(
      dds,
      blind = FALSE
    )
    
    vsd_tumor <- assay(
      vsd
    )[
      ,
      tumor_samples_unique$sample_id,
      drop = FALSE
    ]
    
    cat(
      "\nVST tumor dimensions = ",
      paste(dim(vsd_tumor), collapse = " x "),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 4. REBUILD THE 11-GENE EXPRESSION TABLE
    ############################################################
    
    candidate_genes_11 <- c(
      "MYLK",
      "PRKCB",
      "PDGFD",
      "PRKCA",
      "PDGFC",
      "VEGFD",
      "MET",
      "MYLK4",
      "KDR",
      "EGFR",
      "PDGFRA"
    )
    
    candidate_map_11 <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = candidate_genes_11,
      keytype = "SYMBOL",
      columns = c(
        "SYMBOL",
        "ENSEMBL"
      )
    )
    
    candidate_map_11 <- candidate_map_11[
      !duplicated(candidate_map_11$SYMBOL),
    ]
    
    vsd_ensembl <- sub(
      "\\..*$",
      "",
      rownames(vsd_tumor)
    )
    
    candidate_idx <- match(
      candidate_map_11$ENSEMBL,
      vsd_ensembl
    )
    
    cat(
      "\nMissing candidate genes = ",
      sum(is.na(candidate_idx)),
      "\n",
      sep = ""
    )
    
    candidate_expr <- vsd_tumor[
      candidate_idx,
      ,
      drop = FALSE
    ]
    
    rownames(candidate_expr) <-
      candidate_map_11$SYMBOL
    
    
    ############################################################
    # 5. REBUILD gleason_analysis
    ############################################################
    
    tumor_patient_ids <-
      tumor_samples_unique$patient_id
    
    tumor_sample_ids <-
      tumor_samples_unique$sample_id
    
    gleason_df <- clinical_prad[
      ,
      c(
        "submitter_id",
        "gleason_score"
      )
    ]
    
    gleason_df <- gleason_df[
      match(
        tumor_patient_ids,
        gleason_df$submitter_id
      ),
      ,
      drop = FALSE
    ]
    
    gleason_expr <- t(candidate_expr)
    
    gleason_analysis <- data.frame(
      patient_id = tumor_patient_ids,
      sample_id = tumor_sample_ids,
      gleason_score = gleason_df$gleason_score,
      gleason_expr
    )
    
    # Add age
    age_lookup <- clinical_prad[
      ,
      c(
        "submitter_id",
        "age_at_diagnosis"
      )
    ]
    
    age_lookup$age_years <-
      age_lookup$age_at_diagnosis / 365.25
    
    age_lookup <- age_lookup[
      ,
      c(
        "submitter_id",
        "age_years"
      )
    ]
    
    gleason_analysis <- gleason_analysis %>%
      dplyr::left_join(
        age_lookup,
        by = c(
          "patient_id" =
            "submitter_id"
        )
      )
    
    cat(
      "\n===== RECOVERED gleason_analysis =====\n"
    )
    
    cat(
      "N = ",
      nrow(gleason_analysis),
      "\n",
      sep = ""
    )
    
    cat(
      "Gleason distribution:\n"
    )
    
    print(
      table(
        gleason_analysis$gleason_score,
        useNA = "ifany"
      )
    )
    
    cat(
      "\nAge missing = ",
      sum(
        is.na(
          gleason_analysis$age_years
        )
      ),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 6. RECOVER THE FIVE-GENE DATASET
    ############################################################
    
    candidate_genes_5 <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    cv_source <- gleason_analysis %>%
      dplyr::select(
        patient_id,
        sample_id,
        gleason_score,
        age_years,
        dplyr::all_of(
          candidate_genes_5
        )
      )
    
    cv_data_recovered <- cv_source %>%
      dplyr::filter(
        complete.cases(.)
      )
    
    cat(
      "\n===== RECOVERED FIVE-GENE DATASET =====\n"
    )
    
    cat(
      "N = ",
      nrow(cv_data_recovered),
      "\n",
      sep = ""
    )
    
    print(
      table(
        cv_data_recovered$gleason_score,
        useNA = "ifany"
      )
    )
    
    
    ############################################################
    # 7. CHECK AGAINST OUR FINAL 455-PATIENT COHORT
    ############################################################
    
    cat(
      "\n===== FINAL RECOVERY CHECK =====\n"
    )
    
    cat(
      "Expected N = 455\n"
    )
    
    cat(
      "Recovered N = ",
      nrow(cv_data_recovered),
      "\n",
      sep = ""
    )
    
    if (nrow(cv_data_recovered) != 455) {
      warning(
        paste0(
          "Recovered cohort is ",
          nrow(cv_data_recovered),
          ", not 455. STOP before rebuilding the signature."
        )
      )
    }
    
    
    ############################################################
    # 8. SAVE THE RECOVERED BASE OBJECTS
    ############################################################
    
    save(
      clinical,
      clinical_prad,
      tumor_samples_unique,
      vsd_tumor,
      gleason_analysis,
      cv_data_recovered,
      file =
        "TCGA_PRAD_clean_recovery.RData"
    )
    
    cat(
      "\n===== CLEAN RECOVERY SAVED =====\n"
    )
    
    cat(
      "TCGA_PRAD_clean_recovery.RData\n"
    )
    ############################################################
    # REBUILD FINAL LOCKED 3-GENE SIGNATURE
    # from recovered TCGA data
    ############################################################
    
    library(dplyr)
    library(limma)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    
    ############################################################
    # 1. FINAL LOCKED SIGNATURE PARAMETERS
    ############################################################
    
    final_genes <- c(
      "MYLK",
      "PRKCB",
      "PDGFC"
    )
    
    final_weights <- c(
      MYLK  = 0.5904898,
      PRKCB = 0.5694109,
      PDGFC = 0.5719205
    )
    
    final_means <- c(
      MYLK  = 13.988972,
      PRKCB = 8.563813,
      PDGFC = 10.055332
    )
    
    final_sds <- c(
      MYLK  = 1.4614046,
      PRKCB = 1.2771604,
      PDGFC = 0.9769621
    )
    
    ############################################################
    # 2. CLEAN TCGA SIGNATURE COHORT
    ############################################################
    
    tcga_signature_clean <- cv_data_recovered %>%
      dplyr::select(
        patient_id,
        sample_id,
        gleason_score,
        age_years,
        dplyr::all_of(final_genes)
      ) %>%
      dplyr::filter(
        complete.cases(.)
      )
    
    cat("\n===== FINAL TCGA SIGNATURE COHORT =====\n")
    cat("N =", nrow(tcga_signature_clean), "\n")
    
    ############################################################
    # 3. RECREATE THE LOCKED SIGNATURE
    ############################################################
    
    signature_matrix <- scale(
      tcga_signature_clean[
        ,
        final_genes,
        drop = FALSE
      ],
      center = final_means,
      scale = final_sds
    )
    
    tcga_signature_clean$signature_clean <- as.numeric(
      signature_matrix %*% final_weights
    )
    
    cat("\nSignature summary:\n")
    print(
      summary(
        tcga_signature_clean$signature_clean
      )
    )
    
    ############################################################
    # 4. VERIFY AGAINST THE PREVIOUSLY ESTABLISHED SIGNATURE
    ############################################################
    
    cat("\n===== SIGNATURE PARAMETERS =====\n")
    
    print(
      data.frame(
        Gene = final_genes,
        Weight = as.numeric(final_weights),
        Mean = as.numeric(final_means),
        SD = as.numeric(final_sds)
      ),
      row.names = FALSE
    )
    
    ############################################################
    # 5. RECOVER THE FIVE-GENE PARAMETERS FOR COMPARISON
    ############################################################
    
    five_genes <- c(
      "MYLK",
      "PRKCB",
      "MYLK4",
      "PDGFC",
      "VEGFD"
    )
    
    five_weights <- c(
      MYLK  = 0.56844705,
      PRKCB = 0.53338323,
      MYLK4 = 0.32151403,
      PDGFC = 0.53480845,
      VEGFD = 0.05457954
    )
    
    five_means <- c(
      MYLK  = 13.988972,
      PRKCB = 8.563813,
      MYLK4 = 6.110112,
      PDGFC = 10.055332,
      VEGFD = 4.717336
    )
    
    five_sds <- c(
      MYLK  = 1.4614046,
      PRKCB = 1.2771604,
      MYLK4 = 0.5378170,
      PDGFC = 0.9769621,
      VEGFD = 0.5761807
    )
    
    ############################################################
    # 6. BUILD TRUE SIGNATURE-ASSOCIATED GSEA RANKING
    ############################################################
    
    cat("\n===== BUILDING SIGNATURE-ASSOCIATED GSEA =====\n")
    
    expr_clean <- vsd_tumor[
      ,
      match(
        tcga_signature_clean$sample_id,
        colnames(vsd_tumor)
      ),
      drop = FALSE
    ]
    
    cat(
      "Expression dimensions = ",
      paste(
        dim(expr_clean),
        collapse = " x "
      ),
      "\n",
      sep = ""
    )
    
    stopifnot(
      identical(
        colnames(expr_clean),
        tcga_signature_clean$sample_id
      )
    )
    
    ############################################################
    # Gene-wise association with signature, adjusted for age
    ############################################################
    
    design_clean <- model.matrix(
      ~ signature_clean + age_years,
      data = tcga_signature_clean
    )
    
    fit_clean <- limma::lmFit(
      expr_clean,
      design_clean
    )
    
    fit_clean <- limma::eBayes(
      fit_clean
    )
    
    signature_col <- which(
      colnames(design_clean) ==
        "signature_clean"
    )
    
    signature_t <- fit_clean$t[
      ,
      signature_col
    ]
    
    ############################################################
    # Ensembl -> Entrez
    ############################################################
    
    ensembl_clean <- sub(
      "\\..*$",
      "",
      rownames(expr_clean)
    )
    
    entrez_clean <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_clean,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    keep <- !is.na(
      entrez_clean
    ) &
      !is.na(
        signature_t
      )
    
    signature_gsea_input <- data.frame(
      ENTREZID = entrez_clean[keep],
      t_stat = signature_t[keep],
      stringsAsFactors = FALSE
    )
    
    signature_gsea_input <- signature_gsea_input[
      !duplicated(
        signature_gsea_input$ENTREZID
      ),
    ]
    
    gene_list_signature <- signature_gsea_input$t_stat
    
    names(
      gene_list_signature
    ) <- signature_gsea_input$ENTREZID
    
    gene_list_signature <- sort(
      gene_list_signature,
      decreasing = TRUE
    )
    
    cat(
      "Ranked genes = ",
      length(gene_list_signature),
      "\n",
      sep = ""
    )
    
    ############################################################
    # 7. REMOVE THE THREE SIGNATURE GENES FOR SENSITIVITY
    ############################################################
    
    signature_entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = final_genes,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    signature_entrez <- unname(
      signature_entrez[
        !is.na(signature_entrez)
      ]
    )
    
    gene_list_signature_no3 <-
      gene_list_signature[
        !names(gene_list_signature) %in%
          signature_entrez
      ]
    
    ############################################################
    # 8. SAVE THE RECOVERED FINAL OBJECTS
    ############################################################
    
    save(
      clinical,
      clinical_prad,
      tumor_samples_unique,
      vsd_tumor,
      gleason_analysis,
      cv_data_recovered,
      tcga_signature_clean,
      final_genes,
      final_weights,
      final_means,
      final_sds,
      gene_list_signature,
      gene_list_signature_no3,
      file =
        "TCGA_PRAD_FINAL_RECOVERED_OBJECTS.RData"
    )
    
    cat(
      "\n===== FINAL RECOVERY COMPLETE =====\n"
    )
    
    cat(
      "Saved: TCGA_PRAD_FINAL_RECOVERED_OBJECTS.RData\n"
    )
    ############################################################
    # FINAL SIGNATURE-ASSOCIATED GSEA
    # Recomputed from the recovered final 455-patient dataset
    ############################################################
    
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(patchwork)
    
    set.seed(20260823)
    
    ############################################################
    # 1. PRIMARY GO GSEA
    ############################################################
    
    cat("\n===== FINAL GO GSEA =====\n")
    
    signature_gsea_go_final <- clusterProfiler::gseGO(
      geneList = gene_list_signature,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      eps = 0,
      nPermSimple = 100000,
      verbose = FALSE
    )
    
    signature_go_final <- as.data.frame(
      signature_gsea_go_final
    )
    
    print(
      signature_go_final[
        order(signature_go_final$p.adjust),
        c(
          "Description",
          "NES",
          "pvalue",
          "p.adjust"
        )
      ][
        1:min(15, nrow(signature_go_final)),
        ,
        drop = FALSE
      ]
    )
    
    
    ############################################################
    # 2. PRIMARY KEGG GSEA
    ############################################################
    
    cat("\n===== FINAL KEGG GSEA =====\n")
    
    set.seed(20260823)
    
    signature_gsea_kegg_final <- clusterProfiler::gseKEGG(
      geneList = gene_list_signature,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0,
      nPermSimple = 100000,
      verbose = FALSE
    )
    
    signature_kegg_final <- as.data.frame(
      signature_gsea_kegg_final
    )
    
    print(
      signature_kegg_final[
        order(signature_kegg_final$p.adjust),
        c(
          "Description",
          "NES",
          "pvalue",
          "p.adjust"
        )
      ][
        1:min(15, nrow(signature_kegg_final)),
        ,
        drop = FALSE
      ]
    )
    
    
    ############################################################
    # 3. REMOVE SIGNATURE GENES
    ############################################################
    
    cat("\n===== SIGNATURE-GENE REMOVAL SENSITIVITY =====\n")
    
    signature_gsea_go_no3_final <-
      clusterProfiler::gseGO(
        geneList = gene_list_signature_no3,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 500,
        pAdjustMethod = "BH",
        eps = 0,
        nPermSimple = 100000,
        verbose = FALSE
      )
    
    signature_gsea_kegg_no3_final <-
      clusterProfiler::gseKEGG(
        geneList = gene_list_signature_no3,
        organism = "hsa",
        keyType = "ncbi-geneid",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps = 0,
        nPermSimple = 100000,
        verbose = FALSE
      )
    
    go_no3_final <- as.data.frame(
      signature_gsea_go_no3_final
    )
    
    kegg_no3_final <- as.data.frame(
      signature_gsea_kegg_no3_final
    )
    
    
    ############################################################
    # 4. SENSITIVITY CORRELATIONS
    ############################################################
    
    go_common_final <- intersect(
      signature_go_final$ID,
      go_no3_final$ID
    )
    
    kegg_common_final <- intersect(
      signature_kegg_final$ID,
      kegg_no3_final$ID
    )
    
    go_rho_final <- cor(
      signature_go_final$NES[
        match(
          go_common_final,
          signature_go_final$ID
        )
      ],
      go_no3_final$NES[
        match(
          go_common_final,
          go_no3_final$ID
        )
      ],
      method = "spearman"
    )
    
    kegg_rho_final <- cor(
      signature_kegg_final$NES[
        match(
          kegg_common_final,
          signature_kegg_final$ID
        )
      ],
      kegg_no3_final$NES[
        match(
          kegg_common_final,
          kegg_no3_final$ID
        )
      ],
      method = "spearman"
    )
    
    cat(
      "\nGO sensitivity Spearman rho = ",
      go_rho_final,
      "\n",
      sep = ""
    )
    
    cat(
      "KEGG sensitivity Spearman rho = ",
      kegg_rho_final,
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # 5. PREPARE FIGURE DATA
    ############################################################
    
    go_plot_df <- signature_go_final %>%
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust)
      ) %>%
      dplyr::arrange(
        dplyr::desc(abs(NES))
      ) %>%
      dplyr::slice_head(n = 12) %>%
      dplyr::mutate(
        Description = factor(
          Description,
          levels = rev(Description)
        ),
        neglog10FDR = -log10(p.adjust)
      )
    
    kegg_plot_df <- signature_kegg_final %>%
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust)
      ) %>%
      dplyr::arrange(
        dplyr::desc(abs(NES))
      ) %>%
      dplyr::slice_head(n = 12) %>%
      dplyr::mutate(
        Description = factor(
          Description,
          levels = rev(Description)
        ),
        neglog10FDR = -log10(p.adjust)
      )
    
    
    ############################################################
    # 6. GO PANEL
    ############################################################
    
    go_panel_final <- ggplot(
      go_plot_df,
      aes(
        x = NES,
        y = Description,
        size = neglog10FDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      ) +
      labs(
        title = "A  GO Biological Process",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        size = expression(-log[10] * "(FDR)")
      ) +
      theme_classic(
        base_size = 11
      )
    
    
    ############################################################
    # 7. KEGG PANEL
    ############################################################
    
    kegg_panel_final <- ggplot(
      kegg_plot_df,
      aes(
        x = NES,
        y = Description,
        size = neglog10FDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      ) +
      labs(
        title = "B  KEGG pathways",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        size = expression(-log[10] * "(FDR)")
      ) +
      theme_classic(
        base_size = 11
      )
    
    
    ############################################################
    # 8. FINAL FIGURE 4
    ############################################################
    
    figure4_final <-
      go_panel_final +
      kegg_panel_final +
      plot_annotation(
        title =
          "Pathway associations of the MYLK–PRKCB–PDGFC signature"
      )
    
    print(figure4_final)
    
    
    ############################################################
    # 9. SAVE FIGURE
    ############################################################
    
    ggsave(
      "Figure_4_signature_associated_GSEA_final.png",
      figure4_final,
      width = 12,
      height = 7.5,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      "Figure_4_signature_associated_GSEA_final.pdf",
      figure4_final,
      width = 12,
      height = 7.5,
      units = "in"
    )
    
    
    ############################################################
    # 10. SAVE FINAL GSEA RESULTS
    ############################################################
    
    write.csv(
      signature_go_final,
      "FINAL_signature_associated_GO_GSEA.csv",
      row.names = FALSE
    )
    
    write.csv(
      signature_kegg_final,
      "FINAL_signature_associated_KEGG_GSEA.csv",
      row.names = FALSE
    )
    
    write.csv(
      go_no3_final,
      "FINAL_signature_associated_GO_GSEA_no3.csv",
      row.names = FALSE
    )
    
    write.csv(
      kegg_no3_final,
      "FINAL_signature_associated_KEGG_GSEA_no3.csv",
      row.names = FALSE
    )
    
    
    ############################################################
    # 11. SAVE FINAL COMPLETE WORKSPACE
    ############################################################
    
    save(
      clinical,
      clinical_prad,
      tumor_samples_unique,
      vsd_tumor,
      gleason_analysis,
      cv_data_recovered,
      tcga_signature_clean,
      final_genes,
      final_weights,
      final_means,
      final_sds,
      gene_list_signature,
      gene_list_signature_no3,
      signature_gsea_go_final,
      signature_gsea_kegg_final,
      signature_gsea_go_no3_final,
      signature_gsea_kegg_no3_final,
      signature_go_final,
      signature_kegg_final,
      go_no3_final,
      kegg_no3_final,
      go_rho_final,
      kegg_rho_final,
      figure4_final,
      file =
        "TCGA_PRAD_FINAL_COMPLETE_WORKSPACE.RData"
    )
    
    cat(
      "\n===== FINAL GSEA + WORKSPACE COMPLETE =====\n"
    )
    ############################################################
    # FINAL FIGURE 4 — BALANCED GSEA DISPLAY
    # 6 negative + 6 positive pathways per panel
    ############################################################
    
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    
    ############################################################
    # GO: SELECT BOTH DIRECTIONS
    ############################################################
    
    go_negative <- signature_go_final %>%
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust),
        NES < 0
      ) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 6)
    
    go_positive <- signature_go_final %>%
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust),
        NES > 0
      ) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 6)
    
    go_plot_df <- dplyr::bind_rows(
      go_negative,
      go_positive
    ) %>%
      dplyr::distinct(ID, .keep_all = TRUE) %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(
        Description = factor(
          Description,
          levels = Description
        ),
        neglog10FDR = -log10(p.adjust)
      )
    
    ############################################################
    # KEGG: SELECT BOTH DIRECTIONS
    ############################################################
    
    kegg_negative <- signature_kegg_final %>%
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust),
        NES < 0
      ) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 6)
    
    kegg_positive <- signature_kegg_final %>%
      dplyr::filter(
        !is.na(NES),
        !is.na(p.adjust),
        NES > 0
      ) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 6)
    
    kegg_plot_df <- dplyr::bind_rows(
      kegg_negative,
      kegg_positive
    ) %>%
      dplyr::distinct(ID, .keep_all = TRUE) %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(
        Description = factor(
          Description,
          levels = Description
        ),
        neglog10FDR = -log10(p.adjust)
      )
    
    ############################################################
    # GO PANEL
    ############################################################
    
    go_panel_final <- ggplot(
      go_plot_df,
      aes(
        x = NES,
        y = Description,
        size = neglog10FDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      ) +
      labs(
        title = "A  GO Biological Process",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        size = expression(-log[10] * "(FDR)")
      ) +
      theme_classic(
        base_size = 11
      ) +
      theme(
        plot.title = element_text(face = "bold")
      )
    
    ############################################################
    # KEGG PANEL
    ############################################################
    
    kegg_panel_final <- ggplot(
      kegg_plot_df,
      aes(
        x = NES,
        y = Description,
        size = neglog10FDR
      )
    ) +
      geom_point() +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      ) +
      labs(
        title = "B  KEGG pathways",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        size = expression(-log[10] * "(FDR)")
      ) +
      theme_classic(
        base_size = 11
      ) +
      theme(
        plot.title = element_text(face = "bold")
      )
    
    ############################################################
    # COMBINE
    ############################################################
    
    figure4_final_balanced <-
      go_panel_final +
      kegg_panel_final +
      plot_annotation(
        title =
          "Pathway associations of the MYLK–PRKCB–PDGFC signature"
      )
    
    print(figure4_final_balanced)
    
    ############################################################
    # SAVE FINAL FIGURE
    ############################################################
    
    ggsave(
      "Figure_4_signature_associated_GSEA_final.png",
      figure4_final_balanced,
      width = 12,
      height = 7.5,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      "Figure_4_signature_associated_GSEA_final.pdf",
      figure4_final_balanced,
      width = 12,
      height = 7.5,
      units = "in"
    )
    
    ############################################################
    # SAVE EXACT PATHWAYS SHOWN IN FIGURE
    ############################################################
    
    write.csv(
      go_plot_df,
      "Figure_4_GO_pathways_shown.csv",
      row.names = FALSE
    )
    
    write.csv(
      kegg_plot_df,
      "Figure_4_KEGG_pathways_shown.csv",
      row.names = FALSE
    )
    
    cat("\n===== BALANCED FIGURE 4 COMPLETE =====\n")
    cat("GO pathways shown =", nrow(go_plot_df), "\n")
    cat("KEGG pathways shown =", nrow(kegg_plot_df), "\n")
    ############################################################
    # FINAL SAVE — FIGURE 4 + COMPLETE PROJECT WORKSPACE
    ############################################################
    
    save(
      clinical,
      clinical_prad,
      tumor_samples_unique,
      vsd_tumor,
      gleason_analysis,
      cv_data_recovered,
      tcga_signature_clean,
      final_genes,
      final_weights,
      final_means,
      final_sds,
      gene_list_signature,
      gene_list_signature_no3,
      signature_gsea_go_final,
      signature_gsea_kegg_final,
      signature_gsea_go_no3_final,
      signature_gsea_kegg_no3_final,
      signature_go_final,
      signature_kegg_final,
      go_no3_final,
      kegg_no3_final,
      go_rho_final,
      kegg_rho_final,
      go_plot_df,
      kegg_plot_df,
      figure4_final_balanced,
      file = "TCGA_PRAD_FINAL_COMPLETE_WORKSPACE.RData"
    )
    
    cat("\n===== FINAL WORKSPACE SAVED =====\n")
    
    print(
      file.exists(
        "TCGA_PRAD_FINAL_COMPLETE_WORKSPACE.RData"
      )
    )
    ############################################################
    # FINAL SUPPLEMENTARY TABLES
    #
    # Generates publication-ready CSV files from the locked
    # analysis objects.
    ############################################################
    
    library(dplyr)
    library(survival)
    
    dir.create(
      "Supplementary_Tables",
      showWarnings = FALSE
    )
    
    ############################################################
    # SUPPLEMENTARY TABLE 1
    # ORIGINAL 11-GENE TCGA SCREEN
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 1 =====\n")
    
    candidate_genes_11 <- c(
      "MYLK",
      "PRKCB",
      "PDGFD",
      "PRKCA",
      "PDGFC",
      "VEGFD",
      "MET",
      "MYLK4",
      "KDR",
      "EGFR",
      "PDGFRA"
    )
    
    screen_results <- data.frame(
      Gene = candidate_genes_11,
      Spearman_rho = NA_real_,
      P_value = NA_real_
    )
    
    for (i in seq_along(candidate_genes_11)) {
      
      g <- candidate_genes_11[i]
      
      test <- cor.test(
        gleason_analysis$gleason_score,
        gleason_analysis[[g]],
        method = "spearman",
        exact = FALSE
      )
      
      screen_results$Spearman_rho[i] <-
        unname(test$estimate)
      
      screen_results$P_value[i] <-
        test$p.value
    }
    
    screen_results$FDR <- p.adjust(
      screen_results$P_value,
      method = "BH"
    )
    
    screen_results <- screen_results %>%
      arrange(FDR)
    
    write.csv(
      screen_results,
      "Supplementary_Tables/Supplementary_Table_1_11_gene_screen.csv",
      row.names = FALSE
    )
    
    print(
      screen_results,
      row.names = FALSE
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 2
    # CLEAN 455-PATIENT AGE-ADJUSTED FIVE-GENE MODELS
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 2 =====\n")
    
    five_gene_models <- list()
    
    for (g in five_genes) {
      
      model <- lm(
        gleason_score ~
          get(g) +
          age_years,
        data = cv_data_recovered
      )
      
      co <- summary(model)$coefficients
      
      # Extract gene row safely
      gene_row <- co[
        rownames(co) == "get(g)",
        ,
        drop = FALSE
      ]
      
      # If get(g) doesn't survive as a coefficient name,
      # fit using a temporary variable.
      if (nrow(gene_row) == 0) {
        
        temp_data <- cv_data_recovered %>%
          dplyr::select(
            gleason_score,
            age_years,
            dplyr::all_of(g)
          )
        
        colnames(temp_data)[3] <-
          "gene_expression"
        
        model <- lm(
          gleason_score ~
            gene_expression +
            age_years,
          data = temp_data
        )
        
        co <- summary(model)$coefficients
        
        gene_row <- co[
          "gene_expression",
          ,
          drop = FALSE
        ]
      }
      
      five_gene_models[[g]] <- data.frame(
        Gene = g,
        Beta_Gleason = gene_row["Estimate"],
        SE = gene_row["Std. Error"],
        t_value = gene_row["t value"],
        P_value = gene_row["Pr(>|t|)"]
      )
    }
    
    five_gene_table <- bind_rows(
      five_gene_models
    ) %>%
      mutate(
        FDR = p.adjust(
          P_value,
          method = "BH"
        )
      ) %>%
      arrange(FDR)
    
    write.csv(
      five_gene_table,
      "Supplementary_Tables/Supplementary_Table_2_five_gene_models.csv",
      row.names = FALSE
    )
    
    print(
      five_gene_table,
      row.names = FALSE
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 3
    # FULL GO GSEA
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 3 =====\n")
    
    go_supp <- signature_go_final %>%
      arrange(p.adjust)
    
    write.csv(
      go_supp,
      "Supplementary_Tables/Supplementary_Table_3_full_GO_GSEA.csv",
      row.names = FALSE
    )
    
    cat(
      "GO pathways = ",
      nrow(go_supp),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 4
    # FULL KEGG GSEA
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 4 =====\n")
    
    kegg_supp <- signature_kegg_final %>%
      arrange(p.adjust)
    
    write.csv(
      kegg_supp,
      "Supplementary_Tables/Supplementary_Table_4_full_KEGG_GSEA.csv",
      row.names = FALSE
    )
    
    cat(
      "KEGG pathways = ",
      nrow(kegg_supp),
      "\n",
      sep = ""
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 5
    # GSEA SENSITIVITY — SIGNATURE GENES REMOVED
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 5 =====\n")
    
    go_sensitivity <- go_no3_final %>%
      arrange(p.adjust)
    
    kegg_sensitivity <- kegg_no3_final %>%
      arrange(p.adjust)
    
    write.csv(
      go_sensitivity,
      "Supplementary_Tables/Supplementary_Table_5_GO_sensitivity_no_signature_genes.csv",
      row.names = FALSE
    )
    
    write.csv(
      kegg_sensitivity,
      "Supplementary_Tables/Supplementary_Table_5_KEGG_sensitivity_no_signature_genes.csv",
      row.names = FALSE
    )
    
    ############################################################
    # SUPPLEMENTARY TABLE 6
    # THREE-GENE VS FIVE-GENE VS AGE CROSS-VALIDATION
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 6 =====\n")
    
    cv_final_supp <- read.csv(
      "final_3_vs_5_vs_age_CV_summary.csv"
    )
    
    cv_three_supp <- read.csv(
      "three_gene_CV_summary.csv"
    )
    
    write.csv(
      cv_final_supp,
      "Supplementary_Tables/Supplementary_Table_6_model_comparison.csv",
      row.names = FALSE
    )
    
    write.csv(
      cv_three_supp,
      "Supplementary_Tables/Supplementary_Table_6_all_three_gene_combinations.csv",
      row.names = FALSE
    )
    
    print(
      cv_final_supp,
      row.names = FALSE
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 7
    # FINAL TCGA T-STAGE MODEL
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 7 =====\n")
    
    tcga_T_recovered <- tcga_signature_clean %>%
      dplyr::select(
        patient_id,
        signature_clean,
        gleason_score,
        age_years
      ) %>%
      left_join(
        clinical_prad %>%
          dplyr::select(
            submitter_id,
            ajcc_pathologic_t
          ),
        by = c(
          "patient_id" =
            "submitter_id"
        )
      ) %>%
      mutate(
        T_advanced =
          case_when(
            grepl(
              "^T2",
              ajcc_pathologic_t
            ) ~ 0,
            grepl(
              "^T3|^T4",
              ajcc_pathologic_t
            ) ~ 1,
            TRUE ~ NA_real_
          )
      ) %>%
      filter(
        complete.cases(
          T_advanced,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    tcga_T_final_model <- glm(
      T_advanced ~
        signature_clean +
        gleason_score +
        age_years,
      data = tcga_T_recovered,
      family = binomial
    )
    
    tcga_T_coef <- summary(
      tcga_T_final_model
    )$coefficients
    
    tcga_T_table <- data.frame(
      Term = rownames(tcga_T_coef),
      Estimate = tcga_T_coef[, "Estimate"],
      SE = tcga_T_coef[, "Std. Error"],
      Z = tcga_T_coef[, "z value"],
      P_value = tcga_T_coef[, "Pr(>|z|)"]
    ) %>%
      mutate(
        OR = exp(Estimate),
        CI_lower = exp(
          Estimate -
            1.96 * SE
        ),
        CI_upper = exp(
          Estimate +
            1.96 * SE
        )
      )
    
    write.csv(
      tcga_T_table,
      "Supplementary_Tables/Supplementary_Table_7_TCGA_T_stage_model.csv",
      row.names = FALSE
    )
    
    print(
      tcga_T_table,
      row.names = FALSE
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 8
    # FINAL TCGA SURVIVAL MODELS
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 8 =====\n")
    
    # Reconstruct survival exactly according to the original script:
    # death -> days_to_death
    # alive -> days_to_last_follow_up
    
    survival_clinical <- clinical_prad %>%
      mutate(
        survival_time = ifelse(
          vital_status == "Dead",
          clinical$days_to_death[
            match(
              submitter_id,
              clinical$submitter_id
            )
          ],
          days_to_last_follow_up
        ),
        event = ifelse(
          vital_status == "Dead",
          1,
          0
        )
      )
    
    tcga_survival_clean <- survival_clinical %>%
      dplyr::select(
        submitter_id,
        survival_time,
        event
      ) %>%
      left_join(
        tcga_signature_clean %>%
          dplyr::select(
            patient_id,
            signature_clean,
            gleason_score,
            age_years
          ),
        by = c(
          "submitter_id" =
            "patient_id"
        )
      )
    
    # Unadjusted model
    cox_unadjusted <- survival::coxph(
      survival::Surv(
        survival_time,
        event
      ) ~ signature_clean,
      data = tcga_survival_clean
    )
    
    # Adjusted model
    tcga_survival_complete <- tcga_survival_clean %>%
      filter(
        complete.cases(
          survival_time,
          event,
          signature_clean,
          gleason_score,
          age_years
        )
      )
    
    cox_adjusted <- survival::coxph(
      survival::Surv(
        survival_time,
        event
      ) ~
        scale(signature_clean) +
        gleason_score +
        age_years,
      data = tcga_survival_complete
    )
    
    extract_cox <- function(
    model,
    model_name
    ) {
      
      s <- summary(model)
      
      out <- data.frame(
        Model = model_name,
        Term = rownames(s$coefficients),
        HR = s$coefficients[, "exp(coef)"],
        SE = s$coefficients[, "se(coef)"],
        Z = s$coefficients[, "z"],
        P_value = s$coefficients[, "Pr(>|z|)"],
        CI_lower = s$conf.int[, "lower .95"],
        CI_upper = s$conf.int[, "upper .95"]
      )
      
      out
    }
    
    survival_table <- bind_rows(
      extract_cox(
        cox_unadjusted,
        "Unadjusted"
      ),
      extract_cox(
        cox_adjusted,
        "Adjusted"
      )
    )
    
    write.csv(
      survival_table,
      "Supplementary_Tables/Supplementary_Table_8_survival_models.csv",
      row.names = FALSE
    )
    
    print(
      survival_table,
      row.names = FALSE
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 9
    # FINAL SIGNATURE DEFINITION
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 9 =====\n")
    
    signature_definition <- data.frame(
      Gene = final_genes,
      Weight = as.numeric(final_weights),
      TCGA_mean = as.numeric(final_means),
      TCGA_SD = as.numeric(final_sds)
    )
    
    write.csv(
      signature_definition,
      "Supplementary_Tables/Supplementary_Table_9_signature_definition.csv",
      row.names = FALSE
    )
    
    print(
      signature_definition,
      row.names = FALSE
    )
    
    
    ############################################################
    # SUPPLEMENTARY TABLE 10
    # EXTERNAL VALIDATION SUMMARY
    ############################################################
    
    cat("\n===== SUPPLEMENTARY TABLE 10 =====\n")
    
    external_validation <- data.frame(
      Cohort = c(
        "GSE70768",
        "TCGA-PRAD"
      ),
      N = c(
        111,
        455
      ),
      Signature = c(
        "MYLK–PRKCB–PDGFC",
        "MYLK–PRKCB–PDGFC"
      ),
      Endpoint = c(
        "Gleason score",
        "Gleason score"
      ),
      Odds_Ratio = c(
        0.2256582,
        0.7964547
      ),
      CI_lower = c(
        0.09547773,
        0.712795
      ),
      CI_upper = c(
        0.5333351,
        0.8899334
      ),
      P_value = c(
        0.0006927785,
        5.832804e-05
      )
    )
    
    write.csv(
      external_validation,
      "Supplementary_Tables/Supplementary_Table_10_external_validation.csv",
      row.names = FALSE
    )
    
    print(
      external_validation,
      row.names = FALSE
    )
    
    
    ############################################################
    # FINAL MANIFEST
    ############################################################
    
    cat(
      "\n====================================================\n"
    )
    
    cat(
      "SUPPLEMENTARY TABLE GENERATION COMPLETE\n"
    )
    
    cat(
      "====================================================\n"
    )
    
    print(
      list.files(
        "Supplementary_Tables"
      )
    )