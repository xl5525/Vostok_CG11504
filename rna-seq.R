#DEG from DESeq2
library(DESeq2)
library(tidyverse)
library(ggplot2)
organism = "org.Dm.eg.db"
library(organism, character.only = TRUE)
library(ggrepel)


data=read_tsv("salmon.merged.gene_counts.tsv")
cts=data.frame(data[,c(1,3,4,5,6,7,8,9,10)])
cts=data.frame(cts[,-1],row.names = cts[,1])

cts=round(cts)
sample_list=colnames(cts)
sample_list=substr(sample_list,1,nchar(sample_list)-2)
ncdata=data.frame(sample_list,row.names = colnames(cts))
names(ncdata)=paste("condition")
dds <- DESeqDataSetFromMatrix(countData=cts, colData=ncdata, design= ~condition)
dds$condition <- relevel(dds$condition, ref = "yw")
result<-DESeq(dds)
res_cg11504_rna=data.frame(results(result,contrast=c("condition","cg11504","yw")))

# Keep only rows where the row names contain "FBgn"
res_cg11504_fbgn = res_cg11504_rna[grep("FBgn", row.names(res_cg11504_rna)), ]
res_cg11504_rna_geneSymbol=res_cg11504_fbgn
full_gene_list=bitr(row.names(res_cg11504_rna_geneSymbol), fromType = "FLYBASE", toType = "SYMBOL", OrgDb=organism)
# Make sure the gene ID columns are named consistently
colnames(full_gene_list)[colnames(full_gene_list) == "FLYBASE"] <- "gene_id"

# Merge to include gene symbols
res_cg11504_rna_geneSymbol <- merge(
  res_cg11504_rna_geneSymbol,
  full_gene_list[, c("gene_id", "SYMBOL")],  # Keep only relevant columns
  by.x = "row.names",  # Use the row names from res_cg11504_rna_geneSymbol
  by.y = "gene_id",
  all.x = TRUE  # Keep all entries in res_cg11504_rna_geneSymbol
)
res_cg11504_rna_geneSymbol=data.frame(res_cg11504_rna_geneSymbol[,-1],row.names =res_cg11504_rna_geneSymbol[,1])

#read changed loops
cg_p0.05_loops=read.delim("cg11504_p0.05.bedpe",header=FALSE)
cg_p0.05_genelist=c(cg_p0.05_loops$V8,cg_p0.05_loops$V9,cg_p0.05_loops$V10,cg_p0.05_loops$V11,cg_p0.05_loops$V12,cg_p0.05_loops$V13)
cg_p0.05_genelist=na.omit(unique(cg_p0.05_genelist))
cg_p0.05_genelist_fbgn <-bitr(cg_p0.05_genelist, fromType = "SYMBOL", toType = "FLYBASE", OrgDb=organism)

res_cg11504_rna_labeled_p0.05 <- merge(res_cg11504_fbgn, cg_p0.05_genelist_fbgn, by.x = 0, by.y = "FLYBASE", all.x = TRUE)
res_cg11504_rna_labeled_p0.05 = data.frame(res_cg11504_rna_labeled_p0.05[,-1],row.names = res_cg11504_rna_labeled_p0.05[,1])
res_cg11504_rna_labeled_p0.05$label <- ifelse(res_cg11504_rna_labeled_p0.05$padj < 0.05 & !is.na(res_cg11504_rna_labeled_p0.05$SYMBOL), res_cg11504_rna_labeled_p0.05$SYMBOL, NA)

res_cg11504_rna_labeled_p0.05_p50=res_cg11504_rna_labeled_p0.05[which(-log10(res_cg11504_rna_labeled_p0.05$padj)<50),]
res_cg11504_rna_labeled_p0.05_p50 <- res_cg11504_rna_labeled_p0.05_p50 %>%
  mutate(label_color = case_when(
    label != "" ~ "r",                # Color for non-empty labels
    padj < 0.05 ~ "b",                # Color for significant padj values
    TRUE ~ "g"                        # Default color for unchanged genes
  ))
res_cg11504_rna_labeled_p0.05_p50[which(res_cg11504_rna_labeled_p0.05_p50$label_color=="r"),]


# Create the plot
ggplot(res_cg11504_rna_labeled_p0.05_p50) +
  # Plot all points first
  geom_point(aes(
    x = log2FoldChange, 
    y = -log10(padj), 
    colour = factor(case_when(
      label_color == "b" ~ "Padj < 0.05",  # Significant points
      TRUE ~ "Unchanged genes"  # All other points
    ))
  ), size = 2, alpha = 0.3) +  
  
  # Plot red points associated with reduced loops on top
  geom_point(aes(
    x = log2FoldChange, 
    y = -log10(padj), 
    colour = "Associated with reduced loop"  # Red points
  ), data = subset(res_cg11504_rna_labeled_p0.05_p50, label_color == "r"), size = 2) +  # Filter for red points
  
  # Add text labels
  geom_text_repel(aes(
    x = log2FoldChange, 
    y = -log10(padj), 
    label = label
  ), size = 4, max.overlaps = Inf) +    
  
  # Axis labels and theme
  xlab("log2 fold change") + 
  ylab("-log10 padj") +
  labs(colour = "") +  # Set a meaningful legend title
  theme(
    legend.position = "top",  
    text = element_text(size = 15),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  
  # Customize color scale
  scale_color_manual(values = c(
    "Padj < 0.05" = "#30cbcf",  # Blue for significant points
    "Unchanged genes" = "gray",  # Gray for unchanged genes
    "Associated with reduced loop" = "red"  # Red for labeled points
  ))

#####################Diffbind_plot
res_cg11504_rna_labeled_p0.05_p50_diffbind = res_cg11504_rna_labeled_p0.05[which(-log10(res_cg11504_rna_labeled_p0.05$padj)<50),]
res_cg11504_rna_labeled_p0.05_p50_diffbind <- res_cg11504_rna_labeled_p0.05_p50_diffbind %>%
  mutate(label_color = case_when(
    (label %in% c("CCKLR-17D1", "CCKLR-17D3", "CG11741", "CG12535", "CG14691", "CG15578", "CG31373", "CG42458", "CG44999", "CG45263", "CG4704", "CG9010", "Cbp53E", "Сyp316a1", "Fas2", "Fas3", "Myc", "NK7.1", "NetA", "NetB", "Pvf3", "TkR86C", "dpr5", "kek1", "klg", "mamo")) ~ "r",               # Color for non-empty labels
    padj < 0.05 ~ "b",                # Color for significant padj values
    TRUE ~ "g"                        # Default color for unchanged genes
  ))

ggplot(res_cg11504_rna_labeled_p0.05_p50_diffbind) +
  # Plot all points first
  geom_point(aes(
    x = log2FoldChange, 
    y = -log10(padj), 
    colour = factor(case_when(
      label_color == "b" ~ "Padj < 0.05",  # Significant points
      TRUE ~ "Unchanged genes"  # All other points
    ))
  ), size = 2, alpha = 0.3) +  
  
  # Plot red points associated with reduced loops on top
  geom_point(aes(
    x = log2FoldChange, 
    y = -log10(padj), 
    colour = "Associated with reduced loop and Diffbind"  # Red points
  ), data = subset(res_cg11504_rna_labeled_p0.05_p50_diffbind, label_color == "r"), size = 2) +  # Filter for red points
  
  # Add text labels
  geom_text_repel(data = subset(res_cg11504_rna_labeled_p0.05_p50_diffbind, label_color == "r"), aes(
    x = log2FoldChange, 
    y = -log10(padj), 
    label = label
  ), size = 4, max.overlaps = Inf) +     
  guides(colour = guide_legend(nrow = 2)) +
  # Axis labels and theme
  xlab("log2 fold change") + 
  ylab("-log10 padj") +
  labs(colour = "") +  # Set a meaningful legend title
  theme(
    legend.position = "top",  
    text = element_text(size = 15),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  
  # Customize color scale
  scale_color_manual(values = c(
    "Padj < 0.05" = "#30cbcf",  # Blue for significant points
    "Unchanged genes" = "gray",  # Gray for unchanged genes
    "Associated with reduced loop and Diffbind" = "red"  # Red for labeled points
  ))



###########hbs and sns
# Filter for specific genes (hbs and sns)
genes_to_plot <- c("hbs", "sns")
subset_genes <- res_cg11504_rna_geneSymbol %>%
  filter(SYMBOL %in% genes_to_plot)

# Add a column for gene labels
subset_genes <- subset_genes %>%
  mutate(label = SYMBOL)

# Plot the data
ggplot(subset_genes, aes(x = SYMBOL, y = log2FoldChange, size = -log10(padj), color = padj < 0.05)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(
    title = "Dot Plot for hbs and sns",
    x = "Gene",
    y = "log2 Fold Change",
    size = "-log10(padj)",
    color = "Significant"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  )




