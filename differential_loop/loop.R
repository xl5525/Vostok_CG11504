library(DESeq2)
library(tidyverse)
library(purrr)


###A bigger loop list was used for background
lbe_looplist=read.delim("lbe_converted.bedpe",header=FALSE,row.names = 14, sep='\t')
lbw_looplist=read.delim("lbw_converted.bedpe",header=FALSE,row.names = 14, sep='\t')
elooplist=read.delim("eloops_from_MJ_3.2k.bedpe",header=FALSE,row.names = 7, sep='\t')
wlooplist=read.delim("wing_disc_loop_list.bedpe",header=FALSE,row.names = 7, sep='\t')
merged_looplist=rbind(loop_list,elooplist,wlooplist)
rows_to_remove <- c(row.names(lbe_looplist), row.names(lbw_looplist))
merged_looplist_dedup=merged_looplist[!row.names(merged_looplist) %in% rows_to_remove,]
merged_looplist_dedup$V7=row.names(merged_looplist_dedup)
write.table(merged_looplist_dedup,"merged_loop_list_dedup.bedpe",row.names = FALSE,col.names = FALSE,quote = FALSE, sep='\t')

#############################################
counts_3.2k <- list.files(pattern="*.csv", full.names = TRUE)

cts_list_all=list()
for (i in 1:length(counts_3.2k))
{
  cts_list_all[[i]] = data.frame(read.delim(counts_3.2k[i],sep='\t'))
}
cts <- cts_list_all %>% purrr::reduce(full_join, by = 'loop_id')

cts=data.frame(cts,row.names = 1)
cts=na.omit(cts)
cts=cts[,c(3,4,1,2,5,6)]

sample_list=colnames(cts)
sample_list=substr(sample_list,1,nchar(sample_list)-2)
ncdata=data.frame(sample_list,row.names = colnames(cts))
names(ncdata)=paste("condition")

dds <- DESeqDataSetFromMatrix(countData=cts, colData=ncdata, design= ~condition)

vsd = varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("condition"))

dds$condition <- relevel(dds$condition, ref = "lbrain")

result<-DESeq(dds)

res_cg11504 <- results(result,contrast=c("condition","cg11504","lbrain"))
res_poz <- results(result,contrast=c("condition","poz_lbrain","lbrain"))


#######################################Annotate looplist



closest_bed=read.delim("larval_brain_pt_0.0001_st_0.88_sz_2.0_numbered_annotated_cleaned.txt",header = FALSE)
# Remove complete duplicate rows from closest_bed
closest_bed_cleaned <- unique(closest_bed)

# Order the data frame by column 4 (the loop identifier)
closest_bed_cleaned <- closest_bed_cleaned[order(closest_bed_cleaned[, 4]), ]


# Extract the unique loop names without _A, _B, etc.
loop_names <- unique(gsub("_[A-Za-z]$", "", closest_bed_cleaned$V4))

# Initialize the annotated loop dataframe with loop names
# This creates a data frame with the same number of rows as loop_names
closest_annotated_loop <- data.frame(loop_name = loop_names, stringsAsFactors = FALSE)

# Loop through each row of the closest_annotated_loop data frame
for (i in 1:nrow(closest_annotated_loop)) {
  
  # Modify the matching pattern to ensure it matches exactly (i.e., lbloop1_A or lbloop1_B)
  matched_rows <- closest_bed_cleaned[grepl(paste0("^", closest_annotated_loop$loop_name[i], "_[AB]$"), closest_bed_cleaned$V4), ]
  
  # Extract gene annotations for _A and _B separately, assuming gene annotations are in column 8
  genes_A <- matched_rows[grepl("_A$", matched_rows$V4), 8]
  genes_B <- matched_rows[grepl("_B$", matched_rows$V4), 8]
  
  # Check if there are any genes in the _A category before assigning
  if (length(genes_A) > 0) {
    for (j in 1:length(genes_A)) {
      column_name <- paste0("gene_A_", j)
      # Create the column if it doesn't exist and assign the value
      closest_annotated_loop[i, column_name] <- genes_A[j]
    }
  } else {
    closest_annotated_loop[i, "gene_A_1"] <- NA  # Fill with NA if no _A gene found
  }
  
  # Check if there are any genes in the _B category before assigning
  if (length(genes_B) > 0) {
    for (k in 1:length(genes_B)) {
      column_name <- paste0("gene_B_", k)
      # Create the column if it doesn't exist and assign the value
      closest_annotated_loop[i, column_name] <- genes_B[k]
    }
  } else {
    closest_annotated_loop[i, "gene_B_1"] <- NA  # Fill with NA if no _B gene found
  }
}
row.names(closest_annotated_loop)=closest_annotated_loop$loop_name

res_cg11504_lbloop_annotated=merge(data.frame(res_cg11504),closest_annotated_loop[,c(2,3)],by=0)
res_cg11504_lbloop_annotated=data.frame(res_cg11504_lbloop_annotated[,-1],row.names = res_cg11504_lbloop_annotated$Row.names)

res_poz_lbloop_annotated=merge(data.frame(res_poz),closest_annotated_loop[,c(2,3)],by=0)
res_poz_lbloop_annotated=data.frame(res_poz_lbloop_annotated[,-1],row.names = res_poz_lbloop_annotated$Row.names)

lbloop_list=merge(loop_list,closest_annotated_loop,by=0)
lbloop_list=data.frame(lbloop_list[,-1],row.names = lbloop_list[,1])


##Meta-loops
file_list <- c("cg11504_1.csv", 
               "poz_lbrain_2.csv", 
               "poz_lbrain_1.csv", 
               "lbrain_2.csv", 
               "lbrain_1.csv", 
               "cg11504_2.csv")

file_list <- file.path(dir_path, file_list)

cts_meta <- read.csv(file_list[1],sep='\t')

# Loop through the rest of the files and merge them
for (i in 2:length(file_list)) {
  file_data <- read.csv(file_list[i],sep='\t')
  cts_meta <- left_join(cts_meta, file_data, by = "loop_id")
}

cts_meta=data.frame(cts_meta[,c(6,5,2,7,4,3)],row.names = cts_meta[,1])
colnames(cts_meta)=colnames(cts)
cts_all=rbind(cts,cts_meta)

sample_list_all=colnames(cts_all)
sample_list_all=substr(sample_list_all,1,nchar(sample_list_all)-2)
ncdata_all=data.frame(sample_list_all,row.names = colnames(cts_all))
names(ncdata_all)=paste("condition")

dds_all <- DESeqDataSetFromMatrix(countData=cts_all, colData=ncdata_all, design= ~condition)

vsd = varianceStabilizingTransformation(dds_all)
plotPCA(vsd, intgroup=c("condition"))

dds_all$condition <- relevel(dds_all$condition, ref = "lbrain")

result_all<-DESeq(dds_all)

res_cg11504_all <- results(result_all,contrast=c("condition","cg11504","lbrain"))
res_poz_all <- results(result_all,contrast=c("condition","poz_lbrain","lbrain"))

res_cg11504_all_ordered=data.frame(res_cg11504_all[order(res_cg11504_all$pvalue),])
res_poz_all_ordered=data.frame(res_poz_all[order(res_poz_all$pvalue),])

##Meta only
sample_list_meta=colnames(cts_meta)
sample_list_meta=substr(sample_list_meta,1,nchar(sample_list_meta)-2)
ncdata_meta=data.frame(sample_list_meta,row.names = colnames(cts_meta))
names(ncdata_meta)=paste("condition")

dds_meta <- DESeqDataSetFromMatrix(countData=cts_meta, colData=ncdata_meta, design= ~condition)

vsd = varianceStabilizingTransformation(dds_meta)
plotPCA(vsd, intgroup=c("condition"))

dds_meta$condition <- relevel(dds_meta$condition, ref = "lbrain")

result_meta<-DESeq(dds_meta)

res_cg11504_meta <- results(result_meta,contrast=c("condition","cg11504","lbrain"))
res_poz_meta <- results(result_meta,contrast=c("condition","poz_lbrain","lbrain"))

res_cg11504_meta_ordered=data.frame(res_cg11504_meta[order(res_cg11504_meta$pvalue),])
res_poz_meta_ordered=data.frame(res_poz_meta[order(res_poz_meta$pvalue),])

volcano_plot_differetial_meta=function(df,plot_title){
  df_th=""
  df_th=df$padj<0.05
  df$threshold=df_th
  df_ordered=data.frame(df[order(df$padj), ])
  df_ordered$genelabels=""
  df_ordered$genelabels[which(df_ordered$threshold=="TRUE")] <- rownames(df_ordered)[which(df_ordered$threshold=="TRUE")]
  options(ggrepel.max.overlaps = Inf)
  ggplot(df_ordered) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold),size=8) +
    geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels != "", genelabels, "")), size = 8) +    ggtitle(plot_title) +
    xlab("log2 fold change") + 
    ylab("-log10 padj") +
    theme(legend.position = "none", text = element_text(size = 20),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
}
volcano_plot_differetial_meta(res_poz_meta,"poz meta")
volcano_plot_differetial_meta(res_cg11504_meta,"cg11504 meta")



#####################################p0.05
volcano_plot_differetial_p0.05=function(df,name){
  df_ordered=data.frame(df[order(df$pvalue), ])
  df_ordered$label <- ifelse(1:nrow(df_ordered) <= 30 & df_ordered$gene_A_1 != "" & df_ordered$log2FoldChange<0, 
                             ifelse(df_ordered$gene_A_1 == df_ordered$gene_B_1, 
                                    df_ordered$gene_A_1, 
                                    ifelse(df_ordered$gene_B_1 != "", paste(df_ordered$gene_A_1, df_ordered$gene_B_1, sep = " - "), df_ordered$gene_A_1)), 
                             "")
  options(ggrepel.max.overlaps = Inf)
  ggplot(df_ordered) +
    geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = factor(pvalue < 0.05)),size=2) +
    geom_text_repel(aes(x = log2FoldChange, y = -log10(pvalue), label = label), size = 4) +    
    ggtitle(name) +
    xlab("log2 fold change") + 
    ylab("-log10 pvalue") +
    theme(legend.position = "none", text = element_text(size = 15),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
}
volcano_plot_differetial_p0.05(res_poz_lbloop_annotated,"POZ")
volcano_plot_differetial_p0.05(res_cg11504_lbloop_annotated,"cg11504")

cg11504_sig_loop_list 

# First, find the rows in `res_cg11504` where the padj value is < 0.05
significant_rows <- rownames(res_cg11504)[res_cg11504$padj < 0.05]

# Then, filter these significant rows based on whether they contain "lbloop" in `df`
filtered_rows <- significant_rows[grepl("lbloop", significant_rows)]

# Finally, use the filtered rows to subset `loop_list`
cg11504_sig_loop_list <- loop_list[filtered_rows, ]
cg11504_sig_loop_list$V7=row.names(cg11504_sig_loop_list)
write.table(cg11504_sig_loop_list,"cg11504_sig_loop_list.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# First, find the rows in `res_poz` where the padj value is < 0.05
significant_rows_poz <- rownames(res_poz)[res_poz$padj < 0.05]

# Then, filter these significant rows based on whether they contain "lbloop" in `df`
filtered_rows_poz <- significant_rows_poz[grepl("lbloop", significant_rows_poz)]

# Finally, use the filtered rows to subset `loop_list`
poz_sig_loop_list <- loop_list[filtered_rows_poz, ]
poz_sig_loop_list$V7=row.names(poz_sig_loop_list)
write.table(poz_sig_loop_list,"poz_sig_loop_list.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


cg11504_sig_loop_list_p0.05_rows=rownames(res_cg11504)[which(res_cg11504$pvalue<0.05&res_cg11504$log2FoldChange<0)]
filtered_rows_cg_p0.05 <- cg11504_sig_loop_list_p0.05_rows[grepl("lbloop", cg11504_sig_loop_list_p0.05_rows)]
cg11504_sig_loop_list_p0.05 <- loop_list[filtered_rows_cg_p0.05, ]
cg11504_sig_loop_list_p0.05$V7=row.names(cg11504_sig_loop_list_p0.05)
write.table(cg11504_sig_loop_list_p0.05,"cg11504_sig_loop_list_p0.05.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

poz_sig_loop_list_p0.05_rows=rownames(res_poz)[which(res_poz$pvalue<0.05&res_poz$log2FoldChange<0)]
filtered_rows_poz_p0.05 <- poz_sig_loop_list_p0.05_rows[grepl("lbloop", poz_sig_loop_list_p0.05_rows)]
poz_sig_loop_list_p0.05 <- loop_list[filtered_rows_poz_p0.05, ]
poz_sig_loop_list_p0.05$V7=row.names(poz_sig_loop_list_p0.05)
write.table(poz_sig_loop_list_p0.05,"poz_sig_loop_list_p0.05.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Extract row names where padj >= 0.05 for both datasets
unchanged_poz <- rownames(res_poz)[which(res_poz$padj >= 0.05)]
unchanged_cg11504 <- rownames(res_cg11504)[which(res_cg11504$padj >= 0.05)]

# Find the intersection of the two sets of row names
unchanged_anchors <- intersect(unchanged_poz, unchanged_cg11504)

unchanged_anchor_list = loop_list[unchanged_anchors,]
write.table(unchanged_anchor_list,"unchanged_both_0.05.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

overlapped_p0.05=row.names(poz_sig_loop_list_p0.05)[row.names(poz_sig_loop_list_p0.05) %in% row.names(cg11504_sig_loop_list_p0.05)]
overlapped_p0.05_list=loop_list[overlapped_p0.05,]
overlapped_p0.05_list$V7=row.names(overlapped_p0.05_list)
write.table(overlapped_p0.05_list,"overlapped_p0.05_list.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

unchanged_poz_p0.05 <- rownames(res_poz)[which(res_poz$pvalue >= 0.05)]
unchanged_cg11504_p0.05 <- rownames(res_cg11504)[which(res_cg11504$pvalue >= 0.05)]

# Find the intersection of the two sets of row names
unchanged_anchors_p0.05 <- intersect(unchanged_poz_p0.05, unchanged_cg11504_p0.05)

unchanged_anchor_list_p0.05 = loop_list[unchanged_anchors_p0.05,]
write.table(unchanged_anchor_list,"unchanged_both_p0.05.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

annotated_cg11504_p0.05=annotated_looplist[filtered_rows_cg_p0.05,]
write.table(annotated_cg11504_p0.05,"annotated_cg11504_p0.05.bedpe",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


annotated_poz_p0.05=annotated_looplist[filtered_rows_poz_p0.05,]
