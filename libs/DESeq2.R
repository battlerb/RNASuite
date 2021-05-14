args <- commandArgs(trailingOnly = T)

count_file <- args[1]
sample_file <- args[2]
species <- args[3]
resoud <- args[4]
control_group <- args[5]
treat_group <- args[6]
pvalue <- as.numeric(args[7])
padjust <- as.numeric(args[8])
lfchange <- as.numeric(args[9])
plot_type <- args[10]


library(ggplot2)
library(pheatmap)
library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(pathview)
library(stringr)
library(topGO)
library(KEGG.db)
library(biomaRt)
if(species == 'mmu'){
  library(org.Mm.eg.db)
} else if(species == 'hsa'){
  library(org.Hs.eg.db)
}else if(species == 'ssc'){
  library(org.Ss.eg.db)
}




load_file <- function(count_file, sample_file, resoud, treat_group, control_group){
  countdata <- read.table(count_file, header = T, row.names = 1)
  sample_list <- read.table(sample_file, header = T, sep = ',')
  countdata <- countdata[,(-1:-5)]
  colnames(countdata) <- gsub('_sorted.bam','', colnames(countdata))
  for(i in 1:length(colnames(countdata))){
    colnames(countdata)[i] <- strsplit(colnames(countdata), '.', fixed = T)[[i]][length(strsplit(colnames(countdata), '.', fixed = T)[[i]])]
  }
  treat_num <- length(rownames(subset(sample_list, sample_list$Group==treat_group)))
  control_num <- length(rownames(subset(sample_list, sample_list$Group==control_group)))
  condition <- factor(c(rep(control_group,control_num), rep(treat_group, treat_num)))
  control_id <- subset(sample_list, sample_list$Group==control_group)$SampleID
  treat_id <- subset(sample_list, sample_list$Group==treat_group)$SampleID
  colData <- data.frame(row.names = c(as.character(control_id),as.character(treat_id)), condition)
  countdata <- countdata[,rownames(colData)]
  files <- list(countdata=countdata, conditions=condition, colData=colData)
  return(files)
}


run_deseq2 <- function(countdata, colData, condition, control_group, treat_group){
  dds <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c('condition', control_group, treat_group))
  res_dds <- list(dds, res)
  return(res_dds)
}


annotate_geneIDs <- function(res, species){
  res_table <- as.data.frame(res)
  ens_gene <- as.character(rownames(res_table))
  if(species == 'ssc'){
    for(i in 1:10){
      tryCatch({
        anno_gene<-getBM(attributes = c('ensembl_gene_id','entrezgene_id', 'external_gene_name'),
                         filters = 'ensembl_gene_id',
                         values = ens_gene,
                         mart = useMart('ensembl',dataset = 'sscrofa_gene_ensembl'),
                         useCache = F)
      })
      if(exists('anno_gene')){
        break
      }
    }
    row_num <- match(ens_gene,anno_gene$ensembl_gene_id)
    res_table$gene_symbol <- anno_gene[row_num,3]
    res_table$gene_entrezID <- anno_gene[row_num,2]
  } else {
    if(species == 'mmu'){
      db <- org.Mm.eg.db
    } else if(species == 'hsa'){
      db <- org.Hs.eg.db
    }
    res_table$gene_symbol <- mapIds(x = db,
                                    keys = ens_gene,
                                    keytype = 'ENSEMBL',
                                    column = 'SYMBOL',
                                    multiVals = 'first')
    res_table$gene_entrezID <- mapIds(x = db,
                                      keys = ens_gene,
                                      keytype = 'ENSEMBL',
                                      column = 'ENTREZID',
                                      multiVals = 'first')
  }
  return(res_table)
}


write_tables <- function(res_table, dds, resoud, pvalue, padjust, lfchange){
  sig_count_table <- subset(res_table, 
                            padj <= padjust & abs(log2FoldChange) > lfchange & pvalue <= pvalue)
  write.csv(sig_count_table, file = paste(resoud, 'sig_results.csv', sep = '/'), quote = F)
  write.csv(as.data.frame(counts(dds, normalize = T)), file = paste(resoud, 'normalized_count.csv', sep = '/'), quote = F)
  write.csv(res_table, file = paste(resoud, 'all_results.csv', sep = '/'),quote = F)
  return(sig_count_table)
}


draw_heatmap <- function(dds, sig_result, resoud, plot_type){
  ntd <- normTransform(dds)
  mat <- assay(ntd[row.names(sig_result)])
  heat <- pheatmap(mat = mat,
                   color = colorRampPalette(c("forestgreen", "black", "red"))(255),
                   scale = "row",                                        # Make fonts smaller
                   cellwidth = 600/length(colnames(mat)),                                         # Make the cells wider
                   show_colnames = T,
                   show_rownames = F,
                   border_color = NA,
                   main = 'Heatmap for significant differencial expression genes',
                   fontsize_col = 20,
                   fontsize = 30)
  
  if(plot_type=='png'|plot_type=='both'){
    ggsave(heat, filename = paste(resoud,"HeatMap.png",sep = '/'), width = 20, height = 20)
  } 
  if(plot_type=='pdf'|plot_type=='both'){
    ggsave(heat, filename = paste(resoud,"HeatMap.pdf",sep = '/'), width = 20, height = 20)
  }
}


draw_volcano <- function(results, resoud, plot_type, lfchange, padjust){
  data <- data.frame(gene = row.names(results),
                     pval = -log10(results$padj),
                     lfc  = results$log2FoldChange)
  padjust <- as.numeric(padjust)
  lfchange <- as.numeric(lfchange)
  # Remove any rows that have NA as an entry
  data <- na.omit(data)
  # Color the points which are up or down
  # If fold-change > 0 and pvalue > 1 (Increased significant)
  # If fold-change < 0 and pvalue > 1 (Decreased significant)
  logpadj <- -log10(padjust)
  data <- mutate(data, color = case_when(data$lfc > lfchange & data$pval > logpadj ~ "Increased",
                                         data$lfc < -1*lfchange & data$pval > logpadj ~ "Decreased",
                                         data$pval < logpadj | data$lfc > -1*lfchange & data$lfc < lfchange~ "nonsignificant"))
  
  # Make a basic ggplot2 object with x-y values
  # Add ggplot2 layers
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color)) +
    geom_point(size = .7, alpha = 0.8, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = "#CD4F39", Decreased = "#008B00", nonsignificant = "darkgray")) +
    theme_bw(base_size = 14) +                                                # change overall theme
    theme(legend.position = "right") +                                        # change the legend
    xlab(expression(log[2]("FoldChange"))) +                                  # Change X-Axis label
    ylab(expression(-log[10]("Padj"))) +                                    # Change Y-Axis label
    geom_vline(xintercept = -1*lfchange, color = "darkgrey", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") +
    geom_vline(xintercept = lfchange, color = "darkgrey", linetype = "dashed") +
    geom_hline(yintercept = logpadj, colour = "darkgrey", linetype = "dashed") +  # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") +                                     # Scale yaxis due to large p-values
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(label = "Volcano Plot") +
    theme(plot.title = element_text(hjust = 0.5))
  if(plot_type=='png'|plot_type=='both'){
    ggsave(vol, filename = paste(resoud,"Volcano.png",sep = '/'), width = 10, height = 10)
  } 
  if(plot_type=='pdf'|plot_type=='both'){
    ggsave(vol, filename = paste(resoud,"Volcano.pdf",sep = '/'), width = 10, height = 10)
  }
}

draw_pca <- function(dds, resoud, plot_type){
  ddsMat_rlog <- rlog(dds, blind = FALSE)
  pcaData <- plotPCA(ddsMat_rlog, intgroup=c("condition"), ntop = 500, returnData = T)
  pcaData <- pcaData[order(pcaData$condition,decreasing=F),]
  pcaData$condition <- factor(pcaData$condition)
  pca <- ggplot(pcaData, aes(x = PC1, y = PC2,
                                 color = condition)) +
    geom_point(size = 3, alpha = 1) +  
    geom_text(data = pcaData, aes(x = PC1, y=PC2, label = name), nudge_y = 2)+
    scale_x_continuous(limits = c(-35, 35)) +                # change limits to fix figure dimensions
    scale_y_continuous(limits = c(-35, 35)) +                # change limits to fix figure dimensions
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_blank(),                # 移除绘图区域的主要网格线
          panel.grid.minor = element_blank()) +              # 移除绘图区域的次要网格线
    ggtitle(label = "Principal Component Analysis (PCA)") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
  
  if(plot_type=='png'|plot_type=='both'){
    ggsave(pca, filename = paste(resoud,"PCA.png",sep = '/'), width = 10, height = 10)
  } 
  if(plot_type=='pdf'|plot_type=='both'){
    ggsave(pca, filename = paste(resoud,"PCA.pdf",sep = '/'), width = 10, height = 10)
  }
}


cluster_kegg <- function(enrichoud, sig_result, species, plot_type){
  gene_matrix <- sig_result$log2FoldChange
  names(gene_matrix) <- sig_result$gene_entrezID
  try({KEGG <- enrichKEGG(
    gene          = names(gene_matrix),
    organism   = species,
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.2,
    use_internal_data = T 
  )}, silent = T)

  if(species == 'mmu'){
    db <- org.Mm.eg.db
  } else if(species == 'hsa'){
    db <- org.Hs.eg.db
  }else if(species == 'ssc'){
    db <- org.Ss.eg.db
  }

  try({KEGG_csv <- setReadable(
    KEGG,
    OrgDb = db,
    keyType = "ENTREZID"
  )}, silent = T)

  try({
    KEGG_csv <- as.data.frame(KEGG_csv)
    KEGG_csv$enrich_factor <- 0
    for(i in 1:length(row.names(KEGG_csv))){
      KEGG_csv$enrich_factor[i]<-eval(parse(text=kegg_result$GeneRatio[i]))/eval(parse(text = kegg_result$BgRatio[i]))
      }
  }, silent = T)

  try({write.table(KEGG_csv,
                   file = paste(enrichoud, 'KEGG.tsv', sep = '/'), 
                   quote = F, 
                   row.names = F,
                   sep = '\t')},silent = T)

  try({
  KEGG_dot <- dotplot(KEGG, showCategory = 20, font.size=24) + ggtitle("Dotplot for KEGG")
  if(plot_type=='pdf'|plot_type=='both'){
    ggsave(KEGG_dot, filename = paste(enrichoud,"KEGG.pdf",sep = '/'), width = 24, height = 12)
  } 
  if(plot_type=='png'|plot_type=='both'){
    ggsave(KEGG_dot, filename = paste(enrichoud,"KEGG.png",sep = '/'), width = 24, height = 12)
  }
  },silent = T)
  try({
  kegglist <- as.data.frame(KEGG)
  if(length(kegglist$ID)>=20){
    kegglist <- kegglist$ID[1:20]
  }
  setwd(enrichoud)
  for (i in kegglist) {
    i <- pathview(gene.data  = gene_matrix,
                  pathway.id = i,
                  species    =species,
                  kegg.native = T,
                  same.layer = F)
  }},silent = T)
}

cluster_go <- function(enrichoud, sig_result, species, plot_type){
  gene_matrix <- sig_result$log2FoldChange
  names(gene_matrix) <- sig_result$gene_entrezID
  if(species == 'hsa'){
    db <- org.Hs.eg.db
  } else if(species == 'mmu'){
    db <- org.Mm.eg.db
  } else if(species == 'ssc'){
    db <- org.Ss.eg.db
  }
  tryCatch({go_cc <- enrichGO(gene          = names(gene_matrix),
                              OrgDb         = db,
                              ont           = "CC",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)})
  
  tryCatch({go_bp <- enrichGO(gene          = names(gene_matrix),
                              OrgDb         = db,
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)})
  
  tryCatch({go_mf <- enrichGO(gene          = names(gene_matrix),
                              OrgDb         = db,
                              ont           = "MF",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)})
  
  try({
    pdf(paste(enrichoud,"/GO_BP_graph.pdf",sep = '/'))
    bpgraph <- plotGOgraph(go_bp)
    dev.off()
  },silent =TRUE)
  try({
    pdf(paste(enrichoud,"/GO_CC_graph.pdf",sep = '/'))
    ccgraph <- plotGOgraph(go_cc)
    dev.off()
  },silent =TRUE)
  try({
    pdf(paste(enrichoud,"/GO_MF_graph.pdf",sep = '/'))
    mfgraph <- plotGOgraph(go_mf)
    dev.off()
  },silent =TRUE)
  try({
  bpdot <- dotplot(go_bp, title = "GO Biological Process", showCategory=15, font.size=20)
  bpdot<-bpdot+scale_y_discrete(labels = function(go_bp) str_wrap(go_bp,width = 50))
  if(plot_type=='png'|plot_type=='both'){
    ggsave(bpdot, filename = paste(enrichoud,"/GO_BP_DotPlot.png",sep = ''), width = 24, height = 12)
  } 
  if(plot_type=='pdf'|plot_type=='both'){
      ggsave(bpdot, filename = paste(enrichoud,"/GO_BP_DotPlot.pdf",sep = ''), width = 24, height = 12)
  }
    },silent = T)
  try({
  ccdot <- dotplot(go_cc, title = "GO Cellular Component", showCategory=15, font.size=20)
  ccdot<-ccdot+scale_y_discrete(labels = function(go_cc) str_wrap(go_cc,width = 50))
  if(plot_type=='png'|plot_type=='both'){
    ggsave(ccdot, filename = paste(enrichoud,"/GO_CC_DotPlot.png",sep = ''), width = 24, height = 12)
  } 
  if(plot_type=='pdf'|plot_type=='both'){
    ggsave(ccdot, filename = paste(enrichoud,"/GO_CC_DotPlot.pdf",sep = ''), width = 24, height = 12)
  }
    },silent = T)
  try({
  mfdot <- dotplot(go_mf, title = "GO Molecular Function", showCategory=15, font.size=20)
  mfdot<-mfdot+scale_y_discrete(labels = function(go_mf) str_wrap(go_mf,width = 50))
  if(plot_type=='png'|plot_type=='both'){
    ggsave(mfdot, filename = paste(enrichoud,"/GO_MF_DotPlot.png",sep = ''), width = 24, height = 12)
  } 
  if(plot_type=='pdf'|plot_type=='both'){
    ggsave(mfdot, filename = paste(enrichoud,"/GO_MF_DotPlot.pdf",sep = ''), width = 24, height = 12)
  }
    },silent = T)
  
  # Save result
  try({go_result_cc<-as.data.frame(go_cc)
  go_result_cc$enrich_factor <- 0
    for(i in 1:length(row.names(go_result_cc))){
      go_result_cc$enrich_factor[i]<-eval(parse(text=go_result_cc$GeneRatio[i]))/eval(parse(text = go_result_cc$BgRatio[i]))}}
    ,silent = T)
  try({go_result_bp<-as.data.frame(go_bp)
  go_result_bp$enrich_factor <- 0
    for(i in 1:length(row.names(go_result_bp))){
      go_result_bp$enrich_factor[i]<-eval(parse(text=go_result_bp$GeneRatio[i]))/eval(parse(text = go_result_bp$BgRatio[i]))}}
    ,silent = T)
  try({go_result_mf<-as.data.frame(go_mf)
  go_result_mf$enrich_factor <- 0
    for(i in 1:length(row.names(go_result_mf))){
      go_result_mf$enrich_factor[i]<-eval(parse(text=go_result_mf$GeneRatio[i]))/eval(parse(text = go_result_mf$BgRatio[i]))}}
    ,silent = T)
  
  try({write.csv(x = go_result_cc,
              file      = paste(enrichoud,"/go_cc.tsv",sep = '/'),
              quote     = F,
              row.names = F,
              sep = '\t')}, silent = T)
  try({write.csv(x = go_result_bp,
              file      = paste(enrichoud,"/go_bp.tsv",sep = '/'),
              quote     = F,
              row.names = F,
              sep = '\t')}, silent = T)
  try({write.table(x = go_result_mf,
              file      = paste(enrichoud,"/go_mf.tsv",sep = '/'),
              quote     = F,
              row.names = F,
              sep = '\t')}, silent = T)
}

files_data <- load_file(count_file = count_file,
                           sample_file = sample_file,
                           resoud = resoud,
                        treat_group = treat_group,
                        control_group = control_group)
res_dds <- run_deseq2(countdata = files_data$countdata, condition = files_data$conditions,
                  colData = files_data$colData, control_group = control_group,
                  treat_group = treat_group)
res_table <- annotate_geneIDs(res = res_dds[[2]], species = species)
sig_result <- write_tables(res_table = res_table, dds = res_dds[[1]], resoud = resoud, pvalue = pvalue, padjust = padjust, lfchange = lfchange)
draw_pca(dds = res_dds[[1]], resoud = resoud, plot_type = plot_type)
draw_heatmap(dds = res_dds[[1]], sig_result = sig_result, resoud = resoud, plot_type = plot_type)
draw_volcano(results = res_dds[[2]], resoud = resoud, plot_type = plot_type, lfchange = lfchange, padjust = padjust)
enrichoud <- paste(resoud, 'Enrichment', sep = '/')
cluster_go(enrichoud = enrichoud, sig_result = sig_result, species = species, plot_type = plot_type)
cluster_kegg(enrichoud = enrichoud, sig_result = sig_result, species = species, plot_type = plot_type)
