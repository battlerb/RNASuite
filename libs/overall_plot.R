library(ggplot2)
library(pheatmap)
library(DESeq2)


args <- commandArgs(trailingOnly = T)
count_file <- args[1]
sample_file <- args[2]
output <- args[3]


load_file <- function(count_file, output){
  countdata <- read.table(count_file, header = T, row.names = 1)
  countdata <- countdata[,(-1:-5)]
  colnames(countdata) <- gsub('_sorted.bam','', colnames(countdata))
  for(i in 1:length(colnames(countdata))){
    colnames(countdata)[i] <- strsplit(colnames(countdata), '.', fixed = T)[[i]][length(strsplit(colnames(countdata), '.', fixed = T)[[i]])]
  }
  return(countdata)
}


plot_pca <- function(countdata, sample_list, output){
  sample_list <- read.table(sample_file, header = T, sep = ',')
  rownames(sample_list) <- sample_list[,1]
  condition <- factor(sample_list$Group)
  colData <- data.frame(row.names = rownames(sample_list), condition)
  colData <- data.frame(colData[rownames(colData)[match(rownames(colData), colnames(countdata))],],
                        row.names = rownames(colData)[as.double(match(colnames(countdata), rownames(colData)))])
  colnames(colData) <- 'condition'
  dds <- DESeqDataSetFromMatrix(countdata, colData, ~ condition)
  rld <- rlog(dds)
  PCA <- plotPCA(rld, intgroup=c("condition"), ntop = 500, returnData = T)
  PCA <- PCA[order(PCA$condition,decreasing=F),]
  PCA$condition <- factor(PCA$condition)
  pca <- ggplot(PCA, aes(x = PC1, y = PC2,
                         color = condition)) +
    geom_point(size = 3, alpha = 0.5) +  
    geom_text(data = PCA, aes(x = PC1, y=PC2, label = name),vjust = -1)+
    scale_x_continuous(limits = c(-50, 50)) +                # change limits to fix figure dimensions
    scale_y_continuous(limits = c(-50, 50)) +                # change limits to fix figure dimensions
    theme_bw(base_size = 20) +
    theme(panel.grid.major = element_blank(),                # 移除绘图区域的主要网格线
          panel.grid.minor = element_blank()) +              # 移除绘图区域的次要网格线
    ggtitle(label = "Principal Component Analysis (PCA)") +
    theme(plot.title = element_text(vjust = -1),
          legend.title = element_blank())
  pdf(paste(output,"Overall_PCA.pdf",sep = '/'), paper = 'a4')
  print(pca)
  dev.off()
  png(paste(output,"Overall_PCA.png",sep = '/'), width = 1366, height = 768)
  print(pca)
  dev.off()
 }


plot_heatmap <- function(countdata, output){
  mat <- countdata[rowMeans(countdata)>10,][0:500,]
  mat <- as.data.frame(scale(mat))
  colnames(mat)
  colnames(countdata)
  pdf(paste(output,"Overall_HeatMap.pdf",sep = '/'),paper = 'a4')
  pheatmap(mat = mat,
                   color = colorRampPalette(c("forestgreen", "black", "red"))(255),
                   scale = "row",                                        # Make fonts smaller
                   cellwidth = 200/length(colnames(countdata)),
                   show_rownames = F,
                   border_color = NA,
                   main = 'Heatmap for all samples',
                   fontsize_col = 10,
                   fontsize = 15,
                   show_colnames = T)
  dev.off()
  png(paste(output,"Overall_HeatMap.png",sep = '/'), width = 1366, height = 768)
  pheatmap(mat = mat,
           color = colorRampPalette(c("forestgreen", "black", "red"))(255),
           scale = "row",                                        # Make fonts smaller
           cellwidth = 600/length(colnames(countdata)),
           show_rownames = F,
           border_color = NA,
           main = 'Heatmap for all samples',
           fontsize_col = 10,
           fontsize = 15,
           show_colnames = TRUE)
  dev.off()
}
  
  
  
  #ggsave(heat, filename = paste(output,"Overall_HeatMap.png",sep = '/'), width = 20, height = 20)
  #ggsave(heat, filename = paste(output,"Overall_HeatMap.pdf",sep = '/'), width = 20, height = 20)


countdata <- load_file(count_file = count_file,
                       output = output)
plot_pca(countdata, sample_file, output)
plot_heatmap(countdata, output)