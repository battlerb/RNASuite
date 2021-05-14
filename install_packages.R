options(repos = structure(c(CRAN='http://mirrors.aliyun.com/CRAN/')))
options(BioC_mirror="https://mirrors.nju.edu.cn/bioconductor/")


install.packages('BiocManager')
library(BiocManager)

install(c('ggplot2', 'pheatmap', 'DESeq2', 'dplyr', 'clusterProfiler', 'pathview',
          'stringr', 'topGO', 'KEGG.db', 'org.Mm.eg.db', 'org.Hs.eg.db', 
          'org.Ss.eg.db', 'biomaRt', 'jsonlite'))