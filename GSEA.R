##### Load the required libraries
library(magrittr)
library(reshape2)
library(ggstatsplot)
library(ggsignif)
library(ggpubr)
library(GSEABase)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

rm(list = ls(all.names = TRUE))
gc()

##### My functions
human2mouse <- function(gl) {
  for (i in 1:length(gl)) {
    gl[i] <- paste(
      substr(gl[i],1,1),
      tolower(substr(gl[i],2,nchar(gl[i]))),
      sep = ""
    )
  } 
  gl
}

minmax_22 <- function(data) {
  temp_min <- apply(data, 2, min)
  temp_max <- apply(data, 2, max)
  temp_reduce <- temp_max-temp_min
  result_1 <- sweep(data, 2, temp_min, FUN="-")
  result_2 <- sweep(result_1, 2, temp_reduce, FUN="/")
  result_3 <- sweep(result_2, 2, 4, FUN="*")
  result_4 <- sweep(result_3, 2, 2, FUN="-")
}

##### Extract the CellChat input files from a Scanpy object
library(reticulate)
use_python("/Users/axis/miniconda3/envs/SCSseq/bin/python3.7", required = T)
py_config()

ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("./adata_processed.h5ad")
# access normalized data matrix
data.input <- t(py_to_r(ad_object$raw$X))
rownames(data.input) <- rownames(py_to_r(ad_object$raw$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
head(data.input[,1:5])
# access meta data
meta.data <- py_to_r(ad_object$obs)
meta <- meta.data

##### GSEA
# generate geneset
path.tp <- read.table("./pathways to plot.txt",sep = "\t",header = FALSE)
path.tp <- toupper(path.tp$V1) %>% 
  gsub("-","_",.) %>% 
  gsub(" ","_",.) %>% 
  paste("GOBP",.,sep = "_")
gene.set <- read.gmt("./c5.go.bp.v7.4.symbols.gmt")
gene.set <- gene.set[gene.set$term %in% path.tp,]
# test <- data.frame(unique(gene.set$term))
gene.set$gene <- human2mouse(gene.set$gene)

cell.list <- c('Mature CM', 'Fibroblast', 'Endothelial', 'Pericyte', 
  'Schwann Cell','DC-like Cell', 'Granulocyte', 'Macrophage (M1)', 
  'Macrophage (M2)', 'NK Cell', 'B Cell', 'Erythroid-like')
result <- data.frame(ID = path.tp)
for (cell in cell.list) {
  anno.cell <- meta[meta$cell_type == cell,1] # extract cell barcode
  cell.data <- data.input[,colnames(data.input) %in% anno.cell] # generate matrix of celltype to test
  oth.data <- data.input[,!colnames(data.input) %in% anno.cell] # generate matrix of other celltypes
  t.ls <- c() # t test
  p.ls <- c()
  for (i in 1:nrow(cell.data)) {
    t <- t.test(as.numeric(cell.data[i,]),as.numeric(oth.data[i,]),var.equal=TRUE)
    t.ls <- c(t.ls,t$statistic)
    p.ls <- c(p.ls,t$p.value)
  }
  # generate genelist
  gene.list <- data.frame(
    Symbol = rownames(cell.data),
    T.Value = t.ls,
    P.Value = p.ls
  )
  gene.list <- gene.list[gene.list$P.Value < 0.01,] # screen out deg
  tmp <- gene.list[!duplicated(gene.list$Symbol),] %>% 
    .[order(.$T.Value,decreasing = TRUE),] %>% 
    .[!is.na(.$Symbol),]
  gene.list <- as.numeric(tmp$T.Value) %>% 
    sort(.,decreasing = TRUE)
  names(gene.list) <- tmp$Symbol
  # gsea
  gsea.gobp.res <- GSEA(
    gene.list,
    TERM2GENE = gene.set,
    exponent = 1, 
    minGSSize = 1, 
    maxGSSize = 100,
    pvalueCutoff = 1, 
    pAdjustMethod = "BH",
    seed = TRUE
  )
  
  temp.NES <- gsea.gobp.res@result %>% 
    .[,c(1,5)]
  result <- merge(result,temp.NES,all = TRUE)
  names(result)[ncol(result)] <- cell
}

tp <- result[,-1]
rownames(tp) <- result$ID
library(pheatmap)
library(viridis)
# colorpalette <- colorRampPalette(rev(c("firebrick3","white","navy")))(5)
pheatmap(
  tp,
  color = inferno(10),
  cluster_col = FALSE,
  cluster_rows = FALSE,
  # scale = "row",
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10)

##### Pathway Evaluation
# generate geneset
path.tp <- read.table("./pathways to plot.txt",sep = "\t",header = FALSE)
path.tp <- toupper(path.tp$V1) %>% 
  gsub("-","_",.) %>% 
  gsub(" ","_",.) %>% 
  paste("GOBP",.,sep = "_")
gene.set <- read.gmt("./c5.go.bp.v7.4.symbols.gmt")
gene.set <- gene.set[gene.set$term %in% path.tp,]
# test <- data.frame(unique(gene.set$term))
gene.set$gene <- human2mouse(gene.set$gene)

cell.list <- c('Mature CM', 'Fibroblast', 'Endothelial', 'Pericyte', 
  'Schwann Cell','DC-like Cell', 'Granulocyte', 'Macrophage (M1)', 
  'Macrophage (M2)', 'NK Cell', 'B Cell', 'Erythroid-like')
result <- as.data.frame(array(dim = c(0,3)))
names(result) <- c("Pathway","sample_index","Activity")
for (cell in cell.list) {
  anno.cell <- meta[meta$cell_type == cell,1] # extract cell barcode
  cell.data <- data.input[,colnames(data.input) %in% anno.cell] # generate matrix of celltype to test
  temp.result <- as.data.frame(array(dim = c(0,ncol(cell.data))))
  names(temp.result) <- colnames(cell.data)
  for (path in path.tp) {
    genes.in.path <- gene.set[gene.set$term == path,2]
    temp <- cell.data[rownames(cell.data) %in% genes.in.path,]
    temp <- t(as.data.frame(colMeans(temp)))
    rownames(temp) <- path
    temp.result <- rbind(temp.result,temp)
  }
  temp.result$Pathway <- rownames(temp.result)
  temp.result <- melt(temp.result)
  names(temp.result)[c(2:3)] <- c("sample_index","Activity")
  result <- rbind(result,temp.result)
}
result <- merge(result,meta[,c(1,9)])

for (path in path.tp) {
  tp <- result[result$Pathway == path,]
  g <- ggbetweenstats(
    data = tp,
    x = cell_type,
    y = Activity,
    type = "parametric",
    pairwise.comparisons = FALSE,
    # pairwise.display = "sig",
    # p.adjust.method = "BH",
    # effsize.type = "unbiased",
    results.subtitle = FALSE,
    xlab = "Cell Type",
    ylab = "Pathway Activity",
    title = path,
    tr = 0,
    violin.args = list(width = 0.5,alpha = 0.2),
    # package = "ggsci",
    # palette = "nrc_npg"
  ) +
    theme_classic() +
    theme(
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = NA, color = "black",size = 1,linetype = 1),
      legend.position = "NA"
    )
  ggsave(g,file = paste0(path,'.pdf'),width = 13,height = 7,limitsize = FALSE)
}

# generate matrix for scanpy
library(stringr)
result <- data.frame(sample_index = meta[,1])
for (path in path.tp) {
  genes.in.path <- gene.set[gene.set$term == path,2]
  temp <- data.input[rownames(data.input) %in% genes.in.path,]
  temp <- colMeans(temp)
  
  result$NewPath <- temp # add to meta
  path <- gsub("GOBP_","",path) %>% 
    gsub("_"," ",.) %>% 
    str_to_title(.)
  names(result)[ncol(result)] <- path
}
result[,2:23] <- minmax_22(result[,2:23])
write.csv(result,"pathway_activity.csv",row.names = F)
