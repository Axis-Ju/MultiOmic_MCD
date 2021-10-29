# Transcriptome Data
mmu_heart_trs_path <- c("D:/R/Cardiodev Metabo/mmu_heart_trs/")
mmu_heart_trs_filename <- list.files(mmu_heart_trs_path)
mmu_heart_trs_readname <- c(paste("1W",1:3,sep = "-"),paste("2W",1:3,sep = "-"),
  paste("4W",1:3,sep = "-"),paste("8W",1:3,sep = "-"),
  paste("E10.5",1:3,sep = "-"),paste("E12.5",1:3,sep = "-"),
  paste("E14.5",1:3,sep = "-"),paste("E16.5",1:3,sep = "-"),
  paste("E18.5",1:3,sep = "-"),paste("P1",1:3,sep = "-"))
mmu_heart_trs <- data.frame(seq = seq(1,48820,1))
for (i in 1:length(mmu_heart_trs_filename)) {
  tmp_dir <- paste0(mmu_heart_trs_path,mmu_heart_trs_filename[i])
  tmp_read <- read.table(tmp_dir,header = TRUE, stringsAsFactors = FALSE)
  mmu_heart_trs <- data.frame(mmu_heart_trs,tmp_read[,c(1,7)])# fpkm
}
test <- data.frame()# test if ENSEMBLRIDs correspond
for (i in seq(2,(ncol(mmu_heart_trs)),4)) {
  tmp <- mmu_heart_trs[!(mmu_heart_trs[,i] == mmu_heart_trs[,i+2]),]
  test <- rbind(test,tmp)
}
mmu_heart_trs <- mmu_heart_trs[,-c(1,seq(4,(ncol(mmu_heart_trs)),2))]# assemble expression matrix
mmu_heart_trs <- data.frame(mmu_heart_trs[,1],mmu_heart_trs[,14:ncol(mmu_heart_trs)],mmu_heart_trs[,2:13])# rearrange follow time order
names(mmu_heart_trs) <- c("EnsembleID",paste("E10.5",1:3,sep = "-"),paste("E12.5",1:3,sep = "-"),
  paste("E14.5",1:3,sep = "-"),paste("E16.5",1:3,sep = "-"),paste("E18.5",1:3,sep = "-"),
  paste("P1",1:3,sep = "-"),paste("1W",1:3,sep = "-"),paste("2W",1:3,sep = "-"),
  paste("4W",1:3,sep = "-"),paste("8W",1:3,sep = "-"))# rename

# Proteome Data
mmu_heart_pro <- read.csv("D:/R/Cardiodev Metabo/pro_new.csv",header = TRUE,stringsAsFactors = FALSE)
tmp <- mmu_heart_pro[1,]# split merged gene names by ";"
for (i in 1:nrow(mmu_heart_pro)) {
  a <- data.frame(strsplit(mmu_heart_pro[i,1],split = ";"))
  b <- nrow(a)
  temp <- data.frame(mmu_heart_pro[i,2:ncol(mmu_heart_pro)])
  if (b > 1) {
    for (j in 1:(b-1)) {
      temp <- rbind(temp,mmu_heart_pro[i,2:ncol(mmu_heart_pro)])
    }
  }
  temp <- cbind(a,temp)
  names(temp) <- colnames(mmu_heart_pro)
  tmp <- rbind(tmp,temp)
}
mmu_heart_pro <- tmp[-1,];rm(tmp,i,a,b,temp,j)
mmu_heart_pro <- mmu_heart_pro[!duplicated(mmu_heart_pro[,1]),]

# Phospho-Proteome Data
mmu_heart_phos_raw <- read.csv("D:/R/Cardiodev Metabo/phospho.csv",header = TRUE,stringsAsFactors = FALSE)
mmu_heart_phos_raw <- mmu_heart_phos_raw[mmu_heart_phos_raw$Normalization.method == "divided by protein level",]
mmu_heart_phos <- data.frame(Name = mmu_heart_phos_raw[,3],
  Phosphosite = paste(mmu_heart_phos_raw[,8],
    mmu_heart_phos_raw[,10],sep = ""),
  mmu_heart_phos_raw[,13:52])

# Save Data
save(mmu_heart_trs, mmu_heart_pro, mmu_heart_phos, mmu_heart_phos_raw, file = "../SourceData/omic_data.Rds")