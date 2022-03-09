##pre process/filter



library(ggplot2)

library(pROC)
library(vegan)
library(coin)
library(randomForest)
library(e1071)
library(ROCR)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

read.table("L6_abso.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)->L6_abso

read.table("fmt_reads_sum", header = T, sep="\t", stringsAsFactors = F)->reads_sum

read.table("fmt_meta.txt", header = T, stringsAsFactors = F, sep = "\t", na.strings = c("NA",""))->meta_config_av

read.table("ma_arare",  header = T, stringsAsFactors = F, sep = "\t") -> arare

meta_config <- meta_config_av



##equal between L6_abso and biom summarize-table results
setequal(colnames(L6_abso[,colSums(L6_abso) < 1000]), reads_sum[reads_sum$count < 1000, 'id'])
rm_samples <- colnames(L6_abso[,colSums(L6_abso) < 2500])

setdiff(colnames(L6_abso)[apply(L6_abso, 2, function(x){length(x[x>0]) < 20})], colnames(L6_abso[,colSums(L6_abso) < 2500]))

L6_abso[, !colnames(L6_abso) %in% rm_samples] -> L6_abso_fil


meta_fil_config <- meta_config[meta_config$SRA_Sample %in% colnames(L6_abso_fil) & meta_config$Donor_sra %in% colnames(L6_abso_fil) & meta_config$Previous_sra %in% colnames(L6_abso_fil),]

length(unique(meta_fil_config$ID, meta_fil_config$PRJ))

L6_abso_fil <- L6_abso_fil[, na.omit(unique(c(meta_fil_config$SRA_Sample, meta_fil_config$Donor_sra, meta_fil_config$Previous_sra)))]

#write.table(L6_abso_fil, "L6_absofil.txt", quote=F, sep="\t", row.names = T)

sample_size <- colSums(L6_abso_fil)
L6_rela_fil <- 100000 * t(apply(L6_abso_fil, 1, function(x){x/sample_size}))

head(colSums(L6_rela_fil))

##row means
data.frame(rownames(L6_rela_fil), rowMeans(L6_rela_fil))->L6_rela_fil_mean
colnames(L6_rela_fil_mean)<-c('genus', 'ab')


###filter species

len<-ncol(L6_abso_fil)
L6_abso_fil_sAg <- L6_abso_fil[apply(L6_abso_fil[, 1:len], 1, function(x){ length(x[x > 0]) >= 0.05 * len}),]

nrow(meta_fil_config)->rows 

sample_size <- colSums(L6_abso_fil_sAg)
len<-ncol(L6_rela_fil)
L6_rela_fil_sAg <- L6_rela_fil[apply(L6_rela_fil[, 1:len], 1, function(x){ length(x[x > 0]) >= 0.05 * len}),]


head(colSums(L6_rela_fil_sAg))

Others <- sapply(100000 - as.numeric(as.character(colSums(L6_rela_fil_sAg))), function(x){max(0, x)})
rbind(L6_rela_fil_sAg, Others)->L6_rela_fil_sAg_others

mean(colSums(L6_rela_fil_sAg_others))



A_simp_names <- function(line){
    #line <- as.character(line)
    tmp = line
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__PClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__[Clostridium]'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__P[Clostridium]'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__LClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__RClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Clostridium'){
        line  = 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__EClostridium'
    }
    line = strsplit(line, ";")
    len <- length(line[[1]])
    if(len == 1){
        return(line[[1]])
    }
    if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], 'g__')){
        len = len - 1
        while(len > 1){
            if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], '__')){
                len = len - 1
            }else{
                simp_name = paste(strsplit(line[[1]][len], '__')[[1]][2], '_', sep = '')
                
                break
            }
        }
    }else{
        simp_name = strsplit(line[[1]][len], '__')[[1]][2]
    }
    simp_name = gsub("\\W", ".", simp_name, perl=T)
    simp_name = paste('A_', simp_name, sep = "")
    simp_name
}

P_simp_names <- function(line){
    #line <- as.character(line)
    tmp = line
    tmp = line
    if(line == "k__Bacteria;p__;c__;o__;f__;g__"){return("P_Bacteria")}
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__PClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__[Clostridium]'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__P[Clostridium]'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__LClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__RClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Clostridium'){
        line  = 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__EClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Granulicatella'){
        line = 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Granulicatella[C]'
    }
    line = strsplit(line, ";")
    len <- length(line[[1]]) 
    if(len == 1){
        return(line[[1]])
    }
    if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], 'g__')){
        len = len - 1
        while(len > 1){
            if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], '__')){
                len = len - 1
            }else{
                simp_name = paste(strsplit(line[[1]][len], '__')[[1]][2], '_', sep = '')
                
                break
            }
        }
    }else{
        simp_name = strsplit(line[[1]][len], '__')[[1]][2]
    }
    simp_name = gsub("\\W", ".", simp_name, perl=T)
    simp_name = paste('B_', simp_name, sep = "")
    simp_name
}
D_simp_names <- function(line){
    #line <- as.character(line)
    tmp = line
    tmp = line
    if(line == "k__Bacteria;p__;c__;o__;f__;g__"){return("D_Bacteria")}
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__PClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__[Clostridium]'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__P[Clostridium]'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__LClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__RClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Clostridium'){
        line  = 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__EClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Granulicatella'){
        line = 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Granulicatella[C]'
    }
    line = strsplit(line, ";")
    len <- length(line[[1]]) 
    if(len == 1){
        return(line[[1]])
    }
    if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], 'g__')){
        len = len - 1
        while(len > 1){
            if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], '__')){
                len = len - 1
            }else{
                simp_name = paste(strsplit(line[[1]][len], '__')[[1]][2], '_', sep = '')
                
                break
            }
        }
    }else{
        simp_name = strsplit(line[[1]][len], '__')[[1]][2]
    }
    simp_name = gsub("\\W", ".", simp_name, perl=T)
    simp_name = paste('D_', simp_name, sep = "")
    simp_name
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
    dataframe->Matrix
    #bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
    ncol(Matrix) -> ncols
    bigones <- apply(Matrix,1, function(x){ length(x[x > 0]) > percent*ncols})
    Matrix_1 <- Matrix[bigones,]
    print(percent)
    return(Matrix_1)
}
simp_names <- function(line){
    #line <- as.character(line)
    if(line == 'Others'){return('Others')}
    tmp = line
    
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__PClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__[Clostridium]'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__P[Clostridium]'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__LClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Clostridium'){
        line = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__RClostridium'
    }
    if(line == 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Clostridium'){
        line  = 'k__Bacteria;p__Firmicutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae;g__EClostridium'
    }
    line = strsplit(line, ";")
    len <- length(line[[1]]) 
    if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], 'g__')){
        len = len - 1
        while(len > 1){
            if(endsWith(line[[1]][len], "Other") | endsWith(line[[1]][len], '__')){
                len = len - 1
            }else{
                simp_name = paste(strsplit(line[[1]][len], '__')[[1]][2], '_', sep = '')
                break
            }
        }
    }else{
        simp_name = strsplit(line[[1]][len], '__')[[1]][2]
    }
    
    simp_name
}



dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x<=0.0000001,pseudocount,x))
    
    for(i in 1:matrixColSize) {
        for(j in 1:matrixColSize) { 
            resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                                   as.vector(inMatrix[,j]))
        }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix) 
}



pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
    require(cluster)
    cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
    return(cluster)
}

pam.medoids=function(x,k) {
    require(cluster)
    medoids = as.vector(pam(as.dist(x), k, diss=TRUE)$medoids)
    return(medoids)
}




library(dplyr)
result <- meta_fil_config %>% 
    group_by(Previous_sra, PRJ) %>%
    filter(Time == max(Time)) %>%
    arrange(SRA_Sample, Previous_sra, PRJ, Time)
end_SRA <- as.character(result$SRA_Sample)
meta_fil_config_end <- meta_fil_config[meta_fil_config$SRA_Sample %in% end_SRA, ]


result2 <- meta_fil_config %>% 
    group_by(Previous_sra, PRJ) %>%
    filter(Time == min(Time)) %>%
    arrange(SRA_Sample, Previous_sra, PRJ, Time)
end_SRA2 <- as.character(result2$SRA_Sample)
meta_fil_config_end2 <- meta_fil_config[meta_fil_config$SRA_Sample %in% end_SRA2, ]



###entropy 
library(cluster)
library(clusterSim)
library(ade4)


###pam combined old
pre_data <- L6_rela_fil_sAg_others[,unique(c(meta_fil_config$Previous_sra))]/100000
pre_data_remove = noise.removal(pre_data, percent=0.01)
pre_data.dist=dist.JSD(pre_data_remove)

require(clusterSim)
pre_nclusters=NULL
pre_data.cluster=pam.clustering(pre_data.dist, k=2)
pre_data.medoids = pam.medoids(pre_data.dist, k=2)

pre_obs.pcoa=dudi.pco(pre_data.dist, scannf=F, nf=2)
pre_entro <- cbind(rownames(pre_obs.pcoa$li), paste(rep('before', length(pre_data.cluster)), pre_data.cluster, sep = ''))


don_data <- L6_rela_fil_sAg_others[,unique(c(meta_fil_config$Donor_sra))]/100000
head(colSums(don_data))
don_data_remove = noise.removal(don_data, percent=0.01)
don_data.dist=dist.JSD(don_data_remove)
don_nclusters=NULL
don_data.cluster=pam.clustering(don_data.dist, k=2)
don_data.medoids = pam.medoids(don_data.dist, k=2)

don_obs.pcoa=dudi.pco(don_data.dist, scannf=F, nf=2)
don_entro <- cbind(rownames(don_obs.pcoa$li), paste(rep('donor',length(don_data.cluster)), don_data.cluster,sep=''))


after_data <- L6_rela_fil_sAg_others[,unique(c(meta_fil_config$SRA_Sample))]/100000
head(colSums(after_data))
after_data_remove = noise.removal(after_data, percent=0.01)
after_data.dist=dist.JSD(after_data_remove)
after_nclusters=NULL
after_data.cluster=pam.clustering(after_data.dist, k=3)
after_obs.pcoa=dudi.pco(after_data.dist, scannf=F, nf=2)
after_entro <- cbind(rownames(after_obs.pcoa$li), paste(rep('after', length(after_data.cluster)), after_data.cluster, sep = ''))



tmp<-meta_fil_config
tmp <- merge(tmp, pre_entro, by.x='Previous_sra', by.y='V1')
tmp <- merge(tmp, don_entro, by.x='Donor_sra', by.y='V1')
tmp <- merge(tmp, after_entro, by.x='SRA_Sample', by.y='V1')
colnames(tmp) <- c(colnames(tmp)[1:16], 'pre_entro', 'don_entro', 'after_entro')
meta_fil_config_entro<-tmp

tmp<-meta_fil_config_end
tmp <- merge(tmp, pre_entro, by.x='Previous_sra', by.y='V1')
tmp <- merge(tmp, don_entro, by.x='Donor_sra', by.y='V1')
tmp <- merge(tmp, after_entro, by.x='SRA_Sample', by.y='V1')
colnames(tmp) <- c(colnames(tmp)[1:16], 'pre_entro', 'don_entro', 'after_entro')
meta_fil_config1<-tmp

tmp<-meta_fil_config_end2
tmp <- merge(tmp, pre_entro, by.x='Previous_sra', by.y='V1')
tmp <- merge(tmp, don_entro, by.x='Donor_sra', by.y='V1')
tmp <- merge(tmp, after_entro, by.x='SRA_Sample', by.y='V1')
colnames(tmp) <- c(colnames(tmp)[1:16], 'pre_entro', 'don_entro', 'after_entro')
meta_fil_config2<-tmp
# meta_fil_config2

meta_fil_config1_na <- meta_fil_config1[!meta_fil_config1$postfmt_symptoms %in% c(NA),]
meta_fil_config2_na <- meta_fil_config2[!meta_fil_config2$postfmt_symptoms %in% c(NA),]
L6_rela_fil_sAg_remove <- noise.removal(L6_rela_fil_sAg, percent=0.01)


L6_rela_fil_sAg_remove_simp <- L6_rela_fil_sAg_remove
L6_rela_fil_sAg_remove_simp_name <- sapply(as.character(rownames(L6_rela_fil_sAg_remove)), simp_names)
rownames(L6_rela_fil_sAg_remove_simp) <- c(L6_rela_fil_sAg_remove_simp_name)
