##pre in FMT model selection
source('pre_processing.R')

###Preprocessing
L6_rela_fil_sAg_remove <- noise.removal(L6_rela_fil_sAg, percent=0.05)

meta_fil_config1_na <- meta_fil_config1[!meta_fil_config1$postfmt_symptoms %in% c(NA),]

before_L6 <- L6_rela_fil_sAg_remove[, meta_fil_config1_na[, 'Previous_sra']]
colnames(before_L6) <- meta_fil_config1_na$SRA_Sample

before_L6 <- as.data.frame(before_L6, stringsAsFactors = F)
rownames(before_L6) <- sapply(rownames(before_L6), P_simp_names)

donor_L6 <- (L6_rela_fil_sAg_remove[, meta_fil_config1_na[, 'Donor_sra']])

colnames(donor_L6) <- meta_fil_config1_na$SRA_Sample

donor_L6 <- as.data.frame(donor_L6, stringsAsFactors = F)
rownames(donor_L6) <- sapply(rownames(donor_L6), D_simp_names)

before_donor_L6 <- as.data.frame(cbind(t(donor_L6), t(before_L6)), stringsAsFactors = F)


##add distance
meta_config_in_this<-meta_fil_config1_na
DP_distance <- NULL
len <- nrow(meta_config_in_this)
L6_rela_fil_sAg_others_remove <- noise.removal(L6_rela_fil_sAg_others, percent=0.01)
for(i in 1:len){
    x1 <- as.numeric(as.character(L6_rela_fil_sAg_others_remove[, meta_config_in_this[i, "Donor_sra"]]))
    x2 <- as.numeric(as.character(L6_rela_fil_sAg_others_remove[, meta_config_in_this[i, "Previous_sra"]]))
    DP_distance <- rbind(DP_distance, c(meta_config_in_this[i, "SRA_Sample"], vegdist(rbind(x1, x2), method='bray')))
}

DP_distance_dat <- as.data.frame(DP_distance, stringsAsFactors = F)
colnames(DP_distance_dat) <- c('SRA_Sample', 'distance')
rownames(DP_distance_dat) <- DP_distance_dat$SRA_Sample

##add alpha diversity
tmp_arare <- arare
rownames(tmp_arare) <- tmp_arare$X
# tmp_arare <- tmp_arare[rownames(before_donor_L6), c('X', 'shannon')]
cbind_arare <- cbind(tmp_arare[meta_config_in_this$Donor_sra, c('X', 'shannon')], tmp_arare[meta_config_in_this$Previous_sra, c('X', 'shannon')])
if(unique(cbind_arare[,1] == meta_config_in_this$Donor_sra & cbind_arare[,3] == meta_config_in_this$Previous_sra)==c(TRUE)){
    cn_arare = cbind(meta_config_in_this$SRA_Sample, cbind_arare, meta_config_in_this$postfmt_symptoms, meta_config_in_this$pre_entro, meta_config_in_this$don_entro, meta_config_in_this$PRJ)
    colnames(cn_arare) <- c('names', 'Dnames', 'Darare', 'Pnames', 'Parare', 'symptom', 'pre_entro', 'don_entro', 'PRJ')
}

# tmp <- c_arare[DP_distance_dat$SRA_Sample]
if(unique(DP_distance_dat$SRA_Sample == (cn_arare$names))){
    DP_distance_dat_arare <- cbind(DP_distance_dat, cn_arare)
    # colnames(DP_distance_dat) <- c('SRA_Sample', 'distance', 'arare')
}

DP_distance_dat_arare$distance <- as.numeric(DP_distance_dat_arare$distance)


##combining
# before_donor_L6$B_alpha <- tmp_arare[meta_config_in_this$Previous_sra, c('shannon')]
before_donor_L6$distance <- DP_distance_dat_arare$distance * 100
before_donor_L6$D_alpha <- DP_distance_dat_arare$Darare#tmp_arare[meta_config_in_this$Donor_sra, c('shannon')]

before_donor_L6$DB_alpha <- (DP_distance_dat_arare$Darare + DP_distance_dat_arare$Parare)
# before_donor_L6$D_alpha <- as.numeric((c_arare))


##add enterotype
feature_abun_dat_all <- before_donor_L6
rownames(feature_abun_dat_all) == meta_fil_config1_na$SRA_Sample
feature_abun_dat_all$Disease <- as.numeric(as.factor(as.character(ifelse(meta_fil_config1_na$Diease1 == "CDI", "CDI", "IBD"))))

feature_abun_dat_all$don_entro <- as.numeric(as.factor(as.character(meta_fil_config1_na$don_entro)))*1
feature_abun_dat_all$pre_entro <- as.factor(as.character(meta_fil_config1_na$pre_entro))

# feature_abun_dat$same_entro <- ifelse(feature_abun_dat$don_entro == feature_abun_dat$pre_entro)
feature_abun_dat_all$y <- as.factor(as.character(meta_fil_config1_na$postfmt_symptoms))



##ibd:2
feature_abun_dat <- feature_abun_dat_all[,]
feature_abun_dat_ibd <- feature_abun_dat[feature_abun_dat$Disease %in% c('2'),]
feature_abun_dat_cdi <- feature_abun_dat[!feature_abun_dat$Disease %in% c('2'),]


###validation data
cluster_medoids.JSD <- function(inMatrix, pseudocount=0.000001, k=2, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    # resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    cluster <- rep(0, matrixColSize)
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x<=0.0000001,pseudocount,x))
    
    for(i in 1:(matrixColSize)){
        cluster[i] = ifelse(JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,matrixColSize])) > 
                                JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,matrixColSize-1])), 1, 2)
    }
    return(cluster)
}


##Preprocessing
L6_abso_fil <- noise.removal(L6_abso, 0.01)
val_meta_config1 <- val_meta_config[val_meta_config$Donor_sra %in% colnames(L6_abso), ]

library(dplyr)
result <- val_meta_config1 %>% 
    group_by(Previous_sra, PRJ) %>%
    filter(Time == max(Time)) %>%
    arrange(SRA_Sample, Previous_sra, PRJ, Time)
end_SRA <- as.character(result$SRA_Sample)
val_meta_config2 <- val_meta_config1[val_meta_config1$SRA_Sample %in% end_SRA, ]


val_L6_before <- L6_abso_fil[,val_meta_config2$Previous_sra]

sample_size <- colSums(val_L6_before)
val_L6_before <- 100000 * t(apply(val_L6_before, 1, function(x){x/sample_size}))

colnames(val_L6_before) <- val_meta_config2$SRA_Sample
val_L6_before <- as.data.frame(val_L6_before, stringsAsFactors = F)
rownames(val_L6_before) <- sapply(rownames(val_L6_before), P_simp_names)

val_L6_donor <- L6_abso_fil[,val_meta_config2$Donor_sra]

sample_size <- colSums(val_L6_donor)
val_L6_donor <- 100000 * t(apply(val_L6_donor, 1, function(x){x/sample_size}))

colnames(val_L6_donor) <- val_meta_config2$SRA_Sample
val_L6_donor <- as.data.frame(val_L6_donor, stringsAsFactors = F)
rownames(val_L6_donor) <- sapply(rownames(val_L6_donor), D_simp_names)

val_before_donor_L6 <- as.data.frame(cbind(t(val_L6_donor), t(val_L6_before)), stringsAsFactors = F)


val_pre_data <- L6_abso_fil[,unique(c(val_meta_config2$Previous_sra, pre_data.medoids))]
sample_size <- colSums(val_pre_data)
val_pre_data <- 1 * t(apply(val_pre_data, 1, function(x){x/sample_size}))
val_pre_data_remove = noise.removal(val_pre_data, percent=0.01)

val_pre_data.cluster= cluster_medoids.JSD(val_pre_data_remove, 0.000001, 2)

samples = dim(val_pre_data_remove)[2]
val_meta_config2$pre_entro <- paste(rep('before', length(1:(samples-2))), val_pre_data.cluster[1:(samples-2)], sep = '')

# pre_medoids_dist = val_pre_data.dist[,c(samples-1, samples)]
# require(clusterSim)
# # pre_nclusters=NULL
# val_pre_data.cluster=as.vector(pam(as.dist(val_pre_data.dist), 2, medoids = c(samples, samples-1), diss=TRUE)$clustering)

val_don_data <- L6_abso_fil[,(c(val_meta_config2$Donor_sra, don_data.medoids))]
sample_size <- colSums(val_don_data)
val_don_data <- 1 * t(apply(val_don_data, 1, function(x){x/sample_size}))
val_don_data_remove = noise.removal(val_don_data, percent=0.01)

val_don_data.cluster= cluster_medoids.JSD(val_don_data_remove, 0.000001, 2)
samples = dim(val_don_data_remove)[2]
val_meta_config2$don_entro <- val_don_data.cluster[1:(samples-2)]


##construct dataframe
val_before_donor_L6$Disease <- as.numeric(as.factor(as.character(ifelse(val_meta_config2$Diease1 == "CDI", "CDI", "IBD"))))
val_before_donor_L6$don_entro <- as.numeric(as.factor(as.character(val_meta_config2$don_entro)))*1
# val_before_donor_L6$don_entro <- as.numeric(as.factor(val_meta_config2$don_entro))
pre_entro <- cbind(rownames(pre_obs.pcoa$li), paste(rep('before', length(pre_data.cluster)), pre_data.cluster, sep = ''))


meta_config_in_this<-val_meta_config2

### colnames(meta_config_in_this[])<-c('Donor_sra')
DP_distance <- NULL
len <- nrow(meta_config_in_this)
# val_before_donor_L6 <- noise.removal(val_before_donor_L6, percent=0.01)
for(i in 1:len){
    x1 <- as.numeric(as.character(val_L6_donor[, meta_config_in_this[i, "SRA_Sample"]]))
    x2 <- as.numeric(as.character(val_L6_before[, meta_config_in_this[i, "SRA_Sample"]]))
    DP_distance <- rbind(DP_distance, c(meta_config_in_this[i, "SRA_Sample"], vegdist(rbind(x1, x2), method='bray')))
}

DP_distance_dat <- as.data.frame(DP_distance, stringsAsFactors = F)
colnames(DP_distance_dat) <- c('SRA_Sample', 'distance')
rownames(DP_distance_dat) <- DP_distance_dat$SRA_Sample


tmp_arare <- arare
rownames(tmp_arare) <- tmp_arare$X
# tmp_arare <- tmp_arare[rownames(before_donor_L6), c('X', 'shannon')]
cbind_arare <- cbind(tmp_arare[meta_config_in_this$Donor_sra, c('X', 'shannon')], tmp_arare[meta_config_in_this$Previous_sra, c('X', 'shannon')])
if(unique(cbind_arare[,1] == meta_config_in_this$Donor_sra & cbind_arare[,3] == meta_config_in_this$Previous_sra)==c(TRUE)){
    cn_arare = cbind(meta_config_in_this$SRA_Sample, cbind_arare, meta_config_in_this$postfmt_symptoms, meta_config_in_this$pre_entro, meta_config_in_this$don_entro, meta_config_in_this$PRJ)
    colnames(cn_arare) <- c('names', 'Dnames', 'Darare', 'Pnames', 'Parare', 'symptom', 'pre_entro', 'don_entro', 'PRJ')
}


# tmp <- c_arare[DP_distance_dat$SRA_Sample]
if(unique(DP_distance_dat$SRA_Sample == (cn_arare$names))){
    DP_distance_dat_arare <- cbind(DP_distance_dat, cn_arare)
    # colnames(DP_distance_dat) <- c('SRA_Sample', 'distance', 'arare')
}

DP_distance_dat_arare$distance <- as.numeric(DP_distance_dat_arare$distance)



# val_before_donor_L6$B_alpha <- tmp_arare[meta_config_in_this$Previous_sra, c('shannon')]
val_before_donor_L6$distance <- DP_distance_dat_arare$distance #* 100
val_before_donor_L6$D_alpha <- DP_distance_dat_arare$Darare#tmp_arare[meta_config_in_this$Donor_sra, c('shannon')]

val_before_donor_L6$DB_alpha <- (DP_distance_dat_arare$Darare + DP_distance_dat_arare$Parare) 
# before_donor_L6$D_alpha <- as.numeric((c_arare))



val_before_donor_L6$pre_entro <- (as.factor(as.character(val_meta_config2$pre_entro)))

val_before_donor_L6$y <- as.factor(as.character(val_meta_config2$postfmt_symptoms))

val_feature_data <- val_before_donor_L6



##model evaluation
F1_auc_cal<-function(test_auc, validation_auc){
    F1_auc = 1/(1/test_auc + 1/validation_auc)
    return(F1_auc)
}


mean_impor <- function(impor){
  impor<-sort(rowMeans(impor), decreasing = T)
  impor<-as.data.frame(cbind(names(impor), impor))
  colnames(impor)<-c('id', 'Importance')
  impor$Importance <- as.numeric(as.character(impor$Importance))
  return(impor)
}

naive_rf <- function(training, testing, validation, feature_nb){
    training1 <- training[,feature_nb]
    cols <- ncol(training1)
    rf_naive <- randomForest(y ~ ., data=training1[,], importance=TRUE, proximity=F,
                                 ntree = 500,nPerm=10, type='classification')#-c((cols-2):c(cols-1))
    import_naive <- importance(rf_naive, type=import_type)##
    
    ##testing
    testing1 <- testing[, feature_nb]
    predicted_testing_naive <- predict(rf_naive, testing1[,-c(cols)], type="prob")#c(cols-2):
    testing_prediction_for_roc_curve_naive <- prediction(predicted_testing_naive[,2], testing1[,cols])
    
    auc_rf_naive_testing <- performance(testing_prediction_for_roc_curve_naive, "auc")@y.values[[1]]
    pref_naive <- performance(testing_prediction_for_roc_curve_naive, "tpr", "fpr")
    
    out_tp_naive_testing <- imputate_roc(pref_naive@x.values[[1]], pref_naive@y.values[[1]], seq_roc)
    
    ##validation
    validation1 <- validation[, feature_nb]
    predicted_value_naive <- predict(rf_naive, validation1[,-c(cols)], type="prob")#c(cols-2):
    prediction_for_roc_curve_naive <- prediction(predicted_value_naive[,2], validation1[,cols])
    
    auc_rf_naive <- performance(prediction_for_roc_curve_naive, "auc")@y.values[[1]]
    pref_naive <- performance(prediction_for_roc_curve_naive, "tpr", "fpr")
    
    out_tp_naive <- imputate_roc(pref_naive@x.values[[1]], pref_naive@y.values[[1]], seq_roc)
    
    F1_value <- F1_auc_cal(auc_rf_naive_testing, auc_rf_naive)
    
    out_rf <- list('auc'=c(auc_rf_naive_testing, c(1, out_tp_naive_testing)), 'valid_auc'=c(auc_rf_naive, c(1, out_tp_naive)), 'import_naive'=import_naive, 'F1_value'=F1_value, 'training_dat'=rownames(training), 'predict'= predicted_value_naive[,2],  'predict_name'=rownames(validation1))#pref@x.values[[1]], pref@y.values[[1]]
    return(out_rf)
}

combine_rf <- function(training, testing, validation, feature_b1_d1, feature_b1_d2, feature_b2_d1, feature_b2_d2){
    # cols
    ##
    training_b1 <- training[training[,'pre_entro'] %in% c('before1'), ]#feature_b1
    cols_b1 <- ncol(training_b1)
    
    training_b1_d1 <- training_b1[training_b1[,'don_entro'] %in% c('1', 1),feature_b1_d1]
    rf_b1_d1 <- randomForest(y ~ ., data=training_b1_d1, importance=TRUE, proximity=F,
                                 ntree = 500,nPerm=10, type='classification')
    import_b1_d1 <- importance(rf_b1_d1, type=import_type)##
    
    training_b1_d2 <- training_b1[training_b1 [,'don_entro'] %in% c('2', 2),feature_b1_d2]
    rf_b1_d2 <- randomForest(y ~ ., data=training_b1_d2, importance=TRUE, proximity=F,
                                 ntree = 500,nPerm=10, type='classification')
    import_b1_d2<- importance(rf_b1_d2, type=import_type)
 
    
    training_b2 <- training[training[,'pre_entro'] %in% c('before2'), ]#feature_b2#c(feature_b2_d1[,1:(length(feature_b2_d1)-1)], feature_b2_d2)
    # cols_b2 <- ncol(training_b2)
    # rf_b2 <- randomForest(y ~ ., data=training_b2, importance=TRUE, proximity=F,
                                 # ntree = 500,nPerm=10, type='classification')
    
    training_b2_d1 <- training_b2[training_b2[,'don_entro'] %in% c('1', 1),feature_b2_d1]
    rf_b2_d1 <- randomForest(y ~ ., data=training_b2_d1, importance=TRUE, proximity=F,
                                 ntree = 500,nPerm=10, type='classification')
    import_b2_d1 <- importance(rf_b2_d1, type=import_type)##
    
    training_b2_d2 <- training_b2[training_b2 [,'don_entro'] %in% c('2', 2),feature_b2_d2]
    rf_b2_d2 <- randomForest(y ~ ., data=training_b2_d2, importance=TRUE, proximity=F,
                                 ntree = 500,nPerm=10, type='classification')
    import_b2_d2<- importance(rf_b2_d2, type=import_type)##
    
    #testing
    testing_b1 <- testing[testing[,'pre_entro'] %in% c('before1'), ]#feature_b1
    
    testing_b1_d1 <- testing[testing[,'pre_entro'] %in% c('before1') & testing[,'don_entro'] %in% c('1', 1), feature_b1_d1]
    cols_b1_d1 <- ncol(testing_b1_d1)
    testing_b1_d2 <- testing[testing[,'pre_entro'] %in% c('before1') & testing[,'don_entro'] %in% c('2', 2), feature_b1_d2]
    cols_b1_d2 <- ncol(testing_b1_d2)
    
    predicted_testing_b1_d1 <- predict(rf_b1_d1, testing_b1_d1[,-cols_b1_d1],type="prob")
    predicted_testing_b1_d2 <- predict(rf_b1_d2, testing_b1_d2[,-cols_b1_d2],type="prob")
    
    
    testing_b2 <- testing[testing[,'pre_entro'] %in% c('before2'), ]#feature_b2
    testing_b2_d1 <- testing[testing[,'pre_entro'] %in% c('before2') & testing[,'don_entro'] %in% c('1', 1), feature_b2_d1]
    cols_b2_d1 <- ncol(testing_b2_d1)
    testing_b2_d2 <- testing[testing[,'pre_entro'] %in% c('before2') & testing[,'don_entro'] %in% c('2', 2), feature_b2_d2]
    cols_b2_d2 <- ncol(testing_b2_d2)
    
    predicted_testing_b2_d1 <- predict(rf_b2_d1, testing_b2_d1[,-cols_b2_d1],type="prob")
    predicted_testing_b2_d2 <- predict(rf_b2_d2, testing_b2_d2[,-cols_b2_d2],type="prob")
    
    # predicted_value<-(c(predicted_value_b1, predicted_value_b2))
    prediction_for_roc_curve <- prediction(c(predicted_testing_b1_d1[,2], predicted_testing_b1_d2[,2], predicted_testing_b2_d1[,2], predicted_testing_b2_d2[,2]), c(testing_b1_d1[,cols_b1_d1], testing_b1_d2[,cols_b1_d2], testing_b2_d1[,cols_b2_d1], testing_b2_d2[,cols_b2_d2]))
    
    auc_rf_testing<- performance(prediction_for_roc_curve, "auc")@y.values[[1]]
    
    pref <- performance(prediction_for_roc_curve, "tpr", "fpr")
    
    out_testing_tp <- imputate_roc(pref@x.values[[1]], pref@y.values[[1]], seq_roc)
    
    
    ##validation
    validation_b1 <- validation[validation[,'pre_entro'] %in% c('before1'), ]#feature_b1
    validation_b1_d1 <- validation[validation[,'pre_entro'] %in% c('before1') & validation[,'don_entro'] %in% c('1', 1), feature_b1_d1]
    val_cols1_d1 <- which(colnames(validation_b1_d1) == 'y')
    
    validation_b1_d2 <- validation[validation[,'pre_entro'] %in% c('before1') & validation[,'don_entro'] %in% c('2', 2), feature_b1_d2]
    val_cols1_d2 <- which(colnames(validation_b1_d2) == 'y')
    
    predicted_value_b1_d1 <- predict(rf_b1_d1, validation_b1_d1[,-val_cols1_d1],type="prob")
    predicted_value_b1_d2 <- predict(rf_b1_d2, validation_b1_d2[,-val_cols1_d2],type="prob")
    
    
    validation_b2 <- validation[validation[,'pre_entro'] %in% c('before2'), ]#feature_b2
    validation_b2_d1 <- validation[validation[,'pre_entro'] %in% c('before2') & validation[,'don_entro'] %in% c('1', 1), feature_b2_d1]
    val_cols2_d1 <- which(colnames(validation_b2_d1) == 'y')
    
    validation_b2_d2 <- validation[validation[,'pre_entro'] %in% c('before2') & validation[,'don_entro'] %in% c('2', 2), feature_b2_d2]
    val_cols2_d2 <- which(colnames(validation_b2_d2) == 'y')
    
    predicted_value_b2_d1 <- predict(rf_b2_d1, validation_b2_d1[,-val_cols2_d1],type="prob")
    predicted_value_b2_d2 <- predict(rf_b2_d2, validation_b2_d2[,-val_cols2_d2],type="prob")
    
    # predicted_value<-(c(predicted_value_b1, predicted_value_b2))
    prediction_for_roc_curve <- prediction(c(predicted_value_b1_d1[,2], predicted_value_b1_d2[,2], predicted_value_b2_d1[,2], predicted_value_b2_d2[,2]), c(validation_b1_d1[,val_cols1_d1], validation_b1_d2[,val_cols1_d2], validation_b2_d1[,val_cols2_d1], validation_b2_d2[,val_cols2_d2]))
    
    auc_rf <- performance(prediction_for_roc_curve, "auc")@y.values[[1]]
    
    pref <- performance(prediction_for_roc_curve, "tpr", "fpr")
    
    out_tp <- imputate_roc(pref@x.values[[1]], pref@y.values[[1]], seq_roc)
    
    tmp_perf <- ROCR::performance(prediction_for_roc_curve, "sens", "spec")
    df <- data.frame(cut = tmp_perf@alpha.values[[1]], sens = tmp_perf@x.values[[1]], spec = tmp_perf@y.values[[1]])
    best_value <- c(df[which.max(df$sens + df$spec), "cut"], df[which.min((1-df$sens)*(1-df$sens) + (1-df$spec)*(1-df$spec)), "cut"])
    # min((1 - sensitivities)^2 + (1- specificities)^2)
    #pref@alpha.values[[1]][which.max(out_tp - seq_roc)]
    
    F1_value <- F1_auc_cal(auc_rf_testing, auc_rf)
    

    out_rf <- list('auc'=c(auc_rf_testing, c(1, out_testing_tp)), 'valid_auc'=c(auc_rf, c(1, out_tp)), 'import_b1_d1'=import_b1_d1, 'import_b1_d2'=import_b1_d2, 'import_b2_d1'=import_b2_d1, 'import_b2_d2'=import_b2_d2, 'F1_value'=best_value, 'training_dat'=rownames(training), 'predict'=c(predicted_value_b1_d1[,2], predicted_value_b1_d2[,2], predicted_value_b2_d1[,2], predicted_value_b2_d2[,2]), 'predict_name'=c(rownames(validation_b1_d1), rownames(validation_b1_d2), rownames(validation_b2_d1), rownames(validation_b2_d2)))#pref@x.values[[1]], pref@y.values[[1]]
    return(out_rf)
}

training_rf_n <- function(repeats, feature_abun_dat, validation, feature_nb, final=FALSE){
    if(final == TRUE){
        training <- feature_abun_dat[,]
        testing <- feature_abun_dat[,]
        out <- naive_rf(training, testing, validation, feature_nb)
        return(out)
    }else{
        data_set_size <- floor(nrow(feature_abun_dat)/5)
        indexes <- sample(1:nrow(feature_abun_dat), size = data_set_size)
        training <- feature_abun_dat[-indexes,]
        testing <- feature_abun_dat[indexes,]
        
        out <- naive_rf(training, testing, validation, feature_nb)
        return(out)
    }
}

extract_info <- function(out_auc, repeats, flag){
    out_auc_extract<-list()
    for(i in 1:repeats){
        out_auc_extract[[i]] <- list()
        for(j in 1:flag){
            # out_auc_extract[[i]][[j]] <- list()
            out_auc_extract[[i]][[j]] <- out_auc[(i-1)*flag+j]
        }
    } 
    return(out_auc_extract)
    
}

navie_run<-function(feature_abun_dat, val_feature_data, feature_nb){
    ##filter features
    repeats <- 501
    out_auc <- sapply((1:repeats), training_rf_n, feature_abun_dat=feature_abun_dat, validation=val_feature_data, feature_nb=feature_nb)
    out_auc_1 <- extract_info(out_auc, repeats, 7)
    
    auc <- data.frame(sapply(1:repeats, function(i){out_auc_1[[i]][[1]]}))
    mean_auc <- rowMeans(auc)
    print(mean_auc[1:2])
    
    import_nb <- data.frame(sapply(1:repeats, function(i){out_auc_1[[i]][3]}))
    # import_b2 <- data.frame(sapply(1:repeats, function(i){out_auc_1[[i]][4]}))
    
    mean_impor_nb<-mean_impor(import_nb)
    # mean_impor_b2<-mean_impor(import_b2)
    
    repeats <- 501
    
    ##model train
    out_auc_fil30 <- sapply((1:repeats), training_rf_n, feature_abun_dat=feature_abun_dat, validation=val_feature_data, feature_nb=c(as.character(mean_impor_nb[1:nfeatures, 'id']),'y'))
    
    out_auc_fil30e <- extract_info(out_auc_fil30, repeats, 7)
    
    auc_fil30 <- data.frame(sapply(1:repeats, function(i){out_auc_fil30e[[i]][1]}))
    mean_auc_fil30 <- rowMeans(auc_fil30)
    print(mean_auc_fil30[1:2])
    
    import_nb_fil30 <- data.frame(sapply(1:repeats, function(i){out_auc_fil30e[[i]][3]}))
    # import_b2_fil30 <- data.frame(sapply(1:repeats, function(i){out_auc_fil30e[[i]][4]}))
    
    mean_impor_nb_fil30<-mean_impor(import_nb_fil30)
    # mean_impor_b2_fil30<-mean_impor(import_b2_fil30)
    
    val_auc_fil30 <- data.frame(sapply(1:repeats, function(i){out_auc_fil30e[[i]][2]}))
    mean_val_auc_fil30 <- rowMeans(val_auc_fil30)
    print(mean_val_auc_fil30[1:2])
    F1_fil30 <- data.frame(sapply(1:repeats, function(i){c(out_auc_fil30e[[i]][4])}))
    F1_auc_fil <- cbind(t(F1_fil30), t(auc_fil30[1,]), t(val_auc_fil30[1,]))
    
    
    
    repeats=100
    out_auc_final <- sapply((1:repeats), training_rf_n, feature_abun_dat=feature_abun_dat, validation=val_feature_data, feature_nb=c(as.character(mean_impor_nb[1:nfeatures, 'id']),'y'), final=TRUE)
    
    out_auc_finale <- extract_info(out_auc_final, repeats, 7)
    
    
    import_nb_final <- data.frame(sapply(1:repeats, function(i){out_auc_finale[[i]][3]}))
    
    mean_impor_nb_final<-mean_impor(import_nb_final)
    
    val_auc_final <- data.frame(sapply(1:repeats, function(i){out_auc_finale[[i]][2]}))
    mean_val_auc_final <- rowMeans(val_auc_final)
    print(mean_val_auc_final[1:2])

    
    predict_value <- rowMeans(data.frame(sapply(1:repeats, function(i){out_auc_finale[[i]][6]})))
    predict_value<- data.frame(cbind(names(predict_value), predict_value))
    tmp <- cbind(predict_value, val_feature_data[row.names(predict_value), c('pre_entro', 'don_entro', 'y')])
    
    return(list(mean_impor_nb_fil30, F1_auc_fil, mean_auc_fil30, mean_val_auc_fil30, tmp))
    
}


# cols <- ncol(feature_abun_dat)
training_rf <- function(repeats, feature_abun_dat, validation, feature_b1_d1, feature_b1_d2, feature_b2_d1, feature_b2_d2, final=FALSE){
    if(final == TRUE){
        training <- feature_abun_dat[,]
        testing <- feature_abun_dat[,]
        out <- combine_rf(training, testing, validation, feature_b1_d1, feature_b1_d2, feature_b2_d1, feature_b2_d2)
        return(out)
    }else{
        data_set_size <- floor(nrow(feature_abun_dat)/5)
        indexes <- sample(1:nrow(feature_abun_dat), size = data_set_size)
        training <- feature_abun_dat[-indexes,]
        testing <- feature_abun_dat[indexes,]
        
        out <- combine_rf(training, testing, validation, feature_b1_d1, feature_b1_d2, feature_b2_d1, feature_b2_d2)
        return(out)
    }
}

