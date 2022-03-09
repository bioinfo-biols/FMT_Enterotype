# engraftment ratio calculation R script v1.0
# Author: Ruiqiao He, UCAS
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('coin') (if there are confounder need to be blocked)
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('engraftment_ratio.R')
#
#       result <- engraftment(FMT_abundance, FMT_config, group = 'Outcomes', blocked = 0)
#
#   or:
#       result <- engraftment(FMT_abundance, FMT_config, group = 'Outcomes', blocked='Blocked')
#   
#
#   Input: signature matrix and mixture file, formatted as test file
#   Output: engraft ratio and its p with group (FMT outcomes)
#
#       Options:
#       1) group: the group factor used for hypothesis testing, default: 'Outcomes' as column name
#       2) block: the blocked confounder in hypothesis testing, default 0 for no varible to block
#

#


classify_colonizer <- function(x, abundance){
    index <- match(as.character(unlist(x)), rownames(abundance))
    index_L6 <- abundance[index,]
    
    index_fc<- apply(index_L6, 1, function(y, x){
        out_p <- log((y[3]/y[1])+1);
        out_d <- log((y[3]/y[2])+1);
        c(out_p, out_d, y[3])}, x=x)
    
    quadrant1 <- mean(index_fc[3, index_fc[1,] > index_fc[2,]], trim=0.05)*length(index_fc[3,index_fc[1,] > index_fc[2,]])
    
    quadrant2 <- mean(index_fc[3, index_fc[1,] < index_fc[2,]], trim=0.05)*length(index_fc[3,index_fc[1,] < index_fc[2,]])
    
    -log(quadrant2/quadrant1)
}

##engraftment ratio
cal_fc <- function(FMT_abundance, FMT_config){
    ## remove NA
    abundance = FMT_abundance + 1e-5
    FMT_config_sub <- FMT_config[,c("Before", "Donor", "After")]
    
    fold_change <- apply(FMT_config_sub, 1, classify_colonizer, abundance)

    return(fold_change)

}

engraft_test<-function(fold_change, FMT_config, group, blocked=0){

    
    if(blocked != 0){

        library(coin)
        
        fc_coin <- as.data.frame(cbind(as.numeric(fold_change), (FMT_config[,group]), FMT_config[,blocked]))
        colnames(fc_coin)<- c('fc', 'group', 'block')

        # rownames(fc_coin) <- (FMT_config$SRA_Sample)
        fc_coin$fc <- as.numeric(as.character(fc_coin$fc))
        prj_list <- fc_coin[, 'block']
        list <- NULL
        for(i in unique(prj_list)){
            if (length(prj_list[prj_list %in% i]) > 2){
            list <- c(list, i)
            }
        }
    
        fc_coin_f <- fc_coin[fc_coin[, 'block'] %in% list,]
        
        fc_coin_f$block <- as.factor(fc_coin_f$block)
        fc_coin_f$group <- as.factor(fc_coin_f$group)
        
        pval<-wilcox_test(fc ~ group | block, fc_coin_f)

    }else{

        fc_coin_f <- as.data.frame(cbind(as.numeric(fold_change), (FMT_config[,group])))
        colnames(fc_coin_f)<- c('fc', 'group')

        fc_coin_f$fc <- as.numeric(as.character(fc_coin_f$fc))
        
        pval<-wilcox.test(fc ~ group , fc_coin_f)
    }
    
    return(list( 'engraft_ratio'=fc_coin_f, 'ratio_p'=pval))
}

engraftment <- function(FMT_abundance, FMT_config, group = 'Outcomes', blocked=0){
    fold_change <- cal_fc(FMT_abundance, FMT_config)
    
    output<- engraft_test(fold_change, FMT_config, group, blocked)

    return(output)

}
