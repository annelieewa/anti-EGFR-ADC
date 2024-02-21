####################################################################
# Cetuximab ADC project - Anthony
# Author: Annelie Johansson
# Date: 19 November 2023
###################################################################

require("ggplot2")
require("dplyr")
require("reshape2")
require("gridExtra")
source("theme_manuscript.R")

load("data_for_analysis_all_BC_review_19Nov2023.RData")

genes_of_interest <- c("EGFR", "CCNA1", "CCND1", "CCNE1", "CCNH", "CCNT1",
                       "CDK2", "CDK4", "CDK6", "CDK7", "CDK9",
                       "ESR1", "ERBB2")

length(genes_of_interest) # 13
length(which(rownames(gex_guys) %in% genes_of_interest)) # 13
length(which(rownames(gex_scanb) %in% genes_of_interest)) # 13
length(which(rownames(gex_metabric) %in% genes_of_interest)) # 13
length(which(rownames(gex_tcga) %in% genes_of_interest)) # 13
length(which(rownames(gex_icgc) %in% genes_of_interest)) # 12

##################################################
#### ANALYSIS: EXPRESSION IN TNBC vs non_TNBC ####
##################################################

prepare_data <- function(pat_data, gex_data){
  temp <- as.data.frame(t(gex_data[which(rownames(gex_data) %in% genes_of_interest),]))
  temp <- as.data.frame(temp)
  temp$IHCSubtype2 <- pat_data$IHCSubtype2
  if(any(is.na(temp$IHCSubtype2))) { temp <- temp[-which(is.na(temp$IHCSubtype2)),] }
  return(temp)
}

tmp_plot_guys <- prepare_data(pat_guys, gex_guys)
tmp_plot_scanb <- prepare_data(pat_scanb, gex_scanb)
tmp_plot_metabric <- prepare_data(pat_metabric, gex_metabric)
tmp_plot_tcga <- prepare_data(pat_tcga, gex_tcga)
tmp_plot_icgc <- prepare_data(pat_icgc, gex_icgc)
tmp_plot_icgc$CCNA1 <- NA

tmp_plot_guys$cohort <- rep("Guy's", nrow(tmp_plot_guys))
tmp_plot_scanb$cohort <- rep("SCAN-B", nrow(tmp_plot_scanb))
tmp_plot_metabric$cohort <- rep("METABRIC", nrow(tmp_plot_metabric))
tmp_plot_tcga$cohort <- rep("TCGA", nrow(tmp_plot_tcga))
tmp_plot_icgc$cohort <- rep("ICGC", nrow(tmp_plot_icgc))

pathways <- genes_of_interest
cohorts <- c("Guy's", "SCAN-B","METABRIC", "TCGA", "ICGC")

temp <- rbind(tmp_plot_guys, tmp_plot_scanb, tmp_plot_metabric, tmp_plot_tcga, tmp_plot_icgc)
temp <- melt(temp, id.vars = c("IHCSubtype2", "cohort"))
temp$cohort <- factor(temp$cohort, levels = cohorts)
temp$variable <- factor(temp$variable, levels = pathways)

## add p-values:
dat_text <- dplyr::select(temp, cohort, variable)
dat_text <- unique(dat_text)
dat_text$pval_ER <- rep(NA, nrow(dat_text))
dat_text$pval_HER2 <- rep(NA, nrow(dat_text))
dat_text$variable <- factor(dat_text$variable, levels = pathways)
dat_text$IHCSubtype2 <- rep("TNBC", nrow(dat_text))

for(i in 1:nrow(dat_text)){
  temp_now <- temp[temp$variable == dat_text$variable[i] & temp$cohort == dat_text$cohort[i],]
  if(!all(is.na(temp_now$value))){
    dat_text$pval_ER[i] <- wilcox.test(temp_now$value[temp_now$IHCSubtype2 == "TNBC"],
                                       temp_now$value[temp_now$IHCSubtype2 == "ER+"], alternative = "two.sided")$p.value
    dat_text$pval_HER2[i] <- wilcox.test(temp_now$value[temp_now$IHCSubtype2 == "TNBC"],
                                       temp_now$value[temp_now$IHCSubtype2 == "HER2+"], alternative = "two.sided")$p.value
  }
}
dat_text$pval_ER_text <- ifelse(dat_text$pval_ER < 0.001, "***",
                             ifelse(dat_text$pval_ER < 0.01, "**",
                                    ifelse(dat_text$pval_ER < 0.05, "*", "NS")))
dat_text$pval_HER2_text <- ifelse(dat_text$pval_HER2 < 0.001, "***",
                                ifelse(dat_text$pval_HER2 < 0.01, "**",
                                       ifelse(dat_text$pval_HER2 < 0.05, "*", "NS")))

## set y:
dat_text$y <- NA
for(i in 1:nrow(dat_text)){
  values <- temp$value[which(temp$cohort == dat_text$cohort[i] & temp$variable == dat_text$variable[i])]
  if(!all(is.na(values))){ dat_text$y[i] <- max(values)*1.05 }
}

temp$IHCSubtype2 <- factor(temp$IHCSubtype2, levels = c("TNBC", "HER2+", "ER+"))

### Per gene
for(gene in unique(temp$variable)){
  pdf(paste0("EGFR/Figure_IHC_subtypes_", gene, ".pdf"), width = 2, height = 8)
  g <- ggplot(temp[temp$variable == gene,], aes(x = IHCSubtype2, y = value, fill = IHCSubtype2)) +
    geom_boxplot(show.legend = FALSE) +
    facet_grid(rows = vars(cohort), scales = "free_y") +
    scale_fill_manual(values = c("indianred", "palevioletred", "dodgerblue3")) +
    scale_y_continuous("Gene expression (log2)", expand = c(0.05, 0.05)) +
    xlab("") +
    theme_Publication() + 
    theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1)) +
    geom_text(dat_text[dat_text$variable == gene,], mapping = aes(x = 2, y = y, label = pval_HER2_text)) +
    geom_text(dat_text[dat_text$variable == gene,], mapping = aes(x = 3, y = y, label = pval_ER_text))
  print(g)
  dev.off() 
}


### Extract stats ----

my_fun <- function(x, cohort){
  stats_out <- stats[which(stats$cohort == cohort),]
  rownames(stats_out) <- stats_out$variable
  stats_out_HER2 <- stats_out[genes_of_interest,"pval_HER2"]
  stats_out_ER <- stats_out[genes_of_interest,"pval_ER"]
  names(stats_out_HER2) <- names(stats_out_ER) <- genes_of_interest
  return(list(stats_out_HER2, stats_out_ER))
}

stats <- dplyr::select(dat_text, cohort, variable, pval_HER2, pval_ER)

stats_guys <- my_fun(stats, "Guy's")
stats_scanb <- my_fun(stats, "SCAN-B")
stats_metabric <- my_fun(stats, "METABRIC")
stats_tcga <- my_fun(stats, "TCGA")
stats_icgc <- my_fun(stats, "ICGC")

### HER2
stats <- cbind(stats_guys[[1]], stats_scanb[[1]], stats_metabric[[1]], stats_tcga[[1]], stats_icgc[[1]])
colnames(stats) <- c("Guy's", "SCAN-B", "METABRIC", "TCGA", "ICGC")
stats <- t(apply(stats, 1, function(x) { formatC(x, format = "f", digits = 3)}))
stats[which(stats == "0.000")] <- "<0.001"
stats
#       Guy's    SCAN-B   METABRIC TCGA     ICGC   
# EGFR  "<0.001" "<0.001" "<0.001" "<0.001" "0.016"
# CCNA1 "0.009"  "<0.001" "<0.001" "<0.001" "  NA" 
# CCND1 "<0.001" "<0.001" "0.012"  "<0.001" "0.286"
# CCNE1 "<0.001" "<0.001" "<0.001" "<0.001" "0.428"
# CCNH  "<0.001" "<0.001" "<0.001" "<0.001" "0.654"
# CCNT1 "0.002"  "<0.001" "0.085"  "0.002"  "0.111"
# CDK2  "0.102"  "<0.001" "0.378"  "0.022"  "0.739"
# CDK4  "0.007"  "<0.001" "0.011"  "<0.001" "0.973"
# CDK6  "0.001"  "<0.001" "<0.001" "<0.001" "0.052"
# CDK7  "0.001"  "<0.001" "0.008"  "<0.001" "0.256"
# CDK9  "0.464"  "<0.001" "0.002"  "<0.001" "0.266"
# ESR1  "<0.001" "<0.001" "<0.001" "<0.001" "0.010"
# ERBB2 "0.142"  "<0.001" "<0.001" "<0.001" "0.004"


### ER
stats <- cbind(stats_guys[[2]], stats_scanb[[2]], stats_metabric[[2]], stats_tcga[[2]], stats_icgc[[2]])
colnames(stats) <- c("Guy's", "SCAN-B", "METABRIC", "TCGA", "ICGC")
stats <- t(apply(stats, 1, function(x) { formatC(x, format = "f", digits = 3)}))
stats[which(stats == "0.000")] <- "<0.001"
stats
#       Guy's    SCAN-B   METABRIC TCGA     ICGC    
# EGFR  "0.080"  "<0.001" "<0.001" "<0.001" "<0.001"
# CCNA1 "0.133"  "<0.001" "<0.001" "<0.001" "  NA"  
# CCND1 "0.681"  "<0.001" "<0.001" "<0.001" "<0.001"
# CCNE1 "0.022"  "<0.001" "<0.001" "<0.001" "<0.001"
# CCNH  "0.233"  "<0.001" "<0.001" "<0.001" "<0.001"
# CCNT1 "0.560"  "<0.001" "0.002"  "<0.001" "0.069" 
# CDK2  "0.052"  "<0.001" "0.008"  "<0.001" "<0.001"
# CDK4  "0.083"  "<0.001" "<0.001" "<0.001" "<0.001"
# CDK6  "0.028"  "<0.001" "<0.001" "<0.001" "<0.001"
# CDK7  "0.035"  "<0.001" "<0.001" "<0.001" "0.003" 
# CDK9  "0.093"  "<0.001" "<0.001" "<0.001" "0.753" 
# ESR1  "0.156"  "<0.001" "<0.001" "<0.001" "<0.001"
# ERBB2 "<0.001" "<0.001" "<0.001" "<0.001" "<0.001"


################################################
#### ANALYSIS: EXPRESSION IN PAM50 subtypes ####
################################################

### Missing for ICGC

prepare_data <- function(pat_data, gex_data){
  temp <- as.data.frame(t(gex_data[which(rownames(gex_data) %in% genes_of_interest),]))
  temp <- as.data.frame(temp)
  temp$PAM50 <- pat_data$PAM50
  if(any(is.na(temp$PAM50))) { temp <- temp[-which(is.na(temp$PAM50)),] }
  return(temp)
}

tmp_plot_guys <- prepare_data(pat_guys, gex_guys)
tmp_plot_scanb <- prepare_data(pat_scanb, gex_scanb)
tmp_plot_metabric <- prepare_data(pat_metabric, gex_metabric)
tmp_plot_tcga <- prepare_data(pat_tcga, gex_tcga)

tmp_plot_guys$cohort <- rep("Guy's", nrow(tmp_plot_guys))
tmp_plot_scanb$cohort <- rep("SCAN-B", nrow(tmp_plot_scanb))
tmp_plot_metabric$cohort <- rep("METABRIC", nrow(tmp_plot_metabric))
tmp_plot_tcga$cohort <- rep("TCGA", nrow(tmp_plot_tcga))

pathways <- genes_of_interest
cohorts <- c("Guy's", "SCAN-B","METABRIC", "TCGA")

temp <- rbind(tmp_plot_guys, tmp_plot_scanb,tmp_plot_metabric, tmp_plot_tcga)
temp$PAM50[which(temp$PAM50 == "HER2")] <- "Her2"
temp$PAM50[which(temp$PAM50 == "LuminalA")] <- "LumA"
temp$PAM50[which(temp$PAM50 == "LuminalB")] <- "LumB"
temp$PAM50[which(temp$PAM50 == "Normal-Like")] <- "Normal"

temp <- melt(temp, id.vars = c("PAM50", "cohort"))
temp$cohort <- factor(temp$cohort, levels = cohorts)
temp$variable <- factor(temp$variable, levels = pathways)

## add p-values 
dat_text <- dplyr::select(temp, cohort, variable)
dat_text <- unique(dat_text)
dat_text$pval_her2 <- rep(NA, nrow(dat_text))
dat_text$pval_luma <- rep(NA, nrow(dat_text))
dat_text$pval_lumb <- rep(NA, nrow(dat_text))
dat_text$pval_normal <- rep(NA, nrow(dat_text))
dat_text$variable <- factor(dat_text$variable, levels = pathways)
dat_text$PAM50 <- rep("Basal", nrow(dat_text))
for(i in 1:nrow(dat_text)){
  temp_now <- temp[temp$variable == dat_text$variable[i] & temp$cohort == dat_text$cohort[i],]
  if(!all(is.na(temp_now$value))){
    dat_text$pval_her2[i] <- wilcox.test(temp_now$value[temp_now$PAM50 == "Basal"],
                                         temp_now$value[temp_now$PAM50 == "Her2"], alternative = "two.sided")$p.value
    dat_text$pval_luma[i] <- wilcox.test(temp_now$value[temp_now$PAM50 == "Basal"],
                                         temp_now$value[temp_now$PAM50 == "LumA"], alternative = "two.sided")$p.value
    dat_text$pval_lumb[i] <- wilcox.test(temp_now$value[temp_now$PAM50 == "Basal"],
                                         temp_now$value[temp_now$PAM50 == "LumB"], alternative = "two.sided")$p.value
    dat_text$pval_normal[i] <- wilcox.test(temp_now$value[temp_now$PAM50 == "Basal"],
                                           temp_now$value[temp_now$PAM50 == "Normal"], alternative = "two.sided")$p.value
  }
}
dat_text$pval_text_her2 <- ifelse(dat_text$pval_her2 < 0.001, "***",
                                  ifelse(dat_text$pval_her2 < 0.01, "**", ifelse(dat_text$pval_her2 < 0.05, "*", "")))
dat_text$pval_text_luma <- ifelse(dat_text$pval_luma < 0.001, "***",
                                  ifelse(dat_text$pval_luma < 0.01, "**", ifelse(dat_text$pval_luma < 0.05, "*", "")))
dat_text$pval_text_lumb <- ifelse(dat_text$pval_lumb < 0.001, "***",
                                  ifelse(dat_text$pval_lumb < 0.01, "**", ifelse(dat_text$pval_lumb < 0.05, "*", "")))
dat_text$pval_text_normal <- ifelse(dat_text$pval_normal < 0.001, "***",
                                    ifelse(dat_text$pval_normal < 0.01, "**", ifelse(dat_text$pval_normal < 0.05, "*", "")))


## set y:
dat_text$y <- NA
for(i in 1:nrow(dat_text)){
  values <- temp$value[which(temp$cohort == dat_text$cohort[i] & temp$variable == dat_text$variable[i])]
  if(!all(is.na(values))){ dat_text$y[i] <- max(values)*1.05 }
}

### Per gene
for(gene in unique(temp$variable)){
  pdf(paste0("EGFR/Figure_PAM50_subtypes_", gene, ".pdf"), width = 2.5, height = 8)
  g <- ggplot(temp[temp$variable == gene,], aes(x = PAM50, y = value, fill = PAM50)) +
    geom_boxplot(show.legend = FALSE) +
    facet_grid(rows = vars(cohort), scales = "free_y") +
    scale_fill_manual(values = c("red", "pink", "blue", "lightblue", "green")) +
    scale_y_continuous("Gene expression (log2)", expand = c(0.05, 0.05)) +
    xlab("") +
    theme_Publication() + 
    theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1)) +
    geom_text(dat_text[dat_text$variable == gene,], mapping = aes(x = 2, y = y, label = pval_text_her2)) +
    geom_text(dat_text[dat_text$variable == gene,], mapping = aes(x = 3, y = y, label = pval_text_luma)) +
    geom_text(dat_text[dat_text$variable == gene,], mapping = aes(x = 4, y = y, label = pval_text_lumb)) +
    geom_text(dat_text[dat_text$variable == gene,], mapping = aes(x = 5, y = y, label = pval_text_normal))
  print(g)
  dev.off() 
}


### Extract stats ----
stats <- dplyr::select(dat_text, cohort, variable, pval_her2, pval_luma, pval_lumb, pval_normal)

stats_guys <- stats[which(stats$cohort == "Guy's"),]
rownames(stats_guys) <- stats_guys$variable
stats_guys <- t(stats_guys[c("EGFR", "CCNA1", "CCNE1", "CDK2"),-c(1,2)])
colnames(stats_guys) <- paste0("Guys_", colnames(stats_guys))
stats_scanb <- stats[which(stats$cohort == "SCAN-B"),]
rownames(stats_scanb) <- stats_scanb$variable
stats_scanb <- t(stats_scanb[c("EGFR", "CCNA1", "CCNE1", "CDK2"),-c(1,2)])
colnames(stats_scanb) <- paste0("SCANB_", colnames(stats_scanb))
stats_metabric <- stats[which(stats$cohort == "METABRIC"),]
rownames(stats_metabric) <- stats_metabric$variable
stats_metabric <- t(stats_metabric[c("EGFR", "CCNA1", "CCNE1", "CDK2"),-c(1,2)])
colnames(stats_metabric) <- paste0("METABRIC_", colnames(stats_metabric))
stats_tcga <- stats[which(stats$cohort == "TCGA"),]
rownames(stats_tcga) <- stats_tcga$variable
stats_tcga <- t(stats_tcga[c("EGFR", "CCNA1", "CCNE1", "CDK2"),-c(1,2)])
colnames(stats_tcga) <- paste0("TCGA_", colnames(stats_tcga))

stats <- cbind(stats_guys, stats_scanb, stats_metabric, stats_tcga)
stats <- t(apply(stats, 1, function(x) { formatC(x, format = "f", digits = 3)}))
stats[which(stats == "0.000")] <- "<0.001"

t(stats)
#             pval_her2 pval_luma pval_lumb pval_normal
# Guys_EGFR      "0.001"   "<0.001"  "<0.001"  "0.905"    
# Guys_CCNA1     "<0.001"  "0.004"   "0.001"   "0.334"    
# Guys_CCNE1     "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# Guys_CDK2      "0.077"   "<0.001"  "0.530"   "<0.001"   
# SCANB_EGFR     "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# SCANB_CCNA1    "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# SCANB_CCNE1    "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# SCANB_CDK2     "<0.001"  "<0.001"  "0.511"   "<0.001"   
# METABRIC_EGFR  "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# METABRIC_CCNA1 "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# METABRIC_CCNE1 "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# METABRIC_CDK2  "<0.001"  "<0.001"  "0.292"   "<0.001"   
# TCGA_EGFR      "<0.001"  "<0.001"  "<0.001"  "0.002"    
# TCGA_CCNA1     "<0.001"  "<0.001"  "<0.001"  "0.002"    
# TCGA_CCNE1     "<0.001"  "<0.001"  "<0.001"  "<0.001"   
# TCGA_CDK2      "0.017"   "<0.001"  "0.633"   "<0.001"  

