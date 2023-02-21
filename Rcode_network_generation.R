
# Library used for network generation
library(bnlearn)
library(parallel)
library(Hmisc)
library(xlsx)

# load expression of DE miRNAs and biclust mRNAs for each timepoint
load("miRNA_DE_DEseq2_exp.rda")
load("biclust_mrna_exp.rda")

set.seed(123)
# generate bayesian network using bnlearn package with r =10000 iterations and 48 cores
cl = makeCluster(48)
pb <- txtProgressBar(min = 0, max = length(biclust_mrna_exp), style = 3)
bn_bs_lst <- list()
# For loop to generate network for each timepoint 
for(i in 1:length(biclust_mrna_exp)){
  setTxtProgressBar(pb, i) # Output progress bar to screen
  temp <- rbind(biclust_mrna_exp[[i]], mirna_DE_exp[[i]])
  temp <- data.frame(t(temp))
  bn_bs_lst[[i]]<- boot.strength(temp, algorithm = "hc", R= 10000, cluster = cl) # measuring arc confidence
  save(bn_bs_lst, file = "bn_bs_lst_1k_vst.rda")
}

stopCluster(cl)

bn_bs_lst_1krep <- bn_bs_lst


# extract network in top 25 quantile #
Str_bn_bs_lst_1krep <- list()
for (i in 1:length(bn_bs_lst_1krep)){
  t <- subset(bn_bs_lst_1krep[[i]], bn_bs_lst_1krep[[i]]$strength !=0 & bn_bs_lst_1krep[[i]]$from !=bn_bs_lst_1krep[[i]]$to)
  Str_bn_bs_lst_1krep[[i]] <- subset(t, t$strength >= quantile(t$strength)[4] & t$direction >= quantile(t$direction)[4]) 
}


# number of interactions in each timepoint 
for(i in 1:5){print(dim(Str_bn_bs_lst_1krep[[i]]))}

# Fix naming
for(i in 1:length(Str_bn_bs_lst_1krep)){
  Str_bn_bs_lst_1krep[[i]]$from[grepl("rno.", Str_bn_bs_lst_1krep[[i]]$from)] <-  gsub("\\.","\\-",  Str_bn_bs_lst_1krep[[i]]$from[grepl("rno.", Str_bn_bs_lst_1krep[[i]]$from)])
  Str_bn_bs_lst_1krep[[i]]$to[grepl("rno.", Str_bn_bs_lst_1krep[[i]]$to)]<- gsub("\\.","\\-",  Str_bn_bs_lst_1krep[[i]]$to[grepl("rno.", Str_bn_bs_lst_1krep[[i]]$to)])
}

# Calculate correlations (r) and p-values between DE miRNA and SAMBA biclustered mRNA
cor_temp_lst <- list()
for(i in 1:length(biclust_mrna_exp)){
  cor_temp_lst[[i]]<- rcorr(t(mirna_DE_exp[[i]]), t(biclust_mrna_exp[[i]]), type="spearman")
}

# Replace , with . in names
for(i in 1:length(cor_temp_lst)){
  colnames(cor_temp_lst[[i]]$P)<- gsub(",",".", colnames(cor_temp_lst[[i]]$P))
  rownames(cor_temp_lst[[i]]$P) <- gsub(",",".", rownames(cor_temp_lst[[i]]$P))
  colnames(cor_temp_lst[[i]]$r)<- gsub(",",".", colnames(cor_temp_lst[[i]]$r))
  rownames(cor_temp_lst[[i]]$r) <- gsub(",",".", rownames(cor_temp_lst[[i]]$r))
}

sheets <- c("24hr","72hr","10day","24hAS","1mAS") # To allow naming of sheets

# Add correlations (r), p-values and timepoint variables to network
for(i in 1:length(Str_bn_bs_lst_1krep)){
  for(j in 1:nrow(Str_bn_bs_lst_1krep[[i]])){
    Str_bn_bs_lst_1krep[[i]][j,"r"] <- cor_temp_lst[[i]]$r[Str_bn_bs_lst_1krep[[i]][[j,1]], Str_bn_bs_lst_1krep[[i]][[j,2]]] 
    Str_bn_bs_lst_1krep[[i]][j,"p"] <- cor_temp_lst[[i]]$P[Str_bn_bs_lst_1krep[[i]][[j,1]], Str_bn_bs_lst_1krep[[i]][[j,2]]] 
    Str_bn_bs_lst_1krep[[i]][j,"timepoint"] <-  sheets[i]
  }
}

# Save the final network in rda format
save(Str_bn_bs_lst_1krep, file= "Str_bn_bs_lst_1krep.rda") 

# Save as the network in excel format  
for(i in 1:length(Str_bn_bs_lst_1krep)){
  write.xlsx(Str_bn_bs_lst_1krep[[i]], file="Str_bn_bs_lst_1krep.xlsx", row.names = FALSE, append=TRUE, sheetName=sheets[i])
}



