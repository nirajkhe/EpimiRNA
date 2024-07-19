# EpimiRNA
R code and accompanying data files to reproduce miRNA bayesian network analysis as performed in Khemka et al Sci Rep 14, 15313 (2024). https://doi.org/10.1038/s41598-024-66117-7 

File descriptions:
Network generation code: 
Rcode_network_generation.R : This R script is used to generate the Bayesian network using the provided data. It takes the below miRNA and mRNA expression data as input, and generates the belwo Network fiels as output.

miRNA and mRNA expression data:
miRNA_DE_DEseq2_exp.rda : This RData file contains the miRNA expression data for all significantly differentially expressed miRNAs, with a table for each time point.

biclust_mrna_exp.rda : This RData file contains the IDs and expression of biclustered mRNAs used for network generation.

miRNA_DE_DEseq2.rda: This RData file contains the miRNA differential expression table for each time point compared to the control (10 day).


Network files:
Str_bn_bs_lst_1krep.rda : This RData file contains the Bayesian network-based interactions for all mRNA/miRNA at each time point, along with the correlation and p-value of the correlation.

Str_bn_bs_lst_1krep.xlsx : This Excel file includes Bayesian network-based interactions for all mRNA/miRNA at each time point, with correlation and p-value of the correlation.

