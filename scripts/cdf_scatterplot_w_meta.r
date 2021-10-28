setwd("~/Documents/pipeline_PlasTools/clean_zone/")
library(ggplot2)
library("report")

t_parameters_w_meta_nonredundant<- read.table(file="results/tables_genes_contig/parameters_ORF_aligned_w_meta_non_redundant.tsv", 
                                  header = T )
# cdf_GC_w_meta_non_redundant <- ecdf (t_parameters_w_meta_nonredundant$GC)
# plot( cdf_GC_w_meta_non_redundant, xlab="GC % in the ORF founding a hit " ,ylab = "Cumulative frequency", 
#       main= "CDF GC % in all genes from the whole metagenome \n having at least one match in the plasmid database", 
#       cex.main= 1,
#       col="#69b3a2")
# summary(t_parameters_w_meta_nonredundant)


t_stat_w_meta_nonredundant <- read.table(file="results/tables_genes_contig/table_contigid_nbgenes_stat_aplysina_w_meta_nonredundantdb.tsv", 
                                header = T )
head(t_stat_w_meta_nonredundant)
summary(t_stat_w_meta_nonredundant)

t_stat_w_meta_nonredundant_3genes <- read.table(file="results/tables_genes_contig/table_contigid_nbgenes_stat_aplysina_w_meta_nonredundantdb_min3genes_min1genematch_minproportion033_sortednbgenes.tsv", 
                                             header = T )

summary(t_stat_w_meta_nonredundant_3genes)
max(t_stat_len_w_meta_nonredundant$Length)
cdf_nbgenes <- ecdf (t_stat_len_w_meta_nonredundant$Nb_genes_matching_on_contig)
plot(cdf_nbgenes, xlab="Number of genes (ORF) matching per contig" ,ylab = "Cumulative frequency", 
     main= "CDF nb genes matching per contig in aplysina metagenome \n having at least one match in the plasmid database", 
     cex.main= 1, pch = 19)
mean (t_stat_len_w_meta_nonredundant$Nb_genes_matching_on_contig)
median (t_stat_len_w_meta_nonredundant$Nb_genes_matching_on_contig)

ggplot(t_stat_len_w_meta_nonredundant, aes (x=Proportion_genes.matching, y=Nb_genes_matching_on_contig)) +
  geom_point(color="#6B6B6B", size=2, alpha=0.9) +
  geom_rug(col="steelblue",alpha=0.1, size=1.5)+ 
  ggtitle("Correlation between length of the contig and the proportion of \ngenes matching the plasmid database and the total number gene on the contig")+
  ylab("Nb of genes on the contig matching") +xlab("Proportion of genes matching/total genes on contig ( between 1 and 0) ")
  
 ggplot(t_stat_len_w_meta_nonredundant, aes (x=Nb_genes, y=Length)) +
  geom_point(color="#6B6B6B", size=2, alpha=0.9) +
  geom_rug(col="steelblue",alpha=0.1, size=1.5)+ 
  ggtitle("Correlation between length of the contig and \nthe nb genes on the contig")+
  ylab("Length of the contig (bp)") +xlab("Nb genes on contig")
ggplot(t_stat_len_w_meta_nonredundant, aes (x=Nb_genes_matching_on_contig, y=Length)) +
  geom_point(color="#6B6B6B", size=2, alpha=0.9) +
  geom_rug(col="steelblue",alpha=0.1, size=1.5)+ 
  ggtitle("Correlation between length of the contig and \nthe nb genes matching the database on the contig")+
  ylab("Length of the contig (bp)") +xlab("Nb genes matching on contig")

pearson_test_nbgenes <- cor.test(t_stat_len_w_meta_nonredundant$Nb_genes, t_stat_len_w_meta_nonredundant$Length )
pearson_test_matches <- cor.test(t_stat_len_w_meta_nonredundant$Nb_genes_matching_on_contig, t_stat_len_w_meta_nonredundant$Length )
pearson_test_proportion <- cor.test(t_stat_len_w_meta_nonredundant$Proportion_genes.matching, t_stat_len_w_meta_nonredundant$Length )
pearson_test_nbgenes
pearson_test_matches
pearson_test_proportion

cdf_proportion_w_meta_non_redundant <- ecdf (t_stat_len_w_meta_nonredundant$Proportion_genes.matching)
plot( cdf_proportion_w_meta_non_redundant, xlab="proportion of ORF having a hit on a contig" ,ylab = "Cumulative frequency", 
      main= "CDF proportion of ORF having a hit on a contig \n  from the whole metagenome (having at least one match in the plasmid database)", 
      cex.main= 1,
      col="#69b3a2")

median(t_stat_len_w_meta_nonredundant$Proportion_genes.matching)
mean(t_stat_len_w_meta_nonredundant$Proportion_genes.matching)
#Re 2nd question: if I remember correctly - we (or maybe only I) were interested in knowing 
#how the distribution of hits in the plasmid RefSeq relates to contig length (in ORFs). 
#This issue came up since we have some contigs that have eg 50% or 100% hits in the plasmid DB 
#but these may correspond to contigs of only 2 or 1 ORFs. 
#Another possibility to look into this is to plot a CDF of contig length (in ORFs) for classes of contigs 
#according to the proportion of ORFs having hits in the plasmid DB, e.g., 100-90%, 90-80%, 80-70% etc.

#General distribution, all hits according to the length 
ggplot(t_stat_w_meta_nonredundant_3genes, aes (x=Nb_tot_hits, y=Nb_genes)) +
  geom_point(color="#6B6B6B", size=2, alpha=0.9) +
  ggtitle("Distribution nb genes of the contig and \nthe nb hits general the database on the contig")+
  ylab("Nb genes total") +xlab("Nb hits on contig")
median(t_stat_w_meta_nonredundant_3genes$Nb_tot_hits)

#Among the ORF (distribution)
ggplot(t_parameters_w_meta_nonredundant, aes (x=Nb_hits, y=Length_ORF)) +
  geom_point(color="#6B6B6B", size=2, alpha=0.9) +
  ggtitle("Distribution length of the ORF and \nthe nb hits general the database on the ORF")+
  ylab("Length of the ORF (bp)") +xlab("Nb hits on the ORF")
median(t_parameters_w_meta_nonredundant$Nb_hits)



#Nice blue #2998F9 nice green #69b3a2

#Color PPRmeta prediction and PlasForest (Beate)

