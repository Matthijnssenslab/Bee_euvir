# Run PCoA's, this takes a while.
library(vegan)
library("phyloseq")
library(ggplot2)
library(ape)
library(readr)
library(viridis)
library(dplyr)
library(svglite)

#Regen_coverage matrix (w/o) for consistency with SRA covmat
Flemishbee_df_denovo <- read.delim("../data/Pcoa/Flemishbee_df_denovo_clean.tsv", row.names=1)
Contiglen <- read_delim("../data/Pcoa/sracontig_length.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
Flemishbee_df_norm <- sweep(Flemishbee_df_denovo, 1, Contiglen$X2, '/')
write.csv(Flemishbee_df_norm, file="../data/Pcoa/flembee_df_coverage_lengthnorm.tsv",quote = FALSE)

#Pcoa - Vir
VIRdf <- read.csv("../data/Pcoa/flembee_df_coverage_lengthnorm.tsv", row.names=1)
VIRdf <- t(VIRdf)
samstat <- read.csv("../data/Pcoa/Sample_status.csv", header=FALSE,row.names=1, sep=";")
colnames(samstat) <- c('Status')

vir_distmat <- vegdist(VIRdf, method='bray')
vir_pcoa <- pcoa(vir_distmat)
statcol <- samstat[match(rownames(vir_pcoa$vectors), rownames(samstat)),]
adonis(vir_distmat ~ statcol, permutations = 10000)



#Pcoa - SRA
gzip = gzfile("../data/Pcoa/SRA_df_denovo.tsv.gz")
SRAdf <- read.delim(gzip, row.names=1)
acc_EUtax_Country_RED <- read.csv("../data/Pcoa/acc_EUtax_Country.csv", row.names=1)
# Remove the four contigs listed as reagent contaminants:
SRAdf <- SRAdf[-which(rownames(SRAdf) %in% c("BeeP-25-2013_NODE_2878_length_686_cov_8_561576","BeeP-25-2013_NODE_7671_length_511_cov_8_142857",
                                    "BeeP-35-2013_NODE_6003_length_514_cov_1_240275","DRR029858.Contig_9800_10.032_length_927")), ]

acc_EUtax_Country <- acc_EUtax_Country_RED
acc_EUtax_Country$TaxRed <- NULL
Contiglen <- read_delim("../data/Pcoa/sracontig_length_cleaned.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRAdf_norm <- sweep(SRAdf, 1, Contiglen$X2, '/')
SRAdf_norm <- t(SRAdf_norm)
#Write out. This thing is relatively large.
#write.csv(SRAdf_norm, file="../data/Pcoa/SRA_df_coverage_lengthnorm.tsv",quote = FALSE)

#We now have df's for the sra search, and 1 for the flemish bee samples, but we want them both together.
#This combination was done in jupyter notebook, in the Virussharing notebook.
gzip = gzfile("../data/notebooks_out/flembee_sra_df_coverage_lengthnorm.tsv.gz")
SRAdf_norm <- read.csv(gzip, row.names=1)

sra_distmat <- vegdist(SRAdf_norm, method='bray')
sra_pcoa <- pcoa(sra_distmat)
pcoa_col <- acc_EUtax_Country[match(rownames(sra_pcoa$vectors), rownames(acc_EUtax_Country)),]




# Capscale (individual)
m <- acc_EUtax_Country[match(rownames(SRAdf_norm), rownames(acc_EUtax_Country)),]
m$Country[m$Country == ''] <- NA
#Location
#CAP for single R?.
cap_country <- capscale(SRAdf_norm ~ m[,1], distance = "bray", na.action=na.omit)
anova_country <- anova.cca(cap_country)
cap_tax <- capscale(SRAdf_norm ~ m[,2], distance = "bray", na.action=na.omit)
anova_tax <- anova.cca(cap_tax)
all <- c()
pval_country <- anova_country["Pr(>F)"][[1]][[1]] #get p-value
Fa_country <- anova_country["F"][[1]][[1]] #get F score
r2_country <- RsquareAdj(cap_country)[[1]] #get R2 (effect size)
r2adj_country <- RsquareAdj(cap_country)[[2]] #get adjusted R2 (effect size adjusted by eg number of levels in a factor)
all <- rbind(all,cbind(Fa_country,r2_country,r2adj_country,pval_country))
pval_tax <- anova_tax["Pr(>F)"][[1]][[1]] #get p-value
Fa_tax <- anova_tax["F"][[1]][[1]] #get F score
r2_tax <- RsquareAdj(cap_tax)[[1]] #get R2 (effect size)
r2adj_tax <- RsquareAdj(cap_tax)[[2]] #get adjusted R2 (effect size adjusted by eg number of levels in a factor)
all <- rbind(all,cbind(Fa_tax,r2_tax,r2adj_tax,pval_tax))
FDR = p.adjust(all[,"pval_country"],method="BH")
all = cbind(all,FDR)
colnames(all) <- c("F","r2","r2adj","p-value","FDR")
row.names(all) <-colnames(m)
sig.vars = row.names(all[all[,"FDR"] < 0.1,])
m0 <- na.exclude(m)
SRAdf_norm0 = SRAdf_norm[row.names(SRAdf_norm) %in% row.names(m0),]
#OrdiR2 for both (cumulative R?)
attach(m0)
mod0 = capscale(SRAdf_norm0 ~ 1, distance="bray") # Model with intercept only
mod1 = capscale(SRAdf_norm0 ~ ., data=m0, distance="bray") # Model with all explanatory variables
step.res<-ordiR2step(mod0, scope=formula(mod1), data=m0, direction="forward", Pin = 0.1, R2scope = FALSE, pstep = 1000, trace = T)
stepresanova <- as.data.frame(step.res$anova)
Rdf <- as.data.frame(rbind(c('Country',all[1,3],'Individual', all[1,4])
                           ,c('Country',stepresanova[1 ,1],'Cumulative', stepresanova[1 ,5]),
                           c('EU family',all[2,3],'Individual',all[2,4]),
                           c('EU family',stepresanova[2 ,1],'Cumulative', stepresanova[2 ,5])))
colnames(Rdf) <- c('Variable','Rsq','Rsqtype','P-value')
write.csv(Rdf, "../data/Pcoa/Rsq_pcoa.csv")

################### Plots ################### 
Rdf$Rsqtype <- as.factor(Rdf$Rsqtype)
Rdf$Rsq <- as.numeric(Rdf$Rsq)
ggplot() + 
  geom_bar(stat="identity", aes(y=Rdf$Rsq, x=Rdf$Variable), color=Rdf$Rsqtype, position=position_dodge())

#Get % var explained.
vir_pcoa$values$Relative_eig[1]
vir_pcoa$values$Relative_eig[2]

#Pcoa - Beevir
virpcoaplot <- ggplot() +
  theme_minimal() +
  scale_color_manual(values=c("#2ca02c","#d62728")) +
  geom_point(aes(x=vir_pcoa$vectors[,1], y=vir_pcoa$vectors[,2], col=statcol)) + 
  labs(y= "PCoA2 [8,5%] ", x = "PCoA1 [31,4%] ", color="Sample status")
virpcoaplot
ggsave("../Figures/Pcoa_bel.svg", virpcoaplot,  device="svg", dpi=300)


#Get % var explained.
sra_pcoa$values$Relative_eig[1]
sra_pcoa$values$Relative_eig[2]


#Pcoa - Sra
pcoa_col <- acc_EUtax_Country_RED[match(rownames(sra_pcoa$vectors), rownames(acc_EUtax_Country_RED)),]
srapcoaplot_taxfam <- ggplot() +
  theme_minimal() +
  scale_color_viridis_d() +
  geom_point(aes(x=sra_pcoa$vectors[,1], y=sra_pcoa$vectors[,2], col=pcoa_col$TaxRed)) +
  labs(y="PCoA2 [11.0%] ",x="PCoA1 [16.5%] ", color='Eukaryotic Family')
srapcoaplot_taxfam
ggsave("../Figures/Pcoa_sra_taxfam.svg", srapcoaplot_taxfam,  device="svg", dpi=300)


srapcoaplot_country <- ggplot() +
  theme_minimal() +
  scale_color_viridis_d() +
  geom_point(aes(x=sra_pcoa$vectors[,1], y=sra_pcoa$vectors[,2], col=pcoa_col$Country)) +
  labs(y="PCoA2 [11.0%] ",x="PCoA1 [16.5%] ", color='Country')
srapcoaplot_country
ggsave("../Figures/Pcoa_sra_country.svg", srapcoaplot_country,  device="svg", dpi=300)
