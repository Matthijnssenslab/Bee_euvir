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
Flemishbee_df_denovo <- read.delim("../data/Pcoa/Flemishbee_df_denovo.tsv", row.names=1)
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
acc_EUtax_Country <- acc_EUtax_Country_RED
acc_EUtax_Country$TaxRed <- NULL
Contiglen <- read_delim("../data/Pcoa/sracontig_length.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRAdf_norm <- sweep(SRAdf, 1, Contiglen$X2, '/')
SRAdf_norm <- t(SRAdf_norm)
#Write out. This thing is relatively large.
#write.csv(SRAdf_norm, file="../data/Pcoa/SRA_df_coverage_lengthnorm.tsv",quote = FALSE)

#We now have df's for the sra search, and 1 for the flemish bee samples, but we want them both together.
#This combination was done in jupyter notebook, in the 11_sharing book.
gzip = gzfile("../data/Pcoa/flembee_sra_df_coverage_lengthnorm.tsv.gz")
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
Rdf <- as.data.frame(rbind(c('Country',0.2658770,'Individual', 0.001),c('Country',0.2430856,'Cumulative', 0.002),c('EU family',0.1510465,'Individual',0.001),c('EU family',0.3384194,'Cumulative', 0.002)))
colnames(Rdf) <- c('Variable','Rsq','Rsqtype','P-value')

################### Plots ################### 
ggplot() + 
  geom_bar(stat="identity", aes(y=Rdf$Rsq, x=Rdf$Variable), color=Rdf$Rsqtype, position=position_dodge())

#Pcoa - Beevir

virpcoaplot <- ggplot() +
  theme_minimal() +
  scale_color_manual(values=c("#2ca02c","#d62728")) +
  geom_point(aes(x=vir_pcoa$vectors[,1], y=vir_pcoa$vectors[,2], col=statcol)) + 
  labs(y= "PCoA2", x = "PCoA1", color="Sample status")
virpcoaplot


#Pcoa - Sra
pcoa_col <- acc_EUtax_Country_RED[match(rownames(sra_pcoa$vectors), rownames(acc_EUtax_Country_RED)),]


srapcoaplot_taxfam <- ggplot() +
  theme_minimal() +
  scale_color_viridis_d() +
  geom_point(aes(x=sra_pcoa$vectors[,1], y=sra_pcoa$vectors[,2], col=pcoa_col$TaxRed)) +
  labs(y="PCoA2",x="PCoA1", color='Eukaryotic Family')
srapcoaplot_taxfam



srapcoaplot_country <- ggplot() +
  theme_minimal() +
  scale_color_viridis_d() +
  geom_point(aes(x=sra_pcoa$vectors[,1], y=sra_pcoa$vectors[,2], col=pcoa_col$Country)) +
  labs(y="PCoA2",x="PCoA1", color='Country')
srapcoaplot_country

