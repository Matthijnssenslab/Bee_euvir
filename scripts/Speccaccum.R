library(vegan)

VIRdf <- read.csv("../data/Pcoa/flembee_df_coverage_lengthnorm.tsv", row.names=1)
counts <- t(VIRdf)
sp_con_random <- specaccum(counts, method = "random", ci=2)
plot(sp_con_random)
#Writes sites & richness to df
spmat <- sp_con_random$perm
write.csv(spmat,file ="../data/speccaccummat.csv", quote=F)
