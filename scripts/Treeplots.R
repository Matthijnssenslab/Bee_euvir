library("ggtree")
library("ggplot2")
library("phytools")
library("viridis")
library("svglite")
library("gridExtra")
library("treeio")

setwd("/Users/deboutte/Bioinf/repo/Bee_euvir")


# Bunya -------------------------------------------------------------------
bunya <- read.mrbayes("data/cons_trees/bunya.nex.con.tre")

# Bunya Colors ------------------------------------------------------------
cls <- list(Unclassified = c('KX884775.1_1',
                             'KX884805.1_1',
                             'KX884807.1_1',
                             'MF176242.1_1',
                             'MF176293.1_1',
                             'MF176313.1_1',
                             'MF176329.1_1',
                             'MF176353.1_1',
                             'MF176372.1_1',
                             'MN168168.1_1',
                             'NC_040759.1_1',
                             'NC_032151.1_1'),
            Phasmaviridae = c('KJ434182.1_1',
                              'KM817697.1_1',
                              'KM817698.1_2',
                              'MH822966.1_1',
                              'NC_031307.1_1',
                              'NC_031312.1_2',
                              'NC_034462.1_1',
                              'NC_043642.1_1'),
            Peribunyaviridae = c('KY354237.1_1', 'SRR5109820.Contig_13922_1091.96_length_6449_1'),
            SRA = c('SRR1812786.Contig_29604_62.7809_length_6336_1','SRR2664950.Contig_23840_14.7864_length_6215_1'))
cols <- c("0" ="#7f7f7f", "Unclassified" = "#7f7f7f","SRA" = "#ff7f0e",
          "Phasmaviridae" = "#2ca02c","Peribunyaviridae" = "#2ca02c")

# Bunya Plot --------------------------------------------------------------
bunya@phylo <- midpoint.root(bunya@phylo)
bunya <- groupOTU(bunya, cls)
bunyatree <- ggtree(bunya, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0, 2) +
  geom_treescale() + 
  scale_colour_manual(values = cols)
bunyatree
ggsave("treePlots/bunya.svg",bunyatree, device = "svg", dpi=300)
ggsave("treePlots/bunya.pdf",bunyatree, device = "pdf", dpi=300)

# Orthomyxo_Bee -----------------------------------------------------------
orthobee <- read.mrbayes("data/cons_trees/orthomyxo_bee.nex.con.tre")


# Orthomyxo_Bee Colors ----------------------------------------------------
cls <- list(NODE = c('BP46_PB2PB1PANP','BP47_PB2PB1PANP','BP49_PB2PB1PANP'),
            Influenza = c('AB126193_Influenza_C','FJ969536_Influenza_A','M14880_Influenza_B'),
            Thogotolike = c('GU969313_Dhori','HM627174_Jos','KC506156_Upolu','KC506162_Aransas_Bay','KU708253_Bourbon','Y17873_Thogoto'),
            Quaranja = c('KM114304_Wellfleet_Bay','KX883865_Wuhan_Mosquito_4','KX883876_Jingshan_fly_1','KX883879_Sanxia_water_Strider_3'),
            Isavirus = c('KU587572_Infectious_salmon_anemia'),
            Unclassified = c('KX882061_Rainbow_trout_orthomyxovirus-1','KX882069_Steelhead_trout_orthomyxovirus-1','KX883845_Beihai_orthomyxo-like_1','KX883859_Hubei_earwig_1','KX883884_Hubei_orthoptera_6','KX949591_Sinu','Hubei_orthomyxo-like_2'))
cols <- c("0" ="#7f7f7f", "Unclassified" = "#7f7f7f","NODE" = "#ff7f0e", "Influenza" = "#2ca02c","Thogotolike"="#2ca02c", "Quaranja"="#2ca02c", "Isavirus"="#2ca02c")
cols2 <- c("0"="blue","Unclassified" = "#7f7f7f","NODE" = "#ff7f0e", "Influenza" = "red","Thogotolike"="purple", "Quaranja"="yellow", "Isavirus"="black")

# Orthomyxo_Bee Plot ------------------------------------------------------
orthobee@phylo <- midpoint.root(orthobee@phylo)
orthobee <- groupOTU(orthobee, cls)
orthobeetree <- ggtree(orthobee, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,4) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
orthobeetree
ggsave("treePlots/orthobee.svg",orthobeetree, device = "svg", dpi=300)
ggsave("treePlots/orthobee.pdf",orthobeetree, device = "pdf", dpi=300)


# Orthomyxo_SRA ---------------------------------------------------------------
orthosra <- read.mrbayes("data/cons_trees/orthomyxo_sra.nex.con.tre")

# Orthomyxo_SRA Colors ----------------------------------------------------
cls <- list(NODE = c('BP46_NODE_59_length_2296_cov_1606.73_PB1_1'),
            SRA = c('SRR2064700.Contig_7931_28.6803_length_2416_1',
                    'SRR2916026.Contig_21490_12.8872_length_2343_1',
                    'SRR3993646.Contig_4218_61.3183_length_2345_1',
                    'SRR4301486.Contig_8416_72.6801_length_2424_1',
                    'SRR5885442.Contig_31234_55.8986_length_1924_1',
                    'SRR6001380.Contig_6081_137. 105_length_2095_1',
                    'SRR770031.Contig_2318_75.597_length_2397_1'),
            delta = c('KM392504.1',
                      'KX768818.1'),
            alpha = c('CY035132.1',
                      'CY103874.1',
                      'JQ994252.1',
                      'KJ856190.1'),
            beta = c('CY171757.1',
                     'D00004.1',
                     'EF626634.1',
                     'KX351451.1',
                     'MH684310.1'),
            gamma = c('NC_006308.2'),
            thogoto = c('AF004985.1',
                        'HM627170.1',
                        'KC506157.1',
                        'KC506163.1',
                        'KU708254.1',
                        'KX670390.1',
                        'LC010982.1',
                        'M65866.1',
                        'NC_040731.1'),
            quaranja = c('KM114305.1',
                         'KM817616.1',
                         'KM817617.1',
                         'KM817619.1',
                         'KM817620.1',
                         'KM817621.1',
                         'KM817622.1',
                         'KM817623.1',
                         'KM817624.1',
                         'KM817625.1',
                         'KM817626.1',
                         'KM817627.1',
                         'KX883875.1',
                         'MG770333.1',
                         'MH267793.1',
                         'MH558140.1',
                         'MN053836.1',
                         'MN053837.1',
                         'MN053838.1'),
            unclass = c('JQ928944.1',
                        'KX883844.1',
                        'KX883858.1',
                        'KX883868.1',
                        'KX883873.1',
                        'KX883874.1',
                        'KX883884.1',
                        'KX898491.1',
                        'KX949590.1',
                        'MF176325.1',
                        'MG600038.1',
                        'MG972988.1',
                        'MG972993.1',
                        'MK026596.1',
                        'MK227173.1',
                        'MN167480.1'))
cols <- c("0" ="#7f7f7f", "unclass" = "#7f7f7f","NODE" = "#ff7f0e","SRA"="#1f77b4", "gamma"="#2ca02c","alpha"="#2ca02c","beta"="#2ca02c","delta"="#2ca02c","thogoto"="#2ca02c","quaranja"="#2ca02c")
cols2 <- c("0" ="#7f7f7f", "unclass" = "#7f7f7f","NODE" = "#ff7f0e","SRA"="#1f77b4", "gamma"="red","alpha"="red","beta"="red","delta"="red","thogoto"="blue","quaranja"="purple")

# Orthomyxo_SRA Plot ------------------------------------------------------
orthosra@phylo <- midpoint.root(orthosra@phylo)
orthosra <- groupOTU(orthosra, cls)
orthosratree <- ggtree(orthosra, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,4) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
orthosratree
ggsave("treePlots/orthosra.svg",orthosratree, device = "svg", dpi=300)
ggsave("treePlots/orthosra.pdf",orthosratree, device = "pdf", dpi=300)


# Partiti -----------------------------------------------------------------
partiti <- read.mrbayes("data/cons_trees/partiti.nex.con.tre")

# Partiti Colors ----------------------------------------------------------
cls <- list(NODE = c('BP15_NODE_21_length_1989_cov_11_700214_1',
                     'BP23_NODE_181_length_1578_cov_4_697323',
                     'BP2_NODE_11_length_2442_cov_8_707023',
                     'BP48_NODE_37_length_1861_cov_5_466667_1',
                     'BeeP-11-2013_NODE_907_length_1923_cov_7_401408_1',
                     'BeeP-32-2013_NODE_32_length_2325_cov_12_654359_1',
                     'BeeP-32-2013_NODE_58_length_1715_cov_25_394994_1',
                     'BeeP-33-2013_NODE_5_length_2414_cov_84_842533_1',
                     'BeeP-35-2013_NODE_401_length_2385_cov_62_380849_1',
                     'BeeP-35-2013_NODE_620_length_1910_cov_20_108565_1',
                     'BeeP-35-2013_NODE_899_length_1547_cov_109_939456_2',
                     'BeeP-37-2013_NODE_115_length_2447_cov_879_975527_1',
                     'BeeP-37-2013_NODE_116_length_2446_cov_659_816800_1',
                     'BeeP-37-2013_NODE_164_length_1939_cov_93_179914_1',
                     'BeeP-44-2013_NODE_110_length_1821_cov_96_515482_1',
                     'BeeP-44-2013_NODE_146_length_1607_cov_51_104575_1',
                     'BeeP-44-2013_NODE_94_length_1947_cov_31_678610_1',
                     'BeeP-44-2013_NODE_99_length_1913_cov_780_763617_1',
                     'BeeP-47-2013_NODE_67_length_1854_cov_32_899268_1',
                     'BeeP-47-2013_NODE_88_length_1610_cov_66_399870_1',
                     'BeeP-49-2013_NODE_249_length_1887_cov_3011_785635_1',
                     'BeeP-49-2013_NODE_255_length_1857_cov_503_888202_1',
                     'BeeP-49-2013_NODE_262_length_1819_cov_24_860505_1',
                     'BeeP-49-2013_NODE_307_length_1700_cov_116_550832_1'),
            SRA = c('ERR1354112.Contig_10743_83.2701_length_2186_2',
                    'SRR4002881.Contig_4619_145.602_length_2331_1'),
            Alpha = c('NC_038827.1_1',
                      'NC_038829.1_1',
                      'NC_038831.1_1',
                      'NC_038834.1_1',
                      'NC_038836.1_2'),
            Beta = c('NC_003801.1_1',
                     'NC_007537.1_1',
                     'NC_013109.1_1',
                     'NC_021094.1_1',
                     'NC_021096.1_1',
                     'NC_021098.1_1',
                     'NC_021147.1_1',
                     'NC_028251.1_1',
                     'NC_030882.1_1',
                     'NC_031134.1_1',
                     'NC_035485.1_1',
                     'NC_038837.1_1'),
            Gamma = c('NC_005976.2_1',
                      'NC_040387.1_1'),
            Unclassified = c('NC_013014.1_1',
                             'NC_021873.1_1',
                             'NC_030875.1_1',
                             'NC_030878.1_1',
                             'NC_030889.1_1',
                             'NC_034514.1_1',
                             'NC_040483.1_1',
                             'NC_040496.1_1'))
cols <- c("NODE" = "#ff7f0e","SRA" = "#1f77b4","Alpha" = "#2ca02c","Beta" = "#2ca02c","Gamma" = "#2ca02c","Unclassified" = "#7f7f7f", "0" = "#7f7f7f")
cols2 <- c("NODE" = "#ff7f0e","SRA" = "#1f77b4","Alpha" = "red","Beta" = "blue","Gamma" = "purple","Unclassified" = "#7f7f7f", "0" = "#7f7f7f")

# Partiti Plot ------------------------------------------------------------
partiti@phylo <- midpoint.root(partiti@phylo)
partiti <- groupOTU(partiti, cls)
partititree <- ggtree(partiti, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,4) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
partititree
ggsave("treePlots/partiti.svg",partititree, device = "svg", dpi=300)
ggsave("treePlots/partiti.pdf",partititree, device = "pdf", dpi=300)

# Hepe-astro-NS ---------------------------------------------------------

astrohepe <- read.mrbayes("data/cons_trees/astrohepe.nex.con.tre")

# Hepe-astro-NS Colors ----------------------------------------------------

cls <- list(NODE = c('BeeP_49_2013_NODE_27_length_6327_cov_6_608480_3'),
            Alphatetraviridae = c('EU345431.1_1','KX423453.1_1','NC_001981.1_1','NC_005898.1_1','U18246.1_1','AY594352.1_1'),
            Astroviridae = c('KX907135.1_1','NC_032426.1_1','NC_040647.1_1'),
            Hepeviridae = c('MF190001.1_1','NC_040710.1_1'),
            Bromoviridae = c('NC_009537.1_1'))
cols <- c("NODE" = "#ff7f0e","Alphatetraviridae"="#2ca02c","Astroviridae"="#2ca02c","Hepeviridae"="#2ca02c","Bromoviridae"="#2ca02c")


# Hepe-astro-NS Plot ------------------------------------------------------

astrohepe@phylo <- midpoint.root(astrohepe@phylo)
astrohepe <- groupOTU(astrohepe, cls)
astrohepetree <- ggtree(astrohepe, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,4) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
astrohepetree
ggsave("treePlots/astrohepe.svg",astrohepetree, device = "svg", dpi=300)
ggsave("treePlots/astrohepe.pdf",astrohepetree, device = "pdf", dpi=300)


# SinaiS-AA ---------------------------------------------------------------

sinais <- read.mrbayes("data/cons_trees/sinais-aa.nex.con.tre")


# SinaiS-AA Colors --------------------------------------------------------
cls <- list(NODE = c('BP11_NODE_2_length_5933_cov_8664_753785_2', 
                     'BP35_NODE_8_length_6159_cov_3073_471514_3',
                     'BeeP_34_2013_NODE_11_length_5700_cov_81_138894_3',
                     'BeeP_34_2013_NODE_8_length_5876_cov_117_977410_3'),
            SRA = c('SRR1239309.Contig_10906_56.5422_length_5835_3', 'SRR1239309.Contig_3402_33.8366_length_5852_3', 'SRR1239310.Contig_18500_86.3069_length_5976_3', 'SRR3927501.Contig_20_4621.33_length_5753_3', 'SRR5109829.Contig_5052_19.3642_length_5832_3', 'SRR5117449.Contig_4792_738.395_length_5543_3', 'SRR6031640.Contig_16980_3187.18_length_5924_3', 'SRR6031648.Contig_11953_83.802_length_5935_3', 'SRR6833958.Contig_36829_28.5413_length_5784_3', 'SRR806508.Contig_22923_17.5662_length_5787_3', 'SRR806508.Contig_4561_28.5866_length_5893_3', 'SRR806550.Contig_15712_2565.35_length_5987_3'),
            Sinai1 = c('HQ871931.2_2','KY465697.1_2','KY465698.1_2','KY465699.1_2','KY465700.1_2','KY465701.1_2','KY465702.1_2','KY465703.1_2','KY465704.1_2','KY465705.1_2','LR596015.1_2','NC_035466.1_2'),
            Sinai2 = c('HQ888865.2_2','KY465706.1_2','KY465707.1_2','KY465708.1_2','KY465709.1_2','KY465710.1_2','KY465711.1_2','KY465712.1_2','KY465713.1_2','LR655824.1_2','NC_035467.1_2'),
            Sinai = c('KM886902.1_2','KM886903.1_2','KM886904.1_2','KM886905.1_2','KX883223.1_2','KY465714.1_2','KY465715.1_2','KY465716.1_1','KY465717.1_2','KY465719.1_1','KY465720.1_2','MG918125.1_2','NC_032433.1_2'),
            SinaiTO = c('KY354241.1_2'),
            SinaiNE = c('KY354242.1_2','NC_035113.1_2'),
            SinaiSA1 = c('KY354243.1_2','NC_035111.1_2'),
            SinaiSA2 = c('KY354244.1_2', 'NC_035112.1_2'),
            Sinai3 = c('MH267699.1_2','MH267700.1_2'))

cols <- c("0" ="#7f7f7f", "NODE" = "#ff7f0e","SRA" = "#1f77b4","Sinai1" = "#2ca02c",
          "Sinai2" = "#2ca02c","Sinai" = "#2ca02c","SinaiTO" = "#2ca02c",
          "SinaiNE" = "#2ca02c","SinaiSA1" = "#2ca02c","SinaiSA2" = "#2ca02c",
          "Sinai3" = "#2ca02c")

# SinaiS-AA Plot ----------------------------------------------------------
sinais@phylo <- midpoint.root(sinais@phylo)
sinais <- groupOTU(sinais, cls)
sinaistree <- ggtree(sinais, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,0.5) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
sinaistree
ggsave("treePlots/sinais.svg",sinaistree, device = "svg", dpi=300)
ggsave("treePlots/sinais.pdf",sinaistree, device = "pdf", dpi=300)


# Toti --------------------------------------------------------------------

toti <- read.mrbayes("data/cons_trees/toti.nex.con.tre")


# Toti Colors -------------------------------------------------------------

cls <- list(NODE = c('BP23_NODE_31_length_5100_cov_72_148022_1_',
                     'BeeP_34_2013_NODE_9_length_5861_cov_32_563451_2_'),
            SRA = c('SRR1503019.Contig_22032_542.914_length_5984_1_',
                    'SRR6456322.Contig_9844_124.905_length_7917_1_'),
            victorivirus = c('NC_001963.1_2_',
                             'NC_001964.1_2_',
                             'NC_003607.1_2_',
                             'NC_003876.1_2_',
                             'NC_005074.1_2_',
                             'NC_006367.1_3_',
                             'NC_007523.1_2_',
                             'NC_010246.1_3_',
                             'NC_014823.1_2_',
                             'NC_020997.1_2_',
                             'NC_021565.1_2_',
                             'NC_023547.1_3_',
                             'NC_024151.1_2_',
                             'NC_025366.1_2_',
                             'NC_026140.1_2_',
                             'NC_027209.1_3_',
                             'NC_028477.1_2_',
                             'NC_030224.1_2_',
                             'NC_030392.1_2_',
                             'NC_030867.1_2_',
                             'NC_038928.1_2_',
                             'NC_038929.1_2_',
                             'NC_038930.1_2_',
                             'NC_040653.1_3_',
                             'NC_040793.1_2_'),
            toti = c('NC_009224.1_2_'),
            unclas = c('NC_025218.1_2_',
                       'NC_028948.1_2_',
                       'NC_029989.1_2_',
                       'NC_032424.1_1_',
                       'NC_032806.1_2_',
                       'NC_032819.1.1_',
                       'NC_032851.1_1_',
                       'NC_032931.1.1_',
                       'NC_032948.1_1_',
                       'NC_040530.1_2_'))

cols <- c("0"="#7f7f7f","NODE" = "#ff7f0e","SRA" = "#1f77b4","victorivirus" = "#2ca02c","toti" = "#2ca02c","unclas" = "#7f7f7f")



# Toti Plot ---------------------------------------------------------------
toti@phylo <- midpoint.root(toti@phylo)
toti <- groupOTU(toti, cls)
totitree <- ggtree(toti, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,8) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
totitree
ggsave("treePlots/toti.svg",totitree, device = "svg", dpi=300)
ggsave("treePlots/toti.pdf",totitree, device = "pdf", dpi=300)



# Tymo --------------------------------------------------------------------

tymo <- read.mrbayes("data/cons_trees/tymo.nex.con.tre")


# Tymo Colors -------------------------------------------------------------
cls <- list(NODE = c('BeeP_37_2013_NODE_30_length_5285_cov_819_119432_2_',
                     'BeeP_43_2013_NODE_1_length_6265_cov_122_817873_1_',
                     'BeeP_49_2013_NODE_19_length_8271_cov_90_694044_4_',
                     'BeeP_49_2013_NODE_37_length_5188_cov_28_459010_1_'),
            SRA = c('SRR1015503_Contig_4257_2592_36_length_5110_3_',
                    'SRR1015531_Contig_28871_50_2275_length_5084_1_',
                    'SRR1503108_Contig_2059_227_433_length_6823_1_',
                    'SRR5117444_Contig_2210_2_85334_length_5011_1_'),
            tymoviridae = c('NC_001480_1_1_',
                            'NC_001746_1_1_',
                            'NC_001793_1_1_',
                            'NC_002588_1_1_',
                            'NC_002786_1_1_',
                            'NC_003634_1_1_',
                            'NC_006950_1_1_',
                            'NC_007609_1_1_',
                            'NC_009532_1_1_',
                            'NC_015522_1_1_',
                            'NC_020470_1_1_',
                            'NC_020471_1_1_',
                            'NC_021851_1_1_',
                            'NC_027631_1_1_',
                            'NC_029063_1_1_',
                            'NC_031692_1_1_',
                            'NC_034205_1_1_',
                            'NC_038328_1_1_',
                            'NC_040565_1_1_'),
            Beta = c('NC_001948_1_1_',
                     'NC_003462_2_1_',
                     'NC_003877_1_1_',
                     'NC_005343_1_1_',
                     'NC_008020_1_1_',
                     'NC_009087_2_1_',
                     'NC_009383_1_1_',
                     'NC_009892_1_1_',
                     'NC_009991_1_1_',
                     'NC_011525_1_1_',
                     'NC_012038_1_1_',
                     'NC_013527_1_1_',
                     'NC_014821_1_1_',
                     'NC_017859_1_1_',
                     'NC_018175_1_1_',
                     'NC_018714_1_1_',
                     'NC_023892_1_1_',
                     'NC_025388_1_1_',
                     'NC_026616_2_1_',
                     'NC_027527_1_1_',
                     'NC_028868_1_1_',
                     'NC_028975_1_1_',
                     'NC_029085_1_1_',
                     'NC_029086_1_1_',
                     'NC_029088_1_1_',
                     'NC_031089_1_1_',
                     'NC_035203_1_1_',
                     'NC_038325_1_1_',
                     'NC_038966_1_1_',
                     'NC_040643_1_1_'),
            Alpha = c('NC_003400_1_1_',
                      'NC_028649_1_1_'))

cols <- c("NODE" = "#ff7f0e","SRA" = "#1f77b4","tymoviridae" = "#2ca02c","Beta" = "#2ca02c","Alpha" = "#2ca02c")


# Tymo Plot ---------------------------------------------------------------
tymo@phylo <- midpoint.root(tymo@phylo)
tymo <- groupOTU(tymo, cls)
tymotree <- ggtree(tymo, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,4) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
tymotree
ggsave("treePlots/tymo.svg",tymotree, device = "svg", dpi=300)
ggsave("treePlots/tymo.pdf",tymotree, device = "pdf", dpi=300)



# SinaiS_Reco -------------------------------------------------------------

sinaisReco <- read.mrbayes("data/cons_trees/sinais-reco.nex.con.tre")


# sinaiS_Reco Colors ------------------------------------------------------
cls <- list(NODE = c('BP11_NODE_2_length_5933_cov_8664_753785_3', 
                     'BP35_NODE_8_length_6159_cov_3073_471514_2',
                     'BeeP_34_2013_NODE_11_length_5700_cov_81_138894_2',
                     'BeeP_34_2013_NODE_8_length_5876_cov_117_977410_2'),
            SRA = c('SRR1239309.Contig_10906_56.5422_length_5835_2', 'SRR1239309.Contig_3402_33.8366_length_5852_2',
                    'SRR1239310.Contig_18500_86.3069_length_5976_2', 'SRR3927501.Contig_20_4621.33_length_5753_2', 'SRR5109829.Contig_5052_19.3642_length_5832_2',
                    'SRR5117449.Contig_4792_738.395_length_5543_2', 'SRR6031640.Contig_16980_3187.1', 'SRR6031648.Contig_11953_83.802_length_5935_2',
                    'SRR6833958.Contig_36829_28.5413_length_5784_2', 'SRR806508.Contig_22923_17.5662_length_5787_2', 'SRR806508.Contig_45.18.5866_length_5893_2', 'SRR806550.Contig_15712_2565.35_length_5987_2'),
            Sinai1 = c('HQ871931.2','KY465697.1','KY465698.1','KY465699.1','KY465700.1','KY465701.1','KY465702.1','KY465703.1','KY465704.1','KY465705.1','LR596015.1','NC_035466.1'),
            Sinai2 = c('HQ888865.2','KY465706.1','KY465707.1','KY465708.1','KY465709.1','KY465710.1','KY465711.1','KY465712.1','KY465713.1','LR655824.1','NC_035467.1'),
            Sinai = c('KM886902.1','KM886903.1','KM886904.1','KM886905.1','KX883223.1','KY465714.1','KY465715.1','KY465716.1_1','KY465717.1','KY465719.1_1','KY465720.1','MG918125.1','NC_032433.1'),
            SinaiTO = c('KY354241.1'),
            SinaiNE = c('KY354242.1','NC_035113.1'),
            SinaiSA1 = c('KY354243.1','NC_035111.1'),
            SinaiSA2 = c('KY354244.1', 'NC_035112.1'),
            Sinai3 = c('MH267699.1','MH267700.1'))

cols <- c("0" ="#7f7f7f", "NODE" = "#ff7f0e","SRA" = "#1f77b4","Sinai1" = "#2ca02c",
          "Sinai2" = "#2ca02c","Sinai" = "#2ca02c","SinaiTO" = "#2ca02c",
          "SinaiNE" = "#2ca02c","SinaiSA1" = "#2ca02c","SinaiSA2" = "#2ca02c",
          "Sinai3" = "#2ca02c")




# sinaiS-Reco Plot ------------------------------------------------------
sinaisReco@phylo <- midpoint.root(sinaisReco@phylo)
sinaisReco <- groupOTU(sinaisReco, cls)
sinaisRecotree <- ggtree(sinaisReco, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,2) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
sinaisRecotree
ggsave("treePlots/sinaisReco.svg",sinaisRecotree, device = "svg", dpi=300)
ggsave("treePlots/sinaisReco.pdf",sinaisRecotree, device = "pdf", dpi=300)

# Picorna -----------------------------------------------------------------

picorna <- read.mrbayes("data/cons_trees/picorna.nex.con.tre")


# Picorna Colors ----------------------------------------------------------
#I'm aware this list is too long. However Since I had to rerun trees (putting some stricter inclusion criteria), I cannot be bothered to slim it down.
cls <- list(NODE = c('BP25_NODE_4_length_10404_cov_23.1',
                     'BP36_NODE_2_length_9548_cov_2100_318023_1',
                     'BeeP_11_2013_NODE_163_length_6119_cov_25_630089_1',
                     'BeeP_11_2013_NODE_89_length_9828_cov_350_101733_3',
                     'BeeP_18_2013_NODE.1',
                     'BeeP_35_2013_NODE_.1',
                     'BeeP_35_2013_NODE_60_length_7079_cov_154_4498.1',
                     'BeeP_37_2013_NODE_19_length_6519_cov_37_551382_1'),
            SRA = c('DRR028891.Contig_13953_871.306_length_9532_2',
                    'DRR029070.Contig_14978_423.052_length_5044_1',
                    'DRR029079.Contig_25205_635.734_length_8400_2',
                    'DRR029088.Contig_20659_300.503_length_8750_1',
                    'DRR029097.Contig_45.1',
                    'DRR029626.Contig_27885_122.1.1',
                    'DRR042109.Contig_26746_1597.54_length_8949_1',
                    'DRR042118.Contig_27043_60.265_length_5744_2',
                    'DRR042121.Contig_4205_38.5737_length_6654_1',
                    'DRR042122.Contig_16826_595.626_length_11172_1',
                    'DRR042169.Contig_19034_73266.8_length_9854_1',
                    'DRR042170.Contig_11.1',
                    'DRR129220.Contig_35_1876.38_length_101.1',
                    'ERR1354112.Contig_47798_38.1082_length_10402_1',
                    'ERR1992770.Contig_8425_1079.4_length_5447_1',
                    'ERR528752.Contig_18762_33212.8_length_10106_1',
                    'ERR528754.Contig_3476_73310.3_length_8795_1',
                    'ERR528754.Contig_5610_8040..1',
                    'ERR968678.Contig_2754_32479.1',
                    'SRR1239302.Contig_22774_27701.9_length_7439_1',
                    'SRR1239305.Contig_6874_4817.18_length_9566_2',
                    'SRR1255009.Contig_29506_763.766_length_7867_1',
                    'SRR1255068.Contig_50.1',
                    'SRR1255260.Contig_29135_35127_length_10137_1',
                    'SRR1284947.Contig_37732_40.3568_length_5869_1',
                    'SRR1380984.Contig_11506_100.887_length_6036_1',
                    'SRR1380984.Contig_325_42.4399_length_8247_1',
                    'SRR1502929.Contig_35875_60.96.1',
                    'SRR1503108.Contig_34254_77.5805_length_8369_1',
                    'SRR1503154.Contig_27032_2475.67_length_8873_1',
                    'SRR1503154.Contig_3770_209.723_length_8994_2',
                    'SRR1565472.Contig_10105_23.72_length_6120_2',
                    'SRR1653580.Contig_5344_438778_length_8759_1',
                    'SRR1662295.Contig_1429_18586.2_length_7346_1',
                    'SRR1790682.Contig_13802_38819.9_length_7418_2',
                    'SRR1803085.Contig_16903_23.9063_length_6082_2',
                    'SRR1812779.Contig_15523_639.454_length_5837_1',
                    'SRR1864696.Contig_16496_10.1506_length_5649_1',
                    'SRR1948236.Contig_94.1',
                    'SRR2032844.Contig_14524_34.0605_length_8745_1',
                    'SRR2032844.Contig_29705_493.535_length_9462_2',
                    'SRR2032857.Contig_75.1',
                    'SRR2034262.Contig_4.1',
                    'SRR2080389.Contig_46465_452.775_length_6083_2',
                    'SRR2080396.Contig_81067_13.6737_length_5174_1',
                    'SRR2146250.Contig_21.1',
                    'SRR2146251.Contig_4574_7735.29_length_10115_1',
                    'SRR2177525.Contig_24394_5355.4_length_108.1',
                    'SRR2177537.Contig_25667_9240.42_length_6674_1',
                    'SRR2396634.Contig_5115_45463.4_length_9536_1',
                    'SRR3156045.Contig_17280_34942.2_length_5069_1',
                    'SRR3156047.Contig_21526_11327.6_length_8432_1',
                    'SRR4013500.Contig_35935_1395.48_length_10118_1',
                    'SRR4013516.Contig_14887_4563.27_length_10117_2',
                    'SRR4045764.Contig_28287_21113.9_length_5642_1',
                    'SRR4045850.Contig_190.1',
                    'SRR4084087.Contig_11614_8927.32_length_9557_2',
                    'SRR4250143.Contig_29674_17113.2_length_10102_1',
                    'SRR4250176.Contig_20214_11769.5_length_7583_1',
                    'SRR5109831.Contig_4302_97.6852_length_9578_1',
                    'SRR5109831.Contig_6219_31.8886_length_5533_1',
                    'SRR5116308.Contig_1706_32956.8_length_8385_3',
                    'SRR5117442.Contig_3275_5744.95_length_8643_1',
                    'SRR5117443.Contig_1594_129.633_length_6488_1',
                    'SRR5117443.Contig_4282_118656_length_6666_1',
                    'SRR5117445.Contig_29.1',
                    'SRR5117445.Contig_5646_105.828_length_60.1',
                    'SRR5117445.Contig_7392_2145.29_length_6380_1',
                    'SRR5117445.Contig_7393_2354.5_length_6379_1',
                    'SRR5117445.Contig_7394_1873..1',
                    'SRR5117445.Contig_7595_253.06_length_7010_1_3',
                    'SRR5117445.Contig_8800_43.8666_length_8985_3',
                    'SRR5117446.Contig_3180_534.627_length_9796_1',
                    'SRR5117446.Contig_5028_1577.23_length_9035_1',
                    'SRR5117446.Contig_5480_169.455_length_8182_1',
                    'SRR5117447.Contig_10204_77730.4_length_9212_2',
                    'SRR5117447.Contig_4947_2739.16_length_88.1',
                    'SRR5117447.Contig_9117_25.6389_length_8167_1',
                    'SRR5117449.Contig_6545_449886_length_7855_2',
                    'SRR5117450.Contig_1017_142.726_length_6550_1',
                    'SRR5117450.Contig_2957_74.6619_length_9162_2',
                    'SRR5136449.Contig_93.1',
                    'SRR5387738.Contig_3705_102404_length_5494_2',
                    'SRR5907694.Contig_18862_127.622_length_6708_1',
                    'SRR5907701.Contig_9257_31.1375_length_8376_1',
                    'SRR6001380.Contig_28689_35.9473_length_8774_1',
                    'SRR6047315.Contig_12488_10784.6_length_5117_1',
                    'SRR6788737.Contig_22243_290.138_length_5574_1',
                    'SRR6788742.Contig_21289_172.807_length_9627_1',
                    'SRR6806695.Contig_18575_113.485_length_6775_1',
                    'SRR6806700.Contig_16879_1139_length_5280_1',
                    'SRR6806702.Contig_8773_62.9543_length_7423_1',
                    'SRR6823683.Contig_17167_2220.97_length_8415_1',
                    'SRR6823685.Contig_1826_2787.67_length_84.1',
                    'SRR6833955.Contig_31686_58007.2_length_5475_1',
                    'SRR6833962.Contig_32592_101504_length_10143_1',
                    'SRR7192242.Contig_33466_34.1313_length_5260_1',
                    'SRR806508.Contig_1676_28.6109_length_8157_1',
                    'SRR806709.Contig_20209_5699.86_length_8808_1',
                    'SRR806709.Contig_26036_19626.1',
                    'SRR807359.Contig_2988_193.1'),
            ifla = c('AB070959.1',
                     'AF092924.1',
                     'AF469603.1',
                     'AJ489744.2_Deformed_Wing_Virus_gene_for_polyprotein',
                     'AY251269.2_Varroa_destructor_virus_1',
                     'AY292384.1',
                     'EU035616.1',
                     'GU109335.1',
                     'GU938761.1',
                     'HM067437.1',
                     'HM067438.1',
                     'HM162356.1',
                     'HQ322114.1',
                     'JF440525.1',
                     'JF440526.1',
                     'JQ390591.1',
                     'JQ390592.1',
                     'JQ413340.1',
                     'JX194121.1',
                     'JX270795.1',
                     'JX270796.1',
                     'JX270797.1',
                     'JX270798.1',
                     'JX270799.1',
                     'JX270800.1',
                     'JX878304.1',
                     'JX878305.1',
                     'KC007374.1',
                     'KC285046.1',
                     'KC786222.1',
                     'KC786223.1',
                     'KC786224.1',
                     'KF500002.1',
                     'KF960044.1',
                     'KJ000692.1',
                     'KJ437447.1',
                     'KJ629183.1',
                     'KJ716805.1',
                     'KJ716806.1',
                     'KJ959613.1',
                     'KJ959614.1',
                     'KM495267.1',
                     'KM884990.1',
                     'KM884991.1',
                     'KM884992.1',
                     'KM884993.1',
                     'KM884994.1',
                     'KM884995.1',
                     'KP296800.1',
                     'KP296801.1',
                     'KP296802.1',
                     'KP296803.1',
                     'KT004425.1',
                     'KU645789.1',
                     'KU847397.1',
                     'KX373900.1',
                     'KX580899.1',
                     'KX663835.1',
                     'KX668139.1',
                     'KX668140.1',
                     'KX668141.1',
                     'KX779454.1',
                     'KX783225.1',
                     'KX786162.1',
                     'KX819276.1',
                     'KY243931.1',
                     'KY273489.1',
                     'KY465671.1',
                     'KY465672.1',
                     'KY465673.1',
                     'KY465674.1',
                     'KY465675.1',
                     'KY465676.1',
                     'KY465677.1',
                     'KY465678.1',
                     'KY465679.1',
                     'KY774627.1',
                     'KY774628.1',
                     'KY887697.1',
                     'KY887698.1',
                     'KY887699.1',
                     'KY909333.1',
                     'MF036686.1',
                     'MF346349.1',
                     'MF623170.1',
                     'MF623172.1',
                     'MF770715.1',
                     'MG545286.1',
                     'MG545287.1',
                     'MH069503.1',
                     'MH069504.1',
                     'MH069505.1',
                     'MH069506.1',
                     'MH069507.1',
                     'MH107056.1',
                     'MH165180.1',
                     'MH267695.1',
                     'MH267696.1',
                     'MH267697.1',
                     'MH267698.1',
                     'MH509439.1',
                     'NC_002066.1',
                     'NC_004830.2_Deformed_wing_virus',
                     'NC_014137.1',
                     'NC_023022.1',
                     'NC_031338.1',
                     'NC_031749.1'),
            como = c('AB295643.1',
                     'AB456531.1',
                     'AF394606.1',
                     'AF394608.1',
                     'AY303786.1',
                     'AY744931.1',
                     'AY744932.1',
                     'D00915.1',
                     'EU450837.1',
                     'EU617326.1',
                     'FJ516745.1',
                     'FJ712026.1',
                     'GQ222381.1',
                     'GQ332372.1',
                     'GQ332373.1',
                     'GQ369526.1',
                     'GQ369527.1',
                     'GQ369528.1',
                     'GQ996948.1',
                     'GQ996951.1',
                     'GQ996952.1',
                     'GQ996953.1',
                     'GU562879.1',
                     'GU810903.1',
                     'GU968732.1',
                     'HE613269.1',
                     'HM032712.1',
                     'JF968120.1',
                     'JN391442.1',
                     'JQ975057.1',
                     'JX513889.1',
                     'JX513894.1',
                     'KC138732.1',
                     'KC900162.1',
                     'KP404602.1',
                     'KU522584.1',
                     'KX011072.1',
                     'KX011073.1',
                     'KX011074.1',
                     'KX011075.1',
                     'KX034840.1',
                     'KX034841.1',
                     'KX034842.1',
                     'KX034843.1',
                     'KX034844.1',
                     'KX034845.1',
                     'KX034846.1',
                     'KX034847.1',
                     'KX034848.1',
                     'KX034849.1',
                     'KX034850.1',
                     'KX034851.1',
                     'KX034852.1',
                     'KX034853.1',
                     'KX034854.1',
                     'KX034855.1',
                     'KX034856.1',
                     'KX034857.1',
                     'KX034858.1',
                     'KX034859.1',
                     'KX034860.1',
                     'KX034861.1',
                     'KX034862.1',
                     'KX034863.1',
                     'KX034864.1',
                     'KX034865.1',
                     'KX034866.1',
                     'KX034867.1',
                     'KX034868.1',
                     'KX034869.1',
                     'KX034870.1',
                     'KX034871.1',
                     'KX034872.1',
                     'KX034873.1',
                     'KX034874.1',
                     'KX034875.1',
                     'KX034876.1',
                     'KX034877.1',
                     'KX034878.1',
                     'KX034879.1',
                     'KX034880.1',
                     'KX034881.1',
                     'KX034882.1',
                     'KX034883.1',
                     'KX034884.1',
                     'KX034885.1',
                     'KX034886.1',
                     'KX034887.1',
                     'KX034888.1',
                     'KX034889.1',
                     'KX034890.1',
                     'KX034891.1',
                     'KX034892.1',
                     'KX034893.1',
                     'KX034894.1',
                     'KX034895.1',
                     'KX034896.1',
                     'KX034897.1',
                     'KX034898.1',
                     'KX034899.1',
                     'KX034900.1',
                     'KX034901.1',
                     'KX034902.1',
                     'KX034903.1',
                     'KX034904.1',
                     'KX034905.1',
                     'KY622124.1',
                     'KY701258.1',
                     'MF804979.1',
                     'MH053441.1',
                     'MH383240.1',
                     'MH383242.1',
                     'MH802014.1',
                     'MH802016.1',
                     'MH802018.1',
                     'MH802025.1',
                     'MH802027.1',
                     'MH802030.1',
                     'MH802032.1',
                     'NC_003496.1',
                     'NC_003741.1',
                     'NC_006057.1',
                     'NC_010709.1',
                     'NC_013218.1',
                     'NC_017939.1',
                     'NC_022004.1',
                     'NC_028139.1',
                     'U70866.1',
                     'X64886.1'),
            dicistro = c('AF150629.1',
                         'AF536531.1',
                         'AY275710.1',
                         'EF219380.1',
                         'EU218534.1',
                         'EU224279.1',
                         'EU224280.1',
                         'EU436423.1',
                         'EU436455.1',
                         'EU436456.1',
                         'HQ897161.1',
                         'JQ320375.1',
                         'JX045857.1',
                         'JX045858.1',
                         'JX480861.1',
                         'KC690268.1',
                         'KC690269.1',
                         'KC690270.1',
                         'KF500001.1',
                         'KF956377.1',
                         'KJ817182.1',
                         'KR021407.1',
                         'KX158871.1',
                         'KX421583.1',
                         'KX610809.1',
                         'KX830963.1',
                         'KX883690.1',
                         'KX883929.1',
                         'KX883977.1',
                         'KX884276.1',
                         'KY243933.1',
                         'KY465689.1',
                         'KY465690.1',
                         'KY465691.1',
                         'KY465692.1',
                         'KY465693.1',
                         'KY465694.1',
                         'KY465695.1',
                         'KY465696.1',
                         'LN907586.1',
                         'LN907587.1',
                         'LN907588.1',
                         'LN907589.1',
                         'MF189971.1',
                         'MF458892.1',
                         'MF458893.1',
                         'MF535297.1',
                         'MF795134.1',
                         'MF795135.1',
                         'MG599488.1',
                         'MH188004.1',
                         'NC_002548.1',
                         'NC_004365.1',
                         'NC_004807.1',
                         'NC_023021.1'),
            polycipi = c('EF428566.1',
                         'KX883910.1',
                         'MF041808.1',
                         'MF041809.1',
                         'MF041810.1',
                         'MG676340.1',
                         'MH213247.1',
                         'NC_032978.1',
                         'NC_035450.1',
                         'NC_035455.1',
                         'NC_035457.1',
                         'NC_039236.1'),
            marna = c('KF478836.2_Marine_RNA_virus_SF_3_polyprotein_gene',
                      'NC_043519.1'),
            picornavir = c('MN167485.1'),
            unclas = c('DQ321720.2_Nora_virus',
                       'GQ257737.1',
                       'JX220408.1',
                       'KF242510.1',
                       'KF242511.1',
                       'KP970076.1',
                       'KP970077.1',
                       'KP970078.1',
                       'KP970079.1',
                       'KP970080.1',
                       'KP970081.1',
                       'KP970082.1',
                       'KP970083.1',
                       'KP970084.1',
                       'KP970085.1',
                       'KP970086.1',
                       'KP970087.1',
                       'KP970088.1',
                       'KP970089.1',
                       'KP970090.1',
                       'KP970091.1',
                       'KP970092.1',
                       'KP970093.1',
                       'KP970094.1',
                       'KP970095.1',
                       'KP970096.1',
                       'KP970097.1',
                       'KP970098.1',
                       'KP970099.1',
                       'KP970100.1',
                       'KP970101.1',
                       'KP970102.1',
                       'KP970103.1',
                       'KP970104.1',
                       'KP970105.1',
                       'KP970106.1',
                       'KX883294.1',
                       'KX883645.1',
                       'KX883659.1',
                       'KX883694.1',
                       'KX883930.1',
                       'KX883975.1',
                       'KX884537.1',
                       'KY354240.1',
                       'MF287666.1',
                       'MF287668.1',
                       'MF593921.1',
                       'MF893254.1',
                       'MG995694.1',
                       'MG995696.1',
                       'MG995697.1',
                       'MG995699.1',
                       'MG995707.1',
                       'MG995709.1',
                       'MG995720.1',
                       'MG995723.1',
                       'MG995724.1',
                       'MG995731.1',
                       'MH384279.1',
                       'MH719200.1',
                       'MH899120.1',
                       'MK533158.1',
                       'MK643150.1',
                       'NC_007919.3_Nora_virus',
                       'NC_024487.1',
                       'NC_024488.1',
                       'NC_032904.1',
                       'NC_033094.1',
                       'NC_033455.1'))

cols <- c("NODE" = "#ff7f0e",
          "SRA" = "#1f77b4",
          "ifla" = "#2ca02c",
          "como" = "#2ca02c",
          "dicistro" = "#2ca02c",
          "polycipi" = "#2ca02c",
          "marna"="#2ca02c",
          "picornavir"="#2ca02c",
          "0" = "red",
          "unclas" = "#7f7f7f")
cols2 <- c("NODE" = "#ff7f0e",
           "SRA" = "#1f77b4",
           "ifla" = "red",
           "como" = "purple",
           "dicistro" = "blue",
           "polycipi" = "black",
           "marna"="darkgrey",
           "picornavir"="lightblue",
           "0" = "red",
           "unclas" = "#7f7f7f")
# Picorna Plot ------------------------------------------------------------
picorna@phylo <- midpoint.root(picorna@phylo)
picorna <- groupOTU(picorna, cls)
picornatree <- ggtree(picorna, aes(color=group), size=1) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = prob > 0.8), color='black',size=2) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,3) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
picornatree
ggsave("treePlots/picorna.svg",picornatree, device = "svg", dpi=300)
ggsave("treePlots/picorna.pdf",picornatree, device = "pdf", dpi=300)


# Rhabdo ------------------------------------------------------------------
rhabdo <- read.mrbayes("data/cons_trees/rhabdo.nex.con.tre")

# REST -----------------------------------------------------------





########################################### Picorna:
picorna <- read.beast("../data/MCC/Picorna.mcctree")


picorna <- groupOTU(picorna, cls)
cols <- c("NODE" = "#ff7f0e",
          "SRA" = "#1f77b4",
          "ifla" = "#2ca02c",
          "como" = "#2ca02c",
          "dicistro" = "#2ca02c",
          "polycipi" = "#2ca02c",
          "marna"="#2ca02c",
          "picornavir"="#2ca02c",
          "0" = "#7f7f7f",
          "unclas" = "#7f7f7f")
cols2 <- c("NODE" = "#ff7f0e",
          "SRA" = "#1f77b4",
          "ifla" = "red",
          "como" = "purple",
          "dicistro" = "blue",
          "polycipi" = "black",
          "marna"="darkgrey",
          "picornavir"="lightblue",
          "0" = "#7f7f7f",
          "unclas" = "#7f7f7f")

picornatree <- ggtree(picorna, aes(color=group),size=5) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) + 
  geom_nodepoint(aes(subset = posterior > 0.8), color='black',size=1) + 
  #geom_text(aes(label=node)) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,3) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
picornatree
#ggsave("Treeplots/picorna.svg",picornatree, device = "svg", dpi=380,width = 6.04, height =50,limitsize = FALSE)

#Collapse some nodes:
p2 <- collapse(picornatree, node=547) + geom_point2(aes(subset=(node==547)), shape=17, size=5, fill="brown")
p3 <- collapse(p2, node=515) + geom_point2(aes(subset=(node==515)), shape=17, size=5, fill="brown")
p4 <- collapse(p3, node=830) + geom_point2(aes(subset=(node==830)), shape=17, size=5, fill="brown")
p5 <- collapse(p4, node=744) + geom_point2(aes(subset=(node==744)), shape=17, size=5, fill="brown")
p6 <- collapse(p5, node=674) + geom_point2(aes(subset=(node==674)), shape=17, size=5, fill="brown")
p7 <- collapse(p6, node=659) + geom_point2(aes(subset=(node==387)), shape=17, size=5, fill="brown")
p8 <- collapse(p7, node=892) + geom_point2(aes(subset=(node==892)), shape=17, size=5, fill="brown")
p9 <- collapse(p8, node=940) + geom_point2(aes(subset=(node==940)), shape=17, size=5, fill="brown")
picornatree <- p9
picornatree
#ggsave("Treeplots/picorna.svg",p9, device = "svg", dpi=380,width = 6.04, height =50,limitsize = FALSE)
########################################### Rhabdo:

rhabdo <- read.beast("../data/MCC/Rhabdo.mcctree")
cls <- list(NODE = c('BeeP-34-2013_NODE_2_length_13379_cov_90_712825_1',
                     'BeeP-38-2013_NODE_4_length_14545_cov_165_264238_1'),
            SRA = c('DRR028884.Contig_7883_22.9921_length_11588_6',
                    'DRR029859.Contig_7759_128.983_length_12223_5',
                    'SRR1255153.Contig_26024_304.72_length_12652_7',
                    'SRR1503154.Contig_19807_185.122_length_12925_5',
                    'SRR2105774.Contig_25932_42.3786_length_11912_5',
                    'SRR2649096.Contig_3656_426.308_length_10338_5',
                    'SRR4013517.Contig_19821_37.8069_length_10983_2',
                    'SRR5868360.Contig_3779_58.7904_length_10073_5',
                    'SRR6819909.Contig_3811_325.126_length_12565_5',
                    'SRR7192237.Contig_7933_213.12_length_10496_7',
                    'SRR974922.Contig_11618_267.764_length_12229_5'),
            rhabdo = c('AB009601.1_Rabies_virus_mRNA_for_L_protein_RNA_dependent_RNA_polymerase',
                       'AB011257.1_Rice_yellow_stunt_virus_genomic_RNA',
                       'AB075039.1_Lettuce_big-vein_virus_LBVV-pol_gene_for_L_protein',
                       'AB244418.1_Orchid_fleck_virus_genomic_RNA',
                       'AB362483.1_Rabies_virus_viral_cRNA',
                       'AB516283.1_Rice_yellow_stunt_virus_viral_cRNA',
                       'AB516441.1_Orchid_fleck_virus_viral_cRNA',
                       'AB517660.1_Rabies_virus_viral_cRNA',
                       'AB735628.2_Persimmon_virus_A_viral_cRNA',
                       'AB981663.1_Rabies_virus_RNA',
                       'AB981664.1_Rabies_virus_RNA',
                       'AF418014.1_Australian_bat_lyssavirus',
                       'AF499686.2_Rabies_virus_strain_SRV9',
                       'AJ746199.1_Lettuce_necrotic_yellows_virus_L_gene_for_putative_polymerase',
                       'AJ810083.1_Chandipura_virus_L_gene_for_RNA-dependent_RNA_polymerase',
                       'AJ810084.2_Isfahan_virus_N_gene',
                       'AJ867584.2_Lettuce_necrotic_yellows_virus',
                       'AM689309.3_Drosophila_melanogaster_sigma_virus_AP30_N',
                       'AY618418.1_Maize_mosaic_virus',
                       'AY674964.1_Taro_vein_chlorosis_virus',
                       'DQ186554.1_Iranian_maize_mosaic_nucleorhabdovirus',
                       'EF157976.1_European_bat_lyssavirus_1_isolate_RV9',
                       'EF564174.1_Rabies_virus_strain_CTN181',
                       'EF614260.1_Irkut_virus',
                       'EF614261.1_Khujand_lyssavirus',
                       'EF687738.1_Lettuce_yellow_mottle_virus',
                       'EU259198.1_Lagos_bat_virus_isolate_KE131',
                       'EU293108.1_Lagos_bat_virus_isolate_0406SEN',
                       'EU293109.1_European_bat_lyssavirus_1_isolate_03002FRA',
                       'EU293110.1_Lagos_bat_virus_isolate_8619NGA',
                       'EU293111.1_Rabies_virus_isolate_8764THA',
                       'EU293112.1_European_bat_lyssavirus_1_isolate_8918FRA',
                       'EU293119.1_Duvenhage_virus_isolate_86132SA',
                       'EU293120.1_Duvenhage_virus_isolate_94286SA',
                       'EU293121.1_Rabies_virus_isolate_8743THA',
                       'EU373657.1_Cocal_virus_Indiana_2',
                       'EU623444.1_Duvenhage_virus_isolate_DUVVSA06',
                       'EU626551.1_European_bat_lyssavirus_1_isolate_08120FRA',
                       'EU626552.1_European_bat_lyssavirus_1_isolate_07240FRA',
                       'EU886633.1_Rabies_virus_collection-date_2004_from_Austria_nucleoprotein_N',
                       'FJ712194.1_Rabies_virus_isolate_D02',
                       'FJ712195.1_Rabies_virus_isolate_F02',
                       'FJ866835.1_Rabies_virus_strain_FJ008',
                       'FJ866836.1_Rabies_virus_strain_FJ009',
                       'FJ959397.1_Rabies_virus_strain_CTN-1',
                       'FJ985748.1_Moussa_virus_isolate_C23',
                       'FJ985749.1_Moussa_virus_isolate_D24',
                       'FR751552.4_Eggplant_mottled_dwarf_virus_G_gene_for_glycoprotein_partial_and_L_gene_for_polymerase_partial',
                       'GQ410980.1_Drosophila_affinis_sigma_virus_10_RNA-dependent_RNA_polymerase_L_gene',
                       'GU170201.1_Shimoni_bat_virus',
                       'GU170202.1_Lagos_bat_virus_isolate_KE576',
                       'GU190711.1_Chandipura_virus_isolate_CIN_0728',
                       'GU212856.1_Chandipura_virus_isolate_CIN_0451',
                       'GU212857.1_Chandipura_virus_isolate_CIN_0755',
                       'GU212858.1_Chandipura_virus_isolate_CIN_0327',
                       'GU345746.1_Rabies_virus_isolate_CQ92',
                       'GU345747.1_Rabies_virus_isolate_J',
                       'GU345748.1_Rabies_virus_isolate_SH06',
                       'GU734660.1_Potato_yellow_dwarf_virus_strain_SYDV',
                       'HM627182.1_Farmington_virus',
                       'HM627187.1_Chandipura_virus_Dak_AR_D_111125_nucleocapsid_protein',
                       'HM849039.1_Soybean_cyst_nematode_associated_northern_cereal_mosaic_virus_tegument_protein',
                       'HQ317918.1_Rabies_virus_strain_CTN-1-31',
                       'JF911700.1_Yug_Bogdanovac_virus',
                       'JN710440.1_Lettuce_big-vein_associated_virus_isolate_Ls302_L-protein_gene',
                       'JN786877.1_Rabies_virus_isolate_QS-05',
                       'JN786878.1_Rabies_virus_isolate_QS-BHK-P7',
                       'JN986749.1_Duvenhage_virus_isolate_NL07',
                       'JQ423952.1_Rabies_virus_isolate_BJ2011E',
                       'JQ595336.1_Rabies_virus_isolate_Fl148_RNA_polymerase_L_gene',
                       'JQ595337.1_Rabies_virus_isolate_WA1770_RNA_polymerase_L_gene',
                       'JQ595339.1_Rabies_virus_isolate_MS079_RNA_polymerase_L_gene',
                       'JQ595340.1_Rabies_virus_isolate_TX4904_RNA_polymerase_L_gene',
                       'JQ595342.1_Rabies_virus_isolate_TX4350_RNA_polymerase_L_gene',
                       'JQ595343.1_Rabies_virus_isolate_TX5168_RNA_polymerase_L_gene',
                       'JQ595347.1_Rabies_virus_isolate_TN132_RNA_polymerase_L_gene',
                       'JQ595348.1_Rabies_virus_isolate_ID7198_RNA_polymerase_L_gene',
                       'JQ595349.1_Rabies_virus_isolate_GA7034_RNA_polymerase_L_gene',
                       'JQ595353.1_Rabies_virus_isolate_AZ1838_RNA_polymerase_L_gene',
                       'JQ595354.1_Rabies_virus_isolate_FL1042_RNA_polymerase_L_gene',
                       'JQ595355.1_Rabies_virus_isolate_FL905_RNA_polymerase_L_gene',
                       'JQ595356.1_Rabies_virus_isolate_FL1165_RNA_polymerase_L_gene',
                       'JQ595357.1_Rabies_virus_isolate_FL769_RNA_polymerase_L_gene',
                       'JQ595358.1_Rabies_virus_isolate_ID7282_RNA_polymerase_L_gene',
                       'JQ595359.1_Rabies_virus_isolate_ID7261_RNA_polymerase_L_gene',
                       'JQ595360.1_Rabies_virus_isolate_ID7233_RNA_polymerase_L_gene',
                       'JQ595362.1_Rabies_virus_isolate_NJ1212_RNA_polymerase_L_gene',
                       'JQ595363.1_Rabies_virus_isolate_NJ1049_RNA_polymerase_L_gene',
                       'JQ595364.1_Rabies_virus_isolate_NJ2262_RNA_polymerase_L_gene',
                       'JQ595365.1_Rabies_virus_isolate_TN183_RNA_polymerase_L_gene',
                       'JQ595366.1_Rabies_virus_isolate_TX5751_RNA_polymerase_L_gene',
                       'JQ595367.1_Rabies_virus_isolate_TX5976_RNA_polymerase_L_gene',
                       'JQ595368.1_Rabies_virus_isolate_CA0253_RNA_polymerase_L_gene',
                       'JQ595369.1_Rabies_virus_isolate_WA1087_RNA_polymerase_L_gene',
                       'JQ595370.1_Rabies_virus_isolate_TX5742_RNA_polymerase_L_gene',
                       'JQ595371.1_Rabies_virus_isolate_WA173_RNA_polymerase_L_gene',
                       'JQ595372.1_Rabies_virus_isolate_WA2017_RNA_polymerase_L_gene',
                       'JQ595373.1_Rabies_virus_isolate_WA580_RNA_polymerase_L_gene',
                       'JQ595375.1_Rabies_virus_isolate_WA1185_RNA_polymerase_L_gene',
                       'JQ595376.1_Rabies_virus_isolate_VA1973_RNA_polymerase_L_gene',
                       'JQ595377.1_Rabies_virus_isolate_CA957_RNA_polymerase_L_gene',
                       'JQ595378.1_Rabies_virus_isolate_CA178_RNA_polymerase_L_gene',
                       'JQ595379.1_Rabies_virus_isolate_FL1078_RNA_polymerase_L_gene',
                       'JQ685943.1_Rabies_virus_isolate_A10-0511',
                       'JQ685975.1_Rabies_virus_isolate_MEXSK3636',
                       'JQ730682.1_Rabies_virus_strain_CYN1009D',
                       'JQ944705.1_Rabies_virus_isolate_1350KRA',
                       'JQ970481.1_Rabies_virus_strain_CJS0621D',
                       'JQ970482.1_Rabies_virus_strain_CJS0636D',
                       'JQ970483.1_Rabies_virus_strain_CJS0848D',
                       'JQ970486.1_Rabies_virus_strain_CSD0708D',
                       'JX129232.1_European_bat_lyssavirus_2_isolate_3278_09',
                       'JX129233.1_European_bat_lyssavirus_2_isolate_1985',
                       'JX403908.1_Drosophila_melanogaster_sigma_virus_strain_G0A_38',
                       'JX403909.1_Drosophila_melanogaster_sigma_virus_strain_G0A_239',
                       'JX403910.1_Drosophila_melanogaster_sigma_virus_strain_G0A_242',
                       'JX403911.1_Drosophila_melanogaster_sigma_virus_strain_G0A_325',
                       'JX403912.1_Drosophila_melanogaster_sigma_virus_strain_G0A_415',
                       'JX403913.1_Drosophila_melanogaster_sigma_virus_strain_G0A_440',
                       'JX403914.1_Drosophila_melanogaster_sigma_virus_strain_G0A_438',
                       'JX403915.1_Drosophila_melanogaster_sigma_virus_strain_G0B_38',
                       'JX403916.1_Drosophila_melanogaster_sigma_virus_strain_G0B_239',
                       'JX403918.1_Drosophila_melanogaster_sigma_virus_strain_G0B_325',
                       'JX403919.1_Drosophila_melanogaster_sigma_virus_strain_G0B_415',
                       'JX403920.1_Drosophila_melanogaster_sigma_virus_strain_G0B_438',
                       'JX403922.1_Drosophila_melanogaster_sigma_virus_strain_H_38',
                       'JX403923.1_Drosophila_melanogaster_sigma_virus_strain_H_239',
                       'JX403925.1_Drosophila_melanogaster_sigma_virus_strain_H_325',
                       'JX403929.1_Drosophila_melanogaster_sigma_virus_strain_L_38',
                       'JX403930.1_Drosophila_melanogaster_sigma_virus_strain_L_239',
                       'JX403931.1_Drosophila_melanogaster_sigma_virus_strain_L_242',
                       'JX403932.1_Drosophila_melanogaster_sigma_virus_strain_L_325',
                       'JX403934.1_Drosophila_melanogaster_sigma_virus_strain_L_440',
                       'JX442979.1_Irkut_virus_isolate_IRKV-THChina12',
                       'JX569193.2_American_bat_vesiculovirus_TFFN-2013_isolate_liver2008',
                       'JX901139.1_Lagos_bat_virus_isolate_KE476',
                       'K02378.1_Vesicular_stomatitis_Indiana_virus_polymerase_L_gene',
                       'KC193267.1_Rabies_virus_isolate_CNM1101C',
                       'KC252633.1_Rabies_virus_isolate_CNM1103C',
                       'KC252634.1_Rabies_virus_isolate_CNM1104D',
                       'KC595280.1_Rabies_virus_isolate_RusLipetsk8052f_2011',
                       'KC595281.1_Rabies_virus_isolate_RusLipetsk8053c_2011',
                       'KC595282.1_Rabies_virus_isolate_RusLipetsk8054f_2011',
                       'KC595283.1_Rabies_virus_isolate_RusLipetsk8057f_2011',
                       'KC602379.1_Farmington_virus_strain_CT_114',
                       'KC660078.1_Rabies_virus_strain_Fengtai',
                       'KC905081.1_Eggplant_mottled_dwarf_virus_isolate_Iran/SH-eg',
                       'KF155000.1_Rabies_virus_isolate_RV2516',
                       'KF155003.1_European_bat_lyssavirus_1_isolate_RV20',
                       'KF209276.1_Citrus_leprosis_virus_nuclear_type_isolate_M2345_segment_RNA2',
                       'KF410949.2_Eggplant_mottled_dwarf_virus_RNA-directed_RNA_polymerase_L_mRNA',
                       'KF468772.1_Chandipura_virus_strain_CH112',
                       'KF468773.1_Chandipura_virus_strain_CH157',
                       'KF468774.1_Chandipura_virus_strain_CH256',
                       'KF468775.1_Chandipura_virus_strain_I653514',
                       'KF935252.1_Vesicular_stomatitis_Indiana_virus_isolate_VSVdipP3_nonfunctional_L_protein_gene',
                       'KF947078.1_Spodoptera_frugiperda_rhabdovirus_isolate_Sf',
                       'KJ004416.1_Rabies_virus_isolate_JZ13-Lv',
                       'KJ082087.1_Eggplant_mottled_dwarf_virus_isolate_Agapanthus',
                       'KJ466147.1_Rabies_virus_strain_CTNCEC25',
                       'KJ564280.1_Rabies_virus_isolate_IMDRV-13',
                       'KJ685548.1_Australian_bat_lyssavirus_isolate_H2',
                       'KM016899.1_Rabies_virus_isolate_WQ14-RF',
                       'KM204983.1_Barur_virus_nucleoprotein',
                       'KM204990.1_Muir_Springs_virus_nucleoprotein',
                       'KM205003.1_Harlingen_virus_nucleoprotein',
                       'KM205010.1_Landjia_virus_nucleoprotein',
                       'KM205012.1_Rochambeau_virus_nucleoprotein',
                       'KM205018.1_Bahia_Grande_virus_nucleoprotein',
                       'KM594025.1_Rabies_virus_isolate_IP_5402/07',
                       'KM594039.1_Rabies_virus_isolate_IP_7841/09',
                       'KM817634.1_Sanxia_Water_Strider_Virus_5_strain_SXSSP11_nucleocapsid_N',
                       'KM817636.1_Shayang_Fly_Virus_3_strain_SYY1-1_ORF1_ORF1',
                       'KM817642.1_Tacheng_Tick_Virus_7_strain_TCRP-3_ORF1_ORF1',
                       'KM817643.1_Taishun_Tick_Virus_strain_BL198_nucleocapsid_N',
                       'KM817645.1_Wuhan_Ant_Virus_strain_WHMY02_RNA-dependent_RNA_polymerase_L_gene',
                       'KM817647.1_Wuhan_Fly_Virus_3_strain_SYY2-5_RNA-dependent_RNA_polymerase_L_gene',
                       'KM817649.1_Wuhan_House_Fly_Virus_2_strain_SYY4-5_ORF1_ORF1',
                       'KM817651.1_Wuhan_Insect_virus_5_strain_YCYC02_nucleocapsid_N',
                       'KM817659.1_Wuhan_Mosquito_Virus_9_strain_JX1-13_ORF1_ORF1',
                       'KM823531.1_Datura_yellow_vein_virus',
                       'KP241939.1_European_bat_lyssavirus_1_isolate_RV2416',
                       'KP735609.1_Diachasmimorpha_longicaudata_rhabdovirus_isolate_UGA',
                       'KP997032.1_Rabies_virus_isolate_PO-01_2014_Primorye',
                       'KR230089.1_Rabies_virus_isolate_SXBJ15',
                       'KR230090.1_Rabies_virus_isolate_SXYL15',
                       'KR534252.2_Rabies_lyssavirus_isolate_RV2899',
                       'KR822813.1_Drosophila_busckii_rhabdovirus_nucleocapsid_protein_N',
                       'KR822817.1_Drosophila_subobscura_rhabdovirus_RNA-dependent_RNA_polymerase_L_gene',
                       'KR822823.1_Drosophila_sturtevanti_rhabdovirus_1_RNA-dependent_RNA_polymerase_L_gene',
                       'KR822824.1_Drosophila_algonquin_sigmavirus_RNA-dependent_RNA_polymerase_L_gene',
                       'KR906741.1_Rabies_virus_isolate_RV2504.1_nucleoprotein_N',
                       'KR906746.1_Rabies_virus_isolate_RV2774.1_nucleoprotein_N',
                       'KR906747.1_Rabies_virus_isolate_RV2775.1_nucleoprotein_N',
                       'KR906753.1_Rabies_virus_isolate_RV2783.1_nucleoprotein_N',
                       'KR906757.1_Rabies_virus_isolate_RV2807.1_nucleoprotein_N',
                       'KR906779.1_Rabies_virus_isolate_RV3100.1_nucleoprotein_N',
                       'KR906790.1_Rabies_virus_isolate_RV3140.1_nucleoprotein_N',
                       'KT160285.1_Drosophila_melanogaster_sigmavirus_isolate_5f_nucleocapsid_protein_N_gene',
                       'KT160286.1_Drosophila_melanogaster_sigmavirus_isolate_7c_nucleocapsid_protein_N_gene',
                       'KT160287.1_Drosophila_melanogaster_sigmavirus_isolate_3h_nucleocapsid_protein_N_gene',
                       'KT336435.1_Rabies_virus_isolate_21467',
                       'KT336436.1_Rabies_virus_isolate_33512',
                       'KT336437.1_Rabies_virus_isolate_34312',
                       'KT868953.1_Australian_bat_lyssavirus_isolate_1510',
                       'KT868954.1_Australian_bat_lyssavirus_isolate_2549',
                       'KT868955.1_Australian_bat_lyssavirus_isolate_2566',
                       'KT868956.1_Australian_bat_lyssavirus_isolate_2430',
                       'KU754523.1_Withyham_virus_putative_L_protein_gene',
                       'KU946961.1_Rabies_virus_strain_CTN181-3',
                       'KX148102.1_Rabies_lyssavirus_isolate_91026MEX',
                       'KX148103.1_Rabies_lyssavirus_isolate_87021AFS',
                       'KX148108.1_Rabies_lyssavirus_isolate_11001NEP',
                       'KX148110.1_Rabies_lyssavirus_isolate_91014MEX_nucleoprotein',
                       'KX148111.1_Rabies_lyssavirus_isolate_91015MEX_nucleoprotein',
                       'KX148112.1_Rabies_lyssavirus_isolate_91010MEX_nucleoprotein',
                       'KX148145.1_Rabies_lyssavirus_isolate_86054YOU_nucleoprotein',
                       'KX148177.1_Rabies_lyssavirus_isolate_96317ISR_nucleoprotein',
                       'KX148178.1_Rabies_lyssavirus_isolate_96306ISR_nucleoprotein',
                       'KX148179.1_Rabies_lyssavirus_isolate_96314ISR_nucleoprotein',
                       'KX148180.1_Rabies_lyssavirus_isolate_93030ISR_nucleoprotein',
                       'KX148203.1_Rabies_lyssavirus_isolate_86031MOZ_nucleoprotein',
                       'KX148207.1_Rabies_lyssavirus_isolate_14015ITA_nucleoprotein',
                       'KX148209.1_Rabies_lyssavirus_isolate_04033MAD_nucleoprotein',
                       'KX148210.1_Rabies_lyssavirus_isolate_98002MAD_nucleoprotein',
                       'KX148211.1_Rabies_lyssavirus_isolate_86046MAD_nucleoprotein',
                       'KX148217.1_Rabies_lyssavirus_isolate_86123BRE_nucleoprotein',
                       'KX148218.1_Rabies_lyssavirus_isolate_14016BOT_nucleoprotein',
                       'KX148219.1_Rabies_lyssavirus_isolate_14017BOT_nucleoprotein',
                       'KX148220.1_Rabies_lyssavirus_isolate_15001AFS_nucleoprotein',
                       'KX148221.1_Rabies_lyssavirus_isolate_15003AFS_nucleoprotein',
                       'KX148222.1_Rabies_lyssavirus_isolate_15002AFS_nucleoprotein',
                       'KX148228.1_Rabies_lyssavirus_isolate_99001NEP_nucleoprotein',
                       'KX148245.1_Rabies_lyssavirus_isolate_09029NEP_nucleoprotein',
                       'KX148246.1_Rabies_lyssavirus_isolate_97002IND_nucleoprotein',
                       'KX148247.1_Rabies_lyssavirus_isolate_99009BIR_nucleoprotein',
                       'KX148248.1_Rabies_lyssavirus_isolate_99015BIR_nucleoprotein',
                       'KX148253.1_Rabies_lyssavirus_isolate_99012CBG_nucleoprotein',
                       'KX148255.1_Rabies_lyssavirus_isolate_99010LAO_nucleoprotein',
                       'KX148256.1_Rabies_lyssavirus_isolate_02004LAO_nucleoprotein',
                       'KX148257.1_Rabies_lyssavirus_isolate_02005LAO_nucleoprotein',
                       'KX148258.1_Rabies_lyssavirus_isolate_02003LAO_nucleoprotein',
                       'KX148259.1_Rabies_lyssavirus_isolate_94272PHI_nucleoprotein',
                       'KX148260.1_Rabies_lyssavirus_isolate_04030PHI_nucleoprotein',
                       'KX148261.1_Rabies_lyssavirus_isolate_94281PHI_nucleoprotein',
                       'KX148262.1_Rabies_lyssavirus_isolate_94278PHI_nucleoprotein',
                       'KX148263.1_Rabies_lyssavirus_isolate_94275PHI_nucleoprotein',
                       'KX148267.1_Rabies_lyssavirus_isolate_02050CHI_nucleoprotein',
                       'KX636164.1_Physostegia_chlorotic_mottle_virus_isolate_PV-1182',
                       'KX708503.1_Rabies_lyssavirus_strain_2835MxChisdog14',
                       'KX708504.1_Rabies_lyssavirus_strain_9084MxYucdog08',
                       'KX852388.1_Gata_virus_strain_M4_L_protein',
                       'KX982179.1_Citrus_leprosis_virus_N_strain_ibi1_segment_RNA2',
                       'KX982180.1_Citrus_leprosis_virus_N_strain_srq1_segment_RNA2',
                       'KX982181.1_Citrus_leprosis_virus_N_strain_sbs1_segment_RNA2',
                       'KY026415.1_Rabies_lyssavirus_isolate_NB.2014.0941',
                       'KY075646.1_Tomato_yellow_mottle-associated_virus',
                       'KY210245.1_Rabies_lyssavirus_isolate_rv3047',
                       'KY210252.1_Rabies_lyssavirus_isolate_rv3057',
                       'KY210263.1_Rabies_lyssavirus_isolate_rv3071',
                       'KY210309.1_Rabies_lyssavirus_isolate_rv3152',
                       'KY210311.1_Rabies_lyssavirus_isolate_rv3154',
                       'KY354230.1_Apis_rhabdovirus_1_isolate_AWD-1442_N_protein',
                       'KY354231.1_Apis_rhabdovirus_1_isolate_RI-49_N_protein',
                       'KY354232.1_Apis_rhabdovirus_1_isolate_T-23_N_protein',
                       'KY354233.1_Apis_rhabdovirus_2_isolate_RI-49_N_protein',
                       'KY354234.1_Apis_rhabdovirus_2_isolate_T-12_N_protein',
                       'KY549567.1_Potato_yellow_dwarf_nucleorhabdovirus_strain_CYDV-constricta',
                       'KY688151.1_European_bat_2_lyssavirus_isolate_FI85',
                       'KY700686.1_Citrus_chlorotic_spot_virus_strain_Trs1_segment_RNA2',
                       'KY706238.1_Physostegia_chlorotic_mottle_virus_isolate_PV-0831',
                       'KY751405.1_Citrus_leprosis_virus_N_isolate_MAS1_segment_RNA2',
                       'KY780299.1_Rabies_lyssavirus_isolate_CTMZZ11',
                       'KY810772.1_Cabbage_cytorhabdovirus_1_strain_FERA_050726',
                       'KY859866.1_Physostegia_chlorotic_mottle_virus_isolate_Hesse_HZ15-192',
                       'KY860593.1_Rabies_lyssavirus_isolate_RV2981',
                       'L32603.1_Sonchus_yellow_net_virus_complete_genome_6',
                       'LC222630.1_Orchid_fleck_dichorhavirus_viral_cRNA',
                       'LC270812.1_Menghai_rhabdovirus_RNA',
                       'LM645016.1_Rabies_virus_viral_cRNA_for_Nucleoprotein_N_gene',
                       'LN680656.1_Eggplant_mottled_dwarf_virus',
                       'LN849915.1_Lagos_bat_virus_complete_genome',
                       'LT839608.1_European_bat_lyssavirus_1_isolate_13424_genome_assembly',
                       'LT839609.1_European_bat_lyssavirus_1_isolate_20174_genome_assembly',
                       'LT839610.1_European_bat_lyssavirus_1_isolate_976_genome_assembly',
                       'LT839611.1_European_bat_lyssavirus_1_isolate_5782_genome_assembly',
                       'LT839612.1_European_bat_lyssavirus_1_isolate_5006_genome_assembly',
                       'LT839613.1_European_bat_lyssavirus_1_isolate_13027_genome_assembly',
                       'LT839614.1_European_bat_lyssavirus_1_isolate_5776_genome_assembly',
                       'LT839615.1_European_bat_lyssavirus_1_isolate_13454_genome_assembly',
                       'LT909527.1_Rabies_lyssavirus_isolate_13102_genome_assembly',
                       'LT909530.1_Rabies_lyssavirus_isolate_13251_genome_assembly',
                       'LT909537.1_Rabies_lyssavirus_isolate_12951_genome_assembly',
                       'LT909543.1_Rabies_lyssavirus_isolate_13249_genome_assembly',
                       'LT909544.1_Rabies_lyssavirus_isolate_13212_genome_assembly',
                       'LT909549.1_Rabies_lyssavirus_isolate_34873_genome_assembly',
                       'M87829.1_Sonchus_yellow_net_virus_L_protein_gene_mRNA',
                       'MF079256.1_Ekpoma_virus_isolate_EKV2-sh',
                       'MF102281.1_Maize_Iranian_mosaic_nucleorhabdovirus',
                       'MF114349.1_Apis_rhabdovirus_1_strain_Apis/USA/A',
                       'MF114350.1_Apis_rhabdovirus_1_strain_Bombus/USA/A',
                       'MF114351.1_Apis_rhabdovirus_1_strain_Varroa/Israel/B',
                       'MF187801.1_European_bat_1_lyssavirus_isolate_01018SLO_nucleoprotein',
                       'MF187802.1_European_bat_1_lyssavirus_isolate_02016DEN_nucleoprotein',
                       'MF187803.1_European_bat_1_lyssavirus_isolate_03002FRA_nucleoprotein',
                       'MF187804.1_European_bat_1_lyssavirus_isolate_04032FRA_nucleoprotein',
                       'MF187805.1_European_bat_1_lyssavirus_isolate_05001FRA_nucleoprotein',
                       'MF187806.1_European_bat_1_lyssavirus_isolate_06001FRA_nucleoprotein',
                       'MF187807.1_European_bat_1_lyssavirus_isolate_06002FRA_nucleoprotein',
                       'MF187808.1_European_bat_1_lyssavirus_isolate_07058FRA_nucleoprotein',
                       'MF187809.1_European_bat_1_lyssavirus_isolate_08120FRA_nucleoprotein',
                       'MF187810.1_European_bat_1_lyssavirus_isolate_08341FRA_nucleoprotein',
                       'MF187811.1_European_bat_1_lyssavirus_isolate_09034FRA_nucleoprotein',
                       'MF187812.1_European_bat_1_lyssavirus_isolate_78983_nucleoprotein',
                       'MF187813.1_European_bat_1_lyssavirus_isolate_107251_nucleoprotein',
                       'MF187814.1_European_bat_1_lyssavirus_isolate_113852_nucleoprotein',
                       'MF187815.1_European_bat_1_lyssavirus_isolate_116883_nucleoprotein',
                       'MF187816.1_European_bat_1_lyssavirus_isolate_120914_nucleoprotein',
                       'MF187817.1_European_bat_1_lyssavirus_isolate_121411_nucleoprotein',
                       'MF187818.1_European_bat_1_lyssavirus_isolate_121633_nucleoprotein',
                       'MF187819.1_European_bat_1_lyssavirus_isolate_122154_nucleoprotein',
                       'MF187820.1_European_bat_1_lyssavirus_isolate_122319_nucleoprotein',
                       'MF187821.1_European_bat_1_lyssavirus_isolate_122938_nucleoprotein',
                       'MF187822.1_European_bat_1_lyssavirus_isolate_123008_nucleoprotein',
                       'MF187823.1_European_bat_1_lyssavirus_isolate_123801_nucleoprotein',
                       'MF187824.1_European_bat_1_lyssavirus_isolate_124193_nucleoprotein',
                       'MF187825.1_European_bat_1_lyssavirus_isolate_124345_nucleoprotein',
                       'MF187826.1_European_bat_1_lyssavirus_isolate_124489_nucleoprotein',
                       'MF187827.1_European_bat_1_lyssavirus_isolate_126669_nucleoprotein',
                       'MF187828.1_European_bat_1_lyssavirus_isolate_127051_nucleoprotein',
                       'MF187829.1_European_bat_1_lyssavirus_isolate_127834_nucleoprotein',
                       'MF187830.1_European_bat_1_lyssavirus_isolate_127835_nucleoprotein',
                       'MF187831.1_European_bat_1_lyssavirus_isolate_128210_nucleoprotein',
                       'MF187832.1_European_bat_1_lyssavirus_isolate_128633_nucleoprotein',
                       'MF187833.1_European_bat_1_lyssavirus_isolate_128635_nucleoprotein',
                       'MF187834.1_European_bat_1_lyssavirus_isolate_128637_nucleoprotein',
                       'MF187835.1_European_bat_1_lyssavirus_isolate_128665_nucleoprotein',
                       'MF187836.1_European_bat_1_lyssavirus_isolate_128681_nucleoprotein',
                       'MF187837.1_European_bat_1_lyssavirus_isolate_128683_nucleoprotein',
                       'MF187838.1_European_bat_1_lyssavirus_isolate_128708_nucleoprotein',
                       'MF187839.1_European_bat_1_lyssavirus_isolate_128827_nucleoprotein',
                       'MF187840.1_European_bat_1_lyssavirus_isolate_129051_nucleoprotein',
                       'MF187841.1_European_bat_1_lyssavirus_isolate_129055_nucleoprotein',
                       'MF187842.1_European_bat_1_lyssavirus_isolate_129087_nucleoprotein',
                       'MF187843.1_European_bat_1_lyssavirus_isolate_129090_nucleoprotein',
                       'MF187844.1_European_bat_1_lyssavirus_isolate_129116_nucleoprotein',
                       'MF187845.1_European_bat_1_lyssavirus_isolate_129123_nucleoprotein',
                       'MF187846.1_European_bat_1_lyssavirus_isolate_129246_nucleoprotein',
                       'MF187847.1_European_bat_1_lyssavirus_isolate_129290_nucleoprotein',
                       'MF187848.1_European_bat_1_lyssavirus_isolate_129394_nucleoprotein',
                       'MF187849.1_European_bat_1_lyssavirus_isolate_129396_nucleoprotein',
                       'MF187850.1_European_bat_1_lyssavirus_isolate_129409_nucleoprotein',
                       'MF187851.1_European_bat_1_lyssavirus_isolate_129428_nucleoprotein',
                       'MF187852.1_European_bat_1_lyssavirus_isolate_130544_nucleoprotein',
                       'MF187853.1_European_bat_1_lyssavirus_isolate_130576_nucleoprotein',
                       'MF187854.1_European_bat_1_lyssavirus_isolate_130662_nucleoprotein',
                       'MF187855.1_European_bat_1_lyssavirus_isolate_130904_nucleoprotein',
                       'MF187856.1_European_bat_1_lyssavirus_isolate_131054_nucleoprotein',
                       'MF187857.1_European_bat_1_lyssavirus_isolate_15007FRA_nucleoprotein',
                       'MF187858.1_European_bat_1_lyssavirus_isolate_8615POL_nucleoprotein',
                       'MF187859.1_European_bat_1_lyssavirus_isolate_8918FRA_nucleoprotein',
                       'MF187860.1_European_bat_1_lyssavirus_isolate_9366HOL_nucleoprotein',
                       'MF187861.1_European_bat_1_lyssavirus_isolate_9367HOL_nucleoprotein',
                       'MF187862.1_European_bat_1_lyssavirus_isolate_9376HOL_nucleoprotein',
                       'MF187863.1_European_bat_1_lyssavirus_isolate_9377HOL_nucleoprotein',
                       'MF187864.1_European_bat_1_lyssavirus_isolate_9394POL_nucleoprotein',
                       'MF187865.1_European_bat_1_lyssavirus_isolate_9395GER_nucleoprotein',
                       'MF187866.1_European_bat_1_lyssavirus_isolate_9396GER_nucleoprotein',
                       'MF187867.1_European_bat_1_lyssavirus_isolate_9397RUS_nucleoprotein',
                       'MF187868.1_European_bat_1_lyssavirus_isolate_9399GER_nucleoprotein',
                       'MF187869.1_European_bat_1_lyssavirus_isolate_94113HOL_nucleoprotein',
                       'MF187870.1_European_bat_1_lyssavirus_isolate_94115HOL_nucleoprotein',
                       'MF187871.1_European_bat_1_lyssavirus_isolate_94116HOL_nucleoprotein',
                       'MF187872.1_European_bat_1_lyssavirus_isolate_94285SPA_nucleoprotein',
                       'MF187873.1_European_bat_1_lyssavirus_isolate_9436GER_nucleoprotein',
                       'MF187874.1_European_bat_1_lyssavirus_isolate_9438GER_nucleoprotein',
                       'MF187875.1_European_bat_1_lyssavirus_isolate_9440GER_nucleoprotein',
                       'MF187876.1_European_bat_1_lyssavirus_isolate_9443UKR_nucleoprotein',
                       'MF187877.1_European_bat_1_lyssavirus_isolate_9477GER_nucleoprotein',
                       'MF187878.1_European_bat_1_lyssavirus_isolate_9478HOL_nucleoprotein',
                       'MF187879.1_European_bat_1_lyssavirus_isolate_9480HOL_nucleoprotein',
                       'MF187880.1_European_bat_1_lyssavirus_isolate_9483SPA_nucleoprotein',
                       'MF197741.1_Rabies_virus_isolate_1311200108POL_complete_genome_5',
                       'MF197743.1_Rabies_virus_isolate_1359120810POL_complete_genome_5',
                       'MF197744.1_European_bat_1_lyssavirus_isolate_067N_complete_genome_5',
                       'MF197745.1_European_bat_1_lyssavirus_isolate_019N_complete_genome_5',
                       'MF360790.1_Blacklegged_tick_rhabdovirus-1_isolate_RTS95.16',
                       'MF360791.1_Dog_Tick_rhabdovirus-1_isolate_RTS107',
                       'MF536978.1_Spodoptera_frugiperda_rhabdovirus_isolate_Sf9',
                       'MF536979.1_Spodoptera_frugiperda_rhabdovirus_isolate_Sf21',
                       'MF543022.1_Black_currant_nucleorhabdovirus_1_isolate_Veloy',
                       'MF918568.1_Red_clover_varicosavirus_isolate_HZ2_segment_RNA1',
                       'MF962659.1_American_dog_tick_rhabdovirus-2_isolate_RTS-110_polymerase_gene',
                       'MF975531.1_Chimay_rhabdovirus_isolate_Chimay-1',
                       'MG201920.1_Rabies_lyssavirus_isolate_GXN119_nucleprotein_N',
                       'MG201922.1_Rabies_lyssavirus_isolate_GXBH2011_nucleoprotein_N',
                       'MG385079.1_Grenada_mosquito_rhabdovirus_1',
                       'MG458308.1_Rabies_lyssavirus_isolate_RV1009',
                       'MG458318.1_Rabies_lyssavirus_isolate_RV2481',
                       'MG552608.1_Vesicular_stomatitis_New_Jersey_virus_isolate_NJ0806VCB',
                       'MG562529.1_Rabies_lyssavirus_isolate_NY.2003.3040',
                       'MG562531.1_Rabies_lyssavirus_isolate_NY.2003.7914',
                       'MG562532.1_Rabies_lyssavirus_isolate_NY.2003.7917',
                       'MG562536.1_Rabies_lyssavirus_isolate_NY.2004.1705',
                       'MG562537.1_Rabies_lyssavirus_isolate_NY.2004.2620',
                       'MG562538.1_Rabies_lyssavirus_isolate_NY.2004.4138',
                       'MG562539.1_Rabies_lyssavirus_isolate_NY.2004.5339',
                       'MG562544.1_Rabies_lyssavirus_isolate_NY.2010.0775',
                       'MG562546.1_Rabies_lyssavirus_isolate_NY.2010.5433',
                       'MG562547.1_Rabies_lyssavirus_isolate_NY.2011.0005',
                       'MG562599.1_Rabies_lyssavirus_isolate_VT.2009.0353',
                       'MG562607.1_Rabies_lyssavirus_isolate_VT.2010.0003',
                       'MG600015.1_Fujian_dimarhabdovirus_strain_BHNC4885_nucleoprotein',
                       'MG604920.1_Wheat_yellow_striate_virus_isolate_SX-HC_nucleocapsid_protein',
                       'MG717932.1_Citrus_chlorotic_spot_virus_isolate_Trs2_segment_RNA2',
                       'MG948563.1_Alfalfa_associated_nucleorhabdovirus',
                       'MG983789.1_Clerodendrum_chlorotic_spot_virus_isolate_SBO1_segment_RNA2',
                       'MG983791.1_Clerodendrum_chlorotic_spot_virus_isolate_SPa1_segment_RNA2',
                       'MH102387.1_Rhinolophus_rhabdovirus_DPuer',
                       'MH129615.1_Strawberry_crinkle_cytorhabdovirus_isolate_A',
                       'MH129616.1_Strawberry_crinkle_cytorhabdovirus_isolate_B',
                       'MH213246.1_Linepithema_humile_rhabdo-like_virus_1_putative_capsid',
                       'MH267691.1_Apis_rhabdovirus_1_isolate_ARV-1_MR',
                       'MH267692.1_Apis_rhabdovirus_1_isolate_ARV-1_MS',
                       'MH323437.1_Zhuye_pepper_nucleorhabdovirus_isolate_ZPNu1',
                       'MH349091.1_Lettuce_big-vein_associated_varicosavirus_isolate_xm_L_protein_gene',
                       'MH384277.1_Drosophila_melanogaster_sigmavirus_strain_COFplus27050_nucleocapsid_protein',
                       'MH384306.1_Drosophila_melanogaster_sigmavirus_strain_MARminus23478_nucleocapsid_protein',
                       'MH595622.1_Citrus_chlorotic_spot_virus_isolate_Trs4_segment_RNA2_RNA-dependent_RNA_polymerase_ORF6_gene',
                       'MH620817.1_Lampyris_noctiluca_rhabdo-like_virus_1_isolate_17FIN6_segment_3',
                       'MH688523.1_Taishun_Tick_Virus_strain_17-L2_NP',
                       'MH707450.1_Grenada_mosquito_rhabdovirus_1_isolate_IVRI-2017',
                       'MH778545.1_Apple_rootstock_virus_A',
                       'MH926029.1_Spodoptera_frugiperda_rhabdovirus',
                       'MH926030.1_Spodoptera_frugiperda_rhabdovirus',
                       'MH926031.1_Spodoptera_frugiperda_rhabdovirus',
                       'MK063878.2_Morogoro_maize-associated_virus_isolate_16-0112_nucleocapsid_protein',
                       'MK063879.2_Morogoro_maize-associated_virus_isolate_16-0114_nucleocapsid_protein',
                       'MK112501.2_Morogoro_maize-associated_virus_isolate_16-0121_nucleocapsid_protein',
                       'MK159261.1_Strawberry_associated_virus_1',
                       'MK240091.1_Raspberry_vein_chlorosis_virus_isolate_Hutton_1',
                       'MK257717.1_Raspberry_vein_chlorosis_virus_isolate_Hutton_2',
                       'MK522805.1_Orchid_fleck_dichorhavirus_isolate_CL_segment_RNA2',
                       'MK522807.1_Orchid_fleck_dichorhavirus_isolate_Br_segment_RNA2',
                       'MK540667.1_Rabies_lyssavirus_isolate_NY.1990.2981R',
                       'MK540697.1_Rabies_lyssavirus_isolate_NY.2004.0057S',
                       'MK540698.1_Rabies_lyssavirus_isolate_NY.2004.0123R',
                       'MK540706.1_Rabies_lyssavirus_isolate_NY.2004.3099R',
                       'MK540707.1_Rabies_lyssavirus_isolate_NY.2004.3856R',
                       'MK540709.1_Rabies_lyssavirus_isolate_NY.2004.4735R',
                       'MK540712.1_Rabies_lyssavirus_isolate_NY.2004.6127S',
                       'MK540713.1_Rabies_lyssavirus_isolate_NY.2004.6604S',
                       'MK540714.1_Rabies_lyssavirus_isolate_NY.2004.6605R',
                       'MK540717.1_Rabies_lyssavirus_isolate_NY.2004.7502R',
                       'MK540720.1_Rabies_lyssavirus_isolate_NY.2004.8082S',
                       'MK540721.1_Rabies_lyssavirus_isolate_NY.2004.8258R',
                       'MK540724.1_Rabies_lyssavirus_isolate_NY.2010.1147R',
                       'MK540733.1_Rabies_lyssavirus_isolate_NY.2010.3151R',
                       'MK540735.1_Rabies_lyssavirus_isolate_NY.2010.3434R',
                       'MK540740.1_Rabies_lyssavirus_isolate_NY.2010.5963S',
                       'MK540742.1_Rabies_lyssavirus_isolate_NY.2011.0005S',
                       'MK540747.1_Rabies_lyssavirus_isolate_NY.2011.0141O',
                       'MK540752.1_Rabies_lyssavirus_isolate_NY.2011.0525R',
                       'MK540762.1_Rabies_lyssavirus_isolate_NY.2011.0948O',
                       'MK540767.1_Rabies_lyssavirus_isolate_NY.2011.1517R',
                       'MK540769.1_Rabies_lyssavirus_isolate_NY.2011.2074O',
                       'MK540771.1_Rabies_lyssavirus_isolate_NY.2011.2208R',
                       'MK540773.1_Rabies_lyssavirus_isolate_NY.2011.2320R',
                       'MK540774.1_Rabies_lyssavirus_isolate_NY.2011.2342O',
                       'MK540782.1_Rabies_lyssavirus_isolate_NY.2011.4696S',
                       'MK540785.1_Rabies_lyssavirus_isolate_NY.2011.5062R',
                       'MK540788.1_Rabies_lyssavirus_isolate_NY.2011.5446S',
                       'MK540793.1_Rabies_lyssavirus_isolate_NY.2011.6409S',
                       'MK540890.1_Rabies_lyssavirus_isolate_ON.1990.9340SSK',
                       'MK540903.1_Rabies_lyssavirus_isolate_ON.1993.0215RFX',
                       'MK567666.1_Rabies_lyssavirus_strain_GSQY-Dog-China-2013',
                       'MK577649.1_Rabies_lyssavirus_strain_GSJC-Sheep-China-2015',
                       'MK760677.1_Rabies_lyssavirus_isolate_87-395_nucleoprotein',
                       'MK760709.1_Rabies_lyssavirus_isolate_91-A-57_nucleoprotein',
                       'MK760727.1_Rabies_lyssavirus_isolate_92-343_nucleoprotein',
                       'MK760728.1_Rabies_lyssavirus_isolate_92-350_nucleoprotein',
                       'MK760739.1_Rabies_lyssavirus_isolate_92-75_nucleoprotein',
                       'MK760764.1_Rabies_lyssavirus_isolate_93-A-88_nucleoprotein',
                       'MK828539.1_Maize_mosaic_nucleorhabdovirus',
                       'MK948541.1_Physostegia_chlorotic_mottle_virus_isolate_JKI_ID_31401',
                       'MN013386.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQM-1-3',
                       'MN013387.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQF-1-5',
                       'MN013388.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQF-5',
                       'MN013389.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQM-1-4',
                       'MN013390.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQM-1-5',
                       'MN013391.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQM-5',
                       'MN013392.1_Guadeloupe_Culex_rhabdovirus_strain_2017-PB-CQM-1-1',
                       'MN013393.1_Guadeloupe_Culex_rhabdovirus_strain_2016-Ab-CQF',
                       'NC_001542.1_Rabies_virus',
                       'NC_001615.3_Sonchus_yellow_net_virus_complete_genome_6',
                       'NC_003243.1_Australian_bat_lyssavirus',
                       'NC_003746.1_Rice_yellow_stunt_virus',
                       'NC_005975.1_Maize_mosaic_virus',
                       'NC_006429.1_Mokola_virus',
                       'NC_006942.1_Taro_vein_chlorosis_virus',
                       'NC_007642.1_Lettuce_necrotic_yellows_virus',
                       'NC_008514.1_Siniperca_chuatsi_rhabdovirus',
                       'NC_009527.1_European_bat_lyssavirus_1',
                       'NC_009528.2_European_bat_lyssavirus_2_isolate_RV1333',
                       'NC_009609.1_Orchid_fleck_virus_genomic_RNA',
                       'NC_011532.1_Lettuce_yellow_mottle_virus',
                       'NC_011558.1_Lettuce_big-vein_associated_virus_RNA_1',
                       'NC_016136.1_Potato_yellow_dwarf_virus',
                       'NC_018381.2_Persimmon_virus_A_viral_cRNA',
                       'NC_020804.1_Tibrogargan_virus_strain_CS132',
                       'NC_020805.1_Chandipura_virus_isolate_CIN_0451',
                       'NC_020807.1_Lagos_bat_virus_isolate_0406SEN',
                       'NC_020808.1_Aravan_virus',
                       'NC_020810.1_Duvenhage_virus_isolate_86132SA',
                       'NC_022580.1_Drosophila_obscura_sigma_virus_10A',
                       'NC_022755.1_American_bat_vesiculovirus_TFFN-2013_isolate_liver2008',
                       'NC_024473.1_Vesicular_stomatitis_New_Jersey_virus_isolate_NJ1184HDB',
                       'NC_025251.1_Bokeloh_bat_lyssavirus_isolate_21961',
                       'NC_025253.1_Farmington_virus_strain_CT_114',
                       'NC_025359.1_Moussa_virus_isolate_C23',
                       'NC_025365.1_Shimoni_bat_virus',
                       'NC_025378.1_Yug_Bogdanovac_virus',
                       'NC_025382.1_Spodoptera_frugiperda_rhabdovirus_isolate_Sf',
                       'NC_025385.1_Khujand_lyssavirus',
                       'NC_025389.1_Eggplant_mottled_dwarf_virus_isolate_Agapanthus',
                       'NC_025391.1_Almpiwar_virus_isolate_MRM4059',
                       'NC_025397.1_Coastal_Plains_virus_strain_DPP53',
                       'NC_025399.1_Oak-Vale_virus_strain_CSIRO_1342',
                       'NC_025405.1_Niakha_virus_isolate_DakArD_88909',
                       'NC_025408.1_Lyssavirus_Ozernoe',
                       'NC_028230.1_Inhangapi_virus_strain_BEAR177325',
                       'NC_028231.1_Datura_yellow_vein_virus',
                       'NC_028232.1_Walkabout_Creek_virus_isolate_CS1056',
                       'NC_028234.1_Santa_barbara_virus_strain_AR775619',
                       'NC_028237.2_Alfalfa_dwarf_virus_isolate_Manfredi',
                       'NC_028241.1_Yata_virus_isolate_DakArB_2181',
                       'NC_028255.1_Cocal_virus_Indiana_2',
                       'NC_030451.1_Diachasmimorpha_longicaudata_rhabdovirus_isolate_UGA',
                       'NC_031079.1_Bole_Tick_Virus_2_strain_BL076_nucleocapsid_N',
                       'NC_031083.1_Huangpi_Tick_Virus_3_strain_H124-2_nucleocapsid_N',
                       'NC_031093.1_Sanxia_Water_Strider_Virus_5_strain_SXSSP11_nucleocapsid_N',
                       'NC_031227.1_Wuhan_Insect_virus_5_strain_YCYC02_nucleocapsid_N',
                       'NC_031236.1_Wuhan_Insect_virus_7_strain_WHYC02_nucleocapsid_N',
                       'NC_031272.1_Tacheng_Tick_Virus_7_strain_TCRP-3_ORF1_ORF1',
                       'NC_031273.1_Taishun_Tick_Virus_strain_BL198_nucleocapsid_N',
                       'NC_031276.1_Wuhan_Ant_Virus_strain_WHMY02_RNA-dependent_RNA_polymerase_L_gene',
                       'NC_031283.1_Wuhan_House_Fly_Virus_2_strain_SYY4-5_ORF1_ORF1',
                       'NC_031301.1_Wuhan_Louse_Fly_Virus_5_strain_BFJSC-5_nucleocapsid_N',
                       'NC_031303.1_Wuhan_Mosquito_Virus_9_strain_JX1-13_ORF1_ORF1',
                       'NC_031691.1_Gata_virus_strain_M4_L_protein',
                       'NC_031955.1_Lleida_bat_lyssavirus_isolate_RV3208',
                       'NC_031957.1_UNVERIFIED',
                       'NC_031988.1_Gannoruwa_bat_lyssavirus_isolate_RV3266',
                       'NC_033705.1_Xinzhou_nematode_virus_4_strain_XZSJSC65771_putative_nucleoprotein',
                       'NC_034240.1_Tomato_yellow_mottle-associated_virus',
                       'NC_034529.1_Sena_Madureira_virus_nucleoprotein',
                       'NC_034533.1_Landjia_virus_nucleoprotein',
                       'NC_034534.1_Rochambeau_virus_nucleoprotein',
                       'NC_034535.1_Barur_virus_nucleoprotein',
                       'NC_034545.1_Mount_Elgon_bat_virus_nucleoprotein',
                       'NC_034546.1_Sweetwater_Branch_virus_nucleoprotein',
                       'NC_034548.1_Oita_virus_nucleoprotein',
                       'NC_034549.1_Klamath_virus_nucleoprotein',
                       'NC_034550.1_Chaco_virus_nucleoprotein',
                       'NC_036390.1_Maize_Iranian_mosaic_nucleorhabdovirus',
                       'NC_038236.1_Vesicular_stomatitis_Indiana_virus_strain_98COE',
                       'NC_038278.1_Drosophila_affinis_sigmavirus_nucleocapsid_protein_N',
                       'NC_038280.1_Drosophila_immigrans_sigmavirus_strain_SCM45623_nucleocapsid_protein',
                       'NC_038281.1_Drosophila_melanogaster_sigma_virus_HAP23',
                       'NC_038282.1_Ekpoma_virus_1_isolate_EKV-1',
                       'NC_038283.1_Ekpoma_virus_2_isolate_EKV-2',
                       'NC_038285.1_Carajas_virus_nucleoprotein',
                       'NC_038287.1_Radi_virus_nucleoprotein',
                       'NC_038755.1_Coffee_ringspot_virus_strain_Lavras_segment_RNA2',
                       'NC_039020.1_Kanyawara_virus_isolate_MPK004_nucleoprotein',
                       'NC_039021.1_Beatrice_Hill_virus_isolate_CSIRO_25',
                       'NC_039206.1_Jurona_virus_nucleoprotein',
                       'NC_040602.1_Menghai_rhabdovirus_isolate_Menghai',
                       'NC_043649.1_Clerodendrum_chlorotic_spot_virus_isolate_Prb1_segment_RNA2'),
            nyami = c('FJ554525.1_Midway_virus_isolate_RML47153',
                      'FJ554526.1_Nyamanini_virus_isolate_tick_39',
                      'HM849038.1_Soybean_cyst_nematode_midway_virus_N-protein',
                      'KF530058.1_Sierra_Nevada_virus',
                      'KM817644.1_Wenzhou_Crab_Virus_1_strain_RBX2_nucleocapsid_N',
                      'KX257488.1_Orinoco_virus_strain_UW1',
                      'KX884406.1_Beihai_rhabdo-like_virus_4_strain_BHJP58499_putative_nucleoprotein',
                      'KX884407.1_Beihai_rhabdo-like_virus_5_strain_BHJP63888_putative_nucleoprotein',
                      'KX884408.1_Beihai_rhabdo-like_virus_3_strain_BHNXC39077_putative_nucleoprotein',
                      'MG550265.1_Soybean_cyst_nematode_nyami-like_virus_isolate_6269845N_nucleoprotein',
                      'MG550266.1_Soybean_cyst_nematode_nyami-like_virus_isolate_6269844N_nucleoprotein',
                      'MG550267.1_Soybean_cyst_nematode_nyami-like_virus_isolate_6232814-16N_nucleoprotein',
                      'MH477287.1_Formica_fusca_virus_1',
                      'MK971153.1_San_Jacinto_virus_isolate_DO-200',
                      'NC_012702.1_Midway_virus',
                      'NC_012703.1_Nyamanini_virus',
                      'NC_024376.1_Sierra_Nevada_virus',
                      'NC_024702.1_Soybean_cyst_nematode_midway_virus_N-protein',
                      'NC_031275.1_Wenzhou_Crab_Virus_1_strain_RBX2_nucleocapsid_N',
                      'NC_032543.1_Beihai_rhabdo-like_virus_4_strain_BHJP58499_putative_nucleoprotein',
                      'NC_032544.1_Beihai_rhabdo-like_virus_5_strain_BHJP63888_putative_nucleoprotein',
                      'NC_043485.1_Orinoco_virus_strain_UW1',
                      'NC_043486.1_Beihai_rhabdo-like_virus_3_strain_BHNXC39077_putative_nucleoprotein'),
            lispi = c('KM817632.1_Lishi_Spider_Virus_2_strain_LSZZ-4_ORF1_ORF1',
                      'NC_031259.1_Lishi_Spider_Virus_2_strain_LSZZ-4_ORF1_ORF1'),
            xinmo = c('KM817638.1_Shuangao_Fly_Virus_2_strain_QSA05_RNA-dependent_RNA_polymerase_L_gene',
                      'NC_043484.1_Drosophila_unispina_virus_1_ORF1'),
            mymona = c('KT598228.1_Soybean_leaf-associated_negative-stranded_RNA_virus_3_isolate_SaNSRV1-1',
                       'LC466007.1_Lentinula_edodes_negative-strand_RNA_virus_1_HG3_RNA',
                       'MH648611.1_Botrytis_cinerea_mymonavirus_1_isolate_Ecan17-2'),
            narna = c('MN034260.1_Narnavirus_sp._isolate_H2_Bulk_Litter_12_scaffold_23_hypothetical_protein_H2BulkLitter1223_000001_gene',
                      'MN035123.1_Narnavirus_sp._isolate_H3_Rhizo_Litter_15_scaffold_62_sequence_1',
                      'MN035214.1_Narnavirus_sp._isolate_H2_Bulk_Litter_12_scaffold_176_sequence_1'),
            polycipi = c('MH477288.1_Lasius_neglectus_virus_2'),
            arto = c('KX884420.1_Hubei_rhabdo-like_virus_8_strain_QTM19395_RNA-dependent_RNA_polymerase_gene',
                     'KX884421.1_Hubei_rhabdo-like_virus_6_strain_QTM24798_hypothetical_protein_and_RNA-dependent_RNA_polymerase_genes',
                     'KX884446.1_Hubei_rhabdo-like_virus_5_strain_WHSFII19440_RNA-dependent_RNA_polymerase_gene',
                     'NC_038269.1_Pteromalus_puparum_negative-strand_RNA_virus_1'),
            unclas = c('KM817631.1_Jingshan_Fly_Virus_2_strain_JSY-4_RNA-dependent_RNA_polymerase_L_gene',
                       'KM817633.1_Sanxia_Water_Strider_Virus_4_strain_SXSSP12_ORF1_ORF1',
                       'KM817637.1_Shuangao_Bedbug_Virus_2_strain_SACC-1_ORF1_ORF1',
                       'KM817639.1_Shuangao_Insect_Virus_6_strain_QSA10_RNA-dependent_RNA_polymerase_L_gene',
                       'KM817641.1_Tacheng_Tick_Virus_6_strain_TCRP-2_ORF1_ORF1',
                       'KM817655.1_Wuhan_Louse_Fly_Virus_8_strain_BFJSC-6_RNA-dependent_RNA_polymerase_L_gene',
                       'KX580901.1_Berant_virus_ORFX_gene',
                       'KX884405.1_Beihai_rhabdo-like_virus_6_strain_BHJJX49420_putative_nucleoprotein',
                       'KX884410.1_Beihai_barnacle_virus_8_strain_BHTH14652_hypothetical_protein_1',
                       'KX884411.1_Beihai_barnacle_virus_7_strain_BHTH16013_hypothetical_protein_1',
                       'KX884412.1_Beihai_rhabdo-like_virus_1_strain_BHTSS15727_hypothetical_protein_1',
                       'KX884413.1_Beihai_rhabdo-like_virus_2_strain_BHTSS7258_hypothetical_protein_1',
                       'KX884414.1_Shayang_ascaridia_galli_virus_2_strain_HC21241_hypothetical_protein_1',
                       'KX884417.1_Hubei_rhabdo-like_virus_3_strain_QCM135517_hypothetical_protein_1',
                       'KX884418.1_Hubei_odonate_virus_10_strain_QTM133243_hypothetical_protein_1',
                       'KX884422.1_Hubei_rhabdo-like_virus_1_strain_QTM25107_hypothetical_protein_1',
                       'KX884440.1_Shayang_ascaridia_galli_virus_2_strain_SYJFC48550_putative_nucleoprotein',
                       'KX884441.1_Wuchan_romanomermis_nematode_virus_2_strain_WCLSXC55347_hypothetical_protein_1',
                       'KX884443.1_Hubei_myriapoda_virus_7_strain_WGML146453_hypothetical_protein_1',
                       'KX884447.1_Hubei_rhabdo-like_virus_2_strain_WHSWHC60518_hypothetical_protein_1',
                       'KX884448.1_Hubei_rhabdo-like_virus_9_strain_WHZHC73015_hypothetical_protein_1',
                       'KX884450.1_Wenling_crustacean_virus_10_strain_WLJQ101844_hypothetical_protein_1',
                       'KX884456.1_Wenling_crustacean_virus_11_strain_WLJQ201798_hypothetical_protein_1',
                       'KX884457.1_Wenling_crustacean_virus_12_strain_WLJQ47777_putative_nucleoprotein',
                       'LN713933.1_Black_grass_varicosavirus-like_virus_segment_RNA_1_1',
                       'MF287670.1_Formica_exsecta_virus_4_ORF5_gp5',
                       'MF444285.1_Sclerotinia_sclerotiorum_negative-stranded_RNA_virus_7_isolate_SsNSRV7_RNA-dependent_RNA_polymerase_gene',
                       'MG489864.1_Lepeophtheirus_virus_LS24',
                       'MK026565.1_Fairlight_virus_segment_1_RNA-dependent_RNA_polymerase_and_glycoprotein_genes',
                       'MK026568.1_Quarantine_head_virus_segment_1_glycoprotein_and_RNA-dependent_RNA_polymerase_genes',
                       'MK584853.1_Penicillium_glabrum_negative-stranded_RNA_virus_1_isolate_CREA-VE-21SY1',
                       'MK584858.1_Penicillium_adametzioides_negative-stranded_RNA_virus_1_isolate_CREA-VE-8SY1',
                       'NC_031088.1_Sanxia_Water_Strider_Virus_4_strain_SXSSP12_ORF1_ORF1',
                       'NC_031217.1_Shuangao_Bedbug_Virus_2_strain_SACC-1_ORF1_ORF1',
                       'NC_031270.1_Tacheng_Tick_Virus_6_strain_TCRP-2_ORF1_ORF1',
                       'NC_031304.1_Wuhan_Tick_Virus_1_strain_X78-2_nucleocapsid_N',
                       'NC_032430.1_Beihai_barnacle_virus_8_strain_BHTH14652_hypothetical_protein_1',
                       'NC_032432.1_Beihai_barnacle_virus_7_strain_BHTH16013_hypothetical_protein_1',
                       'NC_032541.1_Beihai_rhabdo-like_virus_6_strain_BHJJX49420_putative_nucleoprotein',
                       'NC_032555.1_Beihai_rhabdo-like_virus_1_strain_BHTSS15727_hypothetical_protein_1',
                       'NC_032558.1_Beihai_rhabdo-like_virus_2_strain_BHTSS7258_hypothetical_protein_1',
                       'NC_032739.1_Wenling_crustacean_virus_10_strain_WLJQ101844_hypothetical_protein_1',
                       'NC_032781.1_Wenling_crustacean_virus_11_strain_WLJQ201798_hypothetical_protein_1',
                       'NC_032793.1_Wenling_crustacean_virus_12_strain_WLJQ47777_putative_nucleoprotein',
                       'NC_032849.1_Hubei_orthoptera_virus_5_strain_ZCM94262_hypothetical_protein_1',
                       'NC_032929.1_Hubei_rhabdo-like_virus_3_strain_QCM135517_hypothetical_protein_1',
                       'NC_032944.1_Hubei_odonate_virus_10_strain_QTM133243_hypothetical_protein_1',
                       'NC_032971.1_Hubei_rhabdo-like_virus_1_strain_QTM25107_hypothetical_protein_1',
                       'NC_033070.1_Hubei_dimarhabdovirus_virus_1_strain_SCM51525_putative_nucleoprotein',
                       'NC_033182.1_Shayang_ascaridia_galli_virus_2_strain_SYJFC48550_putative_nucleoprotein',
                       'NC_033261.1_Hubei_rhabdo-like_virus_2_strain_WHSWHC60518_hypothetical_protein_1',
                       'NC_033267.1_Hubei_rhabdo-like_virus_9_strain_WHZHC73015_hypothetical_protein_1',
                       'NC_033436.1_Wuchan_romanomermis_nematode_virus_2_strain_WCLSXC55347_hypothetical_protein_1',
                       'NC_033703.1_Xinzhou_dimarhabdovirus_virus_1_strain_XZSJSC65538_hypothetical_protein',
                       'NC_035133.1_Culex_mononega-like_virus_2_strain_mos191gb23464'))
cols <- c("NODE" = "#ff7f0e",
          "SRA" = "#1f77b4",
          "rhabdo"="#2ca02c",
          "nyami"="#2ca02c",
          "lispi"="#2ca02c",
          "xinmo"="#2ca02c",
          "mymona"="#2ca02c",
          "narna"="#2ca02c",
          "polycipi" = "#2ca02c",
          "arto" = "#2ca02c",
          "unclas" = "#7f7f7f",
          "0" = "#7f7f7f")

cols2 <- c("NODE" = "#ff7f0e",
          "SRA" = "#1f77b4",
          "rhabdo"="blue",
          "nyami"="red",
          "lispi"="orange",
          "xinmo"="purple",
          "mymona"="green",
          "narna"="black",
          "polycipi" = "orange",
          "arto" = "lightgreen",
          "unclas" = "#7f7f7f",
          "0" = "#7f7f7f")
rhabdo <- groupOTU(rhabdo, cls)

rhabdotree <- ggtree(rhabdo, aes(color=group), size=5) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = posterior > 0.8), color='black',size=1) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,3) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
rhabdotree
#ggsave("Treeplots/rhabdo.svg",rhabdotree, device = "svg", dpi=380,width = 6.04, height =50,limitsize = FALSE)

r1 <- collapse(rhabdotree, node=849) + geom_point2(aes(subset=(node==849)), shape=17, size=5, fill="brown")
r2 <- collapse(r1, node=1212) + geom_point2(aes(subset=(node==1212)), shape=17, size=5, fill="brown")
r3 <- collapse(r2, node=713) + geom_point2(aes(subset=(node==713)), shape=17, size=5, fill="brown")
r4 <- collapse(r3, node=696) + geom_point2(aes(subset=(node==696)), shape=17, size=5, fill="brown")
r5 <- collapse(r4, node=800) + geom_point2(aes(subset=(node==800)), shape=17, size=5, fill="brown")
rhabdotree <- r5
rhabdotree
#ggsave("Treeplots/rhabdo.svg",rhabdotree, device = "svg", dpi=380,width = 6.04, height =50,limitsize = FALSE)

#rhabdofam <- ggtree(rhabdo, aes(color=group)) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) + 
  #geom_nodepoint(aes(subset = posterior > 0.8), color='black') + 
  #geom_text2(aes(label=node)) + 
  #geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  #xlim(0,5) + 
  #geom_treescale() + 
  #scale_colour_manual(values = cols2)
#ggsave("Treeplots/rhabdofam.svg",rhabdo3, device = "svg", dpi=380,width = 10, height=80,limitsize = FALSE)



########################################### Sinai_NStree:
sinaiNS <- read.beast("../data/MCC/SinaiNS.mcctree")
cls <- list(NODE = c('BeeP-49-2013_NODE_27_length_6327_cov_6_608480_3'),
            Alphatetraviridae = c('EU345431.1_1','KX423453.1_1','NC_001981.1_1','NC_005898.1_1','U18246.1_1','AY594352.1_1'),
            Astroviridae = c('KX907135.1_1','NC_032426.1_1','NC_040647.1_1'),
            Hepeviridae = c('MF190001.1_1','NC_040710.1_1'),
            Bromoviridae = c('NC_009537.1_1'))
sinaiNS <- groupOTU(sinaiNS, cls)
cols <- c("NODE" = "#ff7f0e","Alphatetraviridae"="#2ca02c","Astroviridae"="#2ca02c","Hepeviridae"="#2ca02c","Bromoviridae"="#2ca02c")

sinaiNStree <- ggtree(sinaiNS, aes(color=group), size=5) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = posterior > 0.8), color='black',size=1) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,3) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
sinaiNStree
########################################### Sinai_Stree:
sinaiS <- read.beast("MCCtrees/SinaiS.mcctree")
cls <- list(NODE = c('BP11_NODE_2_length_5933_cov_8664_753785_2', 'BP35_NODE_8_length_6159_cov_3073_471514_3', 'BeeP-34-2013_NODE_11_length_5700_cov_81_138894_3', 'BeeP-34-2013_NODE_8_length_5876_cov_117_977410_3'),
            SRA = c('SRR1239309.Contig_10906_56.5422_length_5835_3', 'SRR1239309.Contig_3402_33.8366_length_5852_3', 'SRR1239310.Contig_18500_86.3069_length_5976_3', 'SRR3927501.Contig_20_4621.33_length_5753_3', 'SRR5109829.Contig_5052_19.3642_length_5832_3', 'SRR5117449.Contig_4792_738.395_length_5543_3', 'SRR6031640.Contig_16980_3187.18_length_5924_3', 'SRR6031648.Contig_11953_83.802_length_5935_3', 'SRR6833958.Contig_36829_28.5413_length_5784_3', 'SRR806508.Contig_22923_17.5662_length_5787_3', 'SRR806508.Contig_4561_28.5866_length_5893_3', 'SRR806550.Contig_15712_2565.35_length_5987_3'),
            Sinai1 = c('HQ871931.2_2','KY465697.1_2','KY465698.1_2','KY465699.1_2','KY465700.1_2','KY465701.1_2','KY465702.1_2','KY465703.1_2','KY465704.1_2','KY465705.1_2','LR596015.1_2','NC_035466.1_2'),
            Sinai2 = c('HQ888865.2_2','KY465706.1_2','KY465707.1_2','KY465708.1_2','KY465709.1_2','KY465710.1_2','KY465711.1_2','KY465712.1_2','KY465713.1_2','LR655824.1_2','NC_035467.1_2'),
            Sinai = c('KM886902.1_2','KM886903.1_2','KM886904.1_2','KM886905.1_2','KX883223.1_2','KY465714.1_2','KY465715.1_2','KY465716.1_1','KY465717.1_2','KY465719.1_1','KY465720.1_2','MG918125.1_2','NC_032433.1_2'),
            SinaiTO = c('KY354241.1_2'),
            SinaiNE = c('KY354242.1_2','NC_035113.1_2'),
            SinaiSA1 = c('KY354243.1_2','NC_035111.1_2'),
            SinaiSA2 = c('KY354244.1_2', 'NC_035112.1_2'),
            Sinai3 = c('MH267699.1_2','MH267700.1_2'))

sinaiS <- groupOTU(sinaiS, cls)
cols <- c("NODE" = "#ff7f0e","SRA" = "#1f77b4","Sinai1" = "#2ca02c","Sinai2" = "#2ca02c","Sinai" = "#2ca02c","SinaiTO" = "#2ca02c","SinaiNE" = "#2ca02c","SinaiSA1" = "#2ca02c","SinaiSA2" = "#2ca02c", "Sinai3" = "#2ca02c")

sinaiStree <- ggtree(sinaiS, aes(color=group),size=5) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = posterior > 0.8), color='black',size=1) + 
  #geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,1) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
sinaiStree
########################################### Toti:
toti <- read.beast("MCCtrees/Toti.mcctree")
cls <- list(NODE = c('BP23_NODE_31_length_5100_cov_72_148022_1_20_2503_-1_ID=1_1_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.518',
                     'BeeP-34-2013_NODE_9_length_5861_cov_32_563451_2_4007_5821_1_ID=2_2_partial=00_start_type=ATG_rbs_motif=CCGCC_rbs_spacer=6bp_gc_cont=0.530'),
            SRA = c('SRR1503019.Contig_22032_542.914_length_5984_1_198_2591_-1_ID=3_1_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.427',
                    'SRR6456322.Contig_9844_124.905_length_7917_1_234_5294_1_ID=6_1_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.588'),
            victorivirus = c('NC_001963.1_2_2574_5090_1_ID=9_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.588',
                             'NC_001964.1_2_2658_5135_1_ID=10_2_partial=00_start_type=ATG_rbs_motif=CCGCC_rbs_spacer=12bp_gc_cont=0.597',
                             'NC_003607.2_2_2605_5112_1_ID=25_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.554',
                             'NC_003876.1_2_2603_5080_1_ID=12_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.543',
                             'NC_005074.1_2_2563_5100_1_ID=16_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.514',
                             'NC_006367.1_3_2818_5316_1_ID=22_3_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.562',
                             'NC_007523.1_2_2386_4875_1_ID=23_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.564',
                             'NC_010246.1_3_2638_5130_1_ID=30_3_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=6bp_gc_cont=0.600',
                             'NC_014823.1_2_2604_5126_1_ID=32_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.595',
                             'NC_020997.1_2_2604_5075_1_ID=36_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.510',
                             'NC_021565.1_2_2664_5192_1_ID=39_2_partial=00_start_type=ATG_rbs_motif=CCGCC_rbs_spacer=8bp_gc_cont=0.576',
                             'NC_023547.1_3_2532_5021_1_ID=42_3_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.583',
                             'NC_024151.1_2_2735_5266_1_ID=44_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.567',
                             'NC_025366.1_2_2661_5144_1_ID=46_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.535',
                             'NC_026140.1_2_2713_5928_1_ID=49_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.608',
                             'NC_027209.1_3_2631_5108_1_ID=58_3_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=6bp_gc_cont=0.585',
                             'NC_028477.1_2_2584_5166_1_ID=64_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.613',
                             'NC_030224.1_2_2554_5037_1_ID=73_2_partial=00_start_type=ATG_rbs_motif=CCGCC_rbs_spacer=15bp_gc_cont=0.592',
                             'NC_030392.1_2_2580_5081_1_ID=75_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=3bp_gc_cont=0.480',
                             'NC_030867.1_2_2572_5055_1_ID=77_2_partial=00_start_type=ATG_rbs_motif=CCGCC_rbs_spacer=9bp_gc_cont=0.562',
                             'NC_038928.1_2_2599_5118_1_ID=96_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.565',
                             'NC_038929.1_2_2672_5176_1_ID=97_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=5bp_gc_cont=0.543',
                             'NC_038930.1_2_2568_5051_1_ID=98_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=4bp_gc_cont=0.589',
                             'NC_040653.1_3_2657_5185_1_ID=115_3_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=6bp_gc_cont=0.618',
                             'NC_040793.1_2_2654_5143_1_ID=116_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=6bp_gc_cont=0.547'),
            toti = c('NC_009224.1_2_2631_5147_1_ID=27_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.518'),
            unclas = c('NC_025218.2_2_261_5321_1_ID=52_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.588',
                       'NC_028948.1_2_2600_5077_1_ID=66_2_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.541',
                       'NC_029989.1_2_2652_5129_1_ID=72_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=6bp_gc_cont=0.559',
                       'NC_032424.1_1_1_4680_1_ID=86_1_partial=10_start_type=Edge_rbs_motif=None_rbs_spacer=None_gc_cont=0.523',
                       'NC_032806.1_2_267_5117_1_ID=88_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=5bp_gc_cont=0.543',
                       'NC_032819.1_1_2_4102_1_ID=89_1_partial=10_start_type=Edge_rbs_motif=None_rbs_spacer=None_gc_cont=0.474',
                       'NC_032851.1_1_101_4273_1_ID=87_1_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.421',
                       'NC_032931.1_1_2_4714_1_ID=90_1_partial=10_start_type=Edge_rbs_motif=None_rbs_spacer=None_gc_cont=0.538',
                       'NC_032948.1_1_72_5027_1_ID=91_1_partial=00_start_type=ATG_rbs_motif=None_rbs_spacer=None_gc_cont=0.578',
                       'NC_040530.1_2_2788_6090_1_ID=114_2_partial=00_start_type=ATG_rbs_motif=CCCC_rbs_spacer=9bp_gc_cont=0.610'))


toti <- groupOTU(toti, cls)
cols <- c("NODE" = "#ff7f0e","SRA" = "#1f77b4","victorivirus" = "#2ca02c","toti" = "#2ca02c","unclas" = "#7f7f7f")


totitree <- ggtree(toti, aes(color=group),size=5) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = posterior > 0.8), color='black',size=1) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,5) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
totitree
########################################### Tymo:
tymo <- read.beast("MCCtrees/Tymo.mcctree")
cls <- list(NODE = c('BeeP-37-2013_NODE_30_length_5285_cov_819_119432_2_733_5283_-1_ID=1_2partial=01start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.516',
                     'BeeP-43-2013_NODE_1_length_6265_cov_122_817873_1_2_5476_1_ID=2_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.554',
                     'BeeP-49-2013_NODE_19_length_8271_cov_90_694044_4_2585_8269_-1_ID=3_4partial=01start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.445',
                     'BeeP-49-2013_NODE_37_length_5188_cov_28_459010_1_3_2561_1_ID=4_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.440'),
            SRA = c('SRR1015503.Contig_4257_2592.36_length_5110_3_1670_5110_-1_ID=5_3partial=01start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.465',
                    'SRR1015531.Contig_28871_50.2275_length_5084_1_3_5084_-1_ID=6_1partial=11start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.461',
                    'SRR1503108.Contig_2059_227.433_length_6823_1_52_6726_-1_ID=7_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.557',
                    'SRR5117444.Contig_2210_2.85334_length_5011_1_2_2641_1_ID=8_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.382'),
            tymoviridae = c('NC_001480.1_1_109_5628_1_ID=16_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.540',
                            'NC_001746.1_1_86_5710_1_ID=17_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.544',
                            'NC_001793.1_1_115_6315_1_ID=18_1partial=00start_type=ATGrbs_motif=GTGrbs_spacer=5bpgc_cont=0.631',
                            'NC_002588.1_1_147_5654_1_ID=21_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.567',
                            'NC_002786.1_1_97_6180_1_ID=22_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.624',
                            'NC_003634.1_1_157_5955_1_ID=26_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.537',
                            'NC_006950.1_1_109_6675_1_ID=36_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.576',
                            'NC_007609.1_1_419_5365_1_ID=39_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.494',
                            'NC_009532.1_1_92_5542_1_ID=50_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.598',
                            'NC_015522.1_1_1_6294_1_ID=82_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.617',
                            'NC_020470.1_1_3_5642_1_ID=107_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.498',
                            'NC_020471.1_1_111_5543_1_ID=108_1partial=00start_type=ATGrbs_motif=GTGrbs_spacer=9bpgc_cont=0.502',
                            'NC_021851.1_1_152_5578_1_ID=111_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.501',
                            'NC_027631.1_1_2_5452_1_ID=130_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.595',
                            'NC_029063.1_1_3_6242_1_ID=157_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.598',
                            'NC_031692.1_1_1_6606_1_ID=172_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.650',
                            'NC_034205.1_1_273_6629_1_ID=184_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.586',
                            'NC_038328.1_1_195_6302_1_ID=314_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.585',
                            'NC_040565.1_1_214_6603_1_ID=349_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.488'),
            Beta = c('NC_001948.1_1_62_6547_1_ID=19_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.430',
                     'NC_003462.2_1_60_6611_1_ID=103_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.427',
                     'NC_003877.1_1_74_5962_1_ID=27_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.393',
                     'NC_005343.1_1_61_6039_1_ID=29_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.455',
                     'NC_008020.1_1_75_6089_1_ID=41_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.448',
                     'NC_009087.2_1_1_6249_1_ID=55_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.455',
                     'NC_009383.1_1_2_6034_1_ID=48_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.435',
                     'NC_009892.1_1_353_6442_1_ID=58_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.432',
                     'NC_009991.1_1_1_6219_1_ID=59_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.451',
                     'NC_011525.1_1_75_5417_1_ID=65_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.452',
                     'NC_012038.1_1_59_5929_1_ID=68_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.446',
                     'NC_013527.1_1_65_6058_1_ID=72_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.443',
                     'NC_014821.1_1_61_6615_1_ID=79_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.415',
                     'NC_017859.1_1_82_6015_1_ID=95_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.473',
                     'NC_018175.1_1_71_5920_1_ID=97_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.497',
                     'NC_018714.1_1_61_6612_1_ID=104_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.422',
                     'NC_023892.1_1_64_6036_1_ID=116_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.474',
                     'NC_025388.1_1_62_6130_1_ID=121_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.445',
                     'NC_026616.2_1_65_5791_1_ID=125_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.425',
                     'NC_027527.1_1_67_5793_1_ID=129_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.431',
                     'NC_028868.1_1_2_6079_1_ID=148_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.435',
                     'NC_028975.1_1_2_6130_1_ID=149_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.437',
                     'NC_029085.1_1_61_6048_1_ID=158_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.453',
                     'NC_029086.1_1_61_6033_1_ID=159_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.454',
                     'NC_029088.1_1_2_6010_1_ID=160_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.487',
                     'NC_031089.1_1_86_5992_1_ID=169_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.419',
                     'NC_035203.1_1_2_6451_1_ID=219_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.426',
                     'NC_038325.1_1_93_6026_1_ID=313_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.401',
                     'NC_038966.1_1_62_6172_1_ID=315_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.453',
                     'NC_040643.1_1_3_6476_1_ID=350_1partial=10start_type=Edgerbs_motif=Nonerbs_spacer=Nonegc_cont=0.426'),
            Alpha = c('NC_003400.1_1_82_4977_1_ID=25_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.530',
                      'NC_028649.1_1_81_4856_1_ID=137_1partial=00start_type=ATGrbs_motif=Nonerbs_spacer=Nonegc_cont=0.512'))
tymo <- groupOTU(tymo, cls)
cols <- c("NODE" = "#ff7f0e","SRA" = "#1f77b4","tymoviridae" = "#2ca02c","Beta" = "#2ca02c","Alpha" = "#2ca02c")

tymotree <- ggtree(tymo, aes(color=group),size=5) +
  #geom_nodelab(aes(label=round(as.numeric(posterior), 2)),alpha=0.5,color='black', size=2, nudge_x=0.01) +
  geom_nodepoint(aes(subset = posterior > 0.8), color='black',size=1) + 
  geom_tiplab() + 
  #theme(legend.position = c(0.1,0.9), 
  #      legend.title = element_blank(), # no title
  #      legend.key = element_blank()) +
  xlim(0,5) + 
  geom_treescale() + 
  scale_colour_manual(values = cols)
tymotree






lay <- rbind(c(1,1,2,2,3,4),
             c(1,1,2,2,5,6),
             c(1,1,2,2,7,8),
             c(1,1,2,2,9,10))

testgrid <- grid.arrange(rhabdotree,picornatree,bunyatree,totitree, orthobeetree,orthosratree,sinaiNStree,sinaiStree,partititree,tymotree,layout_matrix =lay)
#ggsave("Treeplots/testgrid_notiplab.svg",testgrid, device = "svg", dpi=300, height=10,width=15, limitsize=FALSE)