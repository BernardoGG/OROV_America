##OROV rates and TMRCAs##
##B Gutierrez##

setwd("/Users/user/Documents/DPhil_Research/OROV_America")
getwd()

#Packages
library("ggplot2")
library("Hmisc")

#Read data files
OROVS <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_S_Skyline12_norate_300m_ALL_log.csv", sep = ";")
OROVSs.t <- OROVS[order(OROVS$treeModel.rootHeight),]
OROVSs.cr <- OROVS[order(OROVS$ucld.mean),]
OROVSs.mr <- OROVS[order(OROVS$meanRate),]

OROVM <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_BEAST_M_nodecal/OROV_M_nodecal_600m_log.csv", sep = ";")
OROVM <- OROVM[-1,] #Delete starting value row
OROVMs.t <- OROVM[order(OROVM$treeModel.rootHeight),]
OROVMs.cr <- OROVM[order(OROVM$ucld.mean),]
OROVMs.mr <- OROVM[order(OROVM$meanRate),]

OROVM1 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_M_Gen1_Skyline12_norates_300m_ALL_log.csv", sep = ";")
OROVM1s.t <- OROVM1[order(OROVM$treeModel.rootHeight),]
OROVM1s.cr <- OROVM1[order(OROVM$ucld.mean),]
OROVM1s.mr <- OROVM1[order(OROVM$meanRate),]

OROVM2 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_M_Gen2_Skyline12_norates_300m_ALL_log.csv", sep = ";")
OROVM2s.t <- OROVM2[order(OROVM$treeModel.rootHeight),]
OROVM2s.cr <- OROVM2[order(OROVM$ucld.mean),]
OROVM2s.mr <- OROVM2[order(OROVM$meanRate),]

OROVL <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_L_Skyline12_norate_300m_ALL_log.csv", sep = ";")
OROVLs.t <- OROVL[order(OROVL$treeModel.rootHeight),]
OROVLs.cr <- OROVL[order(OROVL$ucld.mean),]
OROVLs.mr <- OROVL[order(OROVL$meanRate),]

#Build objects for TMRCAs and rates
#Burn in 10% of lower values per each feature
burnL <- round(length(OROVLs.t$treeModel.rootHeight)*0.05)
OROVL.tmrca <- OROVLs.t$treeModel.rootHeight[(burnL):(length(OROVLs.t$treeModel.rootHeight)-burnL)]
OROVL.ucld <- OROVLs.cr$ucld.mean[(burnL):(length(OROVLs.cr$ucld.mean)-burnL)]
OROVL.mean <- OROVLs.mr$meanRate[(burnL):(length(OROVLs.mr$meanRate)-burnL)]

burnM <- round(length(OROVMs.t$treeModel.rootHeight)*0.05)
OROVM.tmrca <- OROVMs.t$treeModel.rootHeight[(burnM):(length(OROVMs.t$treeModel.rootHeight)-burnM)]
OROVM.ucld <- OROVMs.cr$ucld.mean[(burnM):(length(OROVMs.cr$ucld.mean)-burnM)]
OROVM.mean <- OROVMs.mr$meanRate[(burnM):(length(OROVMs.mr$meanRate)-burnM)]

burnM1 <- round(length(OROVM1s.t$treeModel.rootHeight)*0.05)
OROVM1.tmrca <- OROVM1s.t$treeModel.rootHeight[(burnM1):(length(OROVM1s.t$treeModel.rootHeight)-burnM1)]
OROVM1.ucld <- OROVM1s.cr$ucld.mean[(burnM1):(length(OROVM1s.cr$ucld.mean)-burnM1)]
OROVM1.mean <- OROVM1s.mr$meanRate[(burnM1):(length(OROVM1s.mr$meanRate)-burnM1)]

burnM2 <- round(length(OROVM2s.t$treeModel.rootHeight)*0.05)
OROVM2.tmrca <- OROVM2s.t$treeModel.rootHeight[(burnM2):(length(OROVM2s.t$treeModel.rootHeight)-burnM2)]
OROVM2.ucld <- OROVM2s.cr$ucld.mean[(burnM2):(length(OROVM2s.cr$ucld.mean)-burnM2)]
OROVM2.mean <- OROVM2s.mr$meanRate[(burnM2):(length(OROVM2s.mr$meanRate)-burnM2)]

burnS <- round(length(OROVSs.t$treeModel.rootHeight)*0.05)
OROVS.tmrca <- OROVSs.t$treeModel.rootHeight[(burnS):(length(OROVSs.t$treeModel.rootHeight)-burnS)]
OROVS.ucld <- OROVSs.cr$ucld.mean[(burnS):(length(OROVSs.cr$ucld.mean)-burnS)]
OROVS.mean <- OROVSs.mr$meanRate[(burnS):(length(OROVSs.mr$meanRate)-burnS)]

#Figures
#Create data frames for tree TMRCAs (log 10 values) and rates
OROVtmrca <- data.frame(Segment = c(rep("L", length(OROVL.tmrca)),
                                    rep("M", length(OROVM.tmrca)),
                                    rep("M1", length(OROVM1.tmrca)),
                                    rep("M2", length(OROVM2.tmrca)),
                                    rep("S", length(OROVS.tmrca))),
                        TMRCA = c((OROVL.tmrca),
                                   (OROVM.tmrca),
                                   (OROVM1.tmrca),
                                   (OROVM2.tmrca),
                                   (OROVS.tmrca))
              )

OROV.UCLDrates <- data.frame(Segment = c(rep("L", length(OROVL.ucld)),
                                    rep("M", length(OROVM.ucld)),
                                    rep("M1", length(OROVM1.ucld)),
                                    rep("M2", length(OROVM2.ucld)),
                                    rep("S", length(OROVS.ucld))),
                        Subst.Rate = c(OROVL.ucld,
                                       OROVM.ucld,
                                       OROVM1.ucld,
                                       OROVM2.ucld,
                                       OROVS.ucld)
              )

OROV.UCLDrates2 <- data.frame(Genotype = c(rep("L", length(OROVL.ucld)),
                                        rep("M", length(OROVM.ucld)),
                                        rep("M - Genotype 1", length(OROVM1.ucld)),
                                        rep("M - Genotype 2", length(OROVM2.ucld)),
                                        rep("S", length(OROVS.ucld))),
                              Segment = c(rep("L", length(OROVL.ucld)),
                                        rep("M", length(OROVM.ucld)+length(OROVM1.ucld)+length(OROVM2.ucld)),
                                        rep("S", length(OROVS.ucld))),
                              Partition = c(rep("Full", length(OROVL.ucld)+length(OROVM.ucld)),
                                          rep("Split", length(OROVM1.ucld)+length(OROVM2.ucld)),
                                          rep("Full", length(OROVS.ucld))),
                              Log.Subst.Rate = c(log10(OROVL.ucld),
                                  log10(OROVM.ucld),
                                  log10(OROVM1.ucld),
                                  log10(OROVM2.ucld),
                                  log10(OROVS.ucld))
              )

#Violin and/or boxplots
TMRCAviolplot <- ggplot(OROVtmrca, aes(x = Segment, y = TMRCA, fill = Segment)) + geom_violin() +
                  scale_fill_brewer(palette = "Blues") + theme_bw() +
                  stat_summary(fun.data = mean_sdl, geom = "pointrange")
TMRCAviolplot

TMRCAboxplot <- ggplot(OROVtmrca, aes(x = Segment, y = TMRCA)) + geom_boxplot()
TMRCAboxplot

RATESviolplot <- ggplot(OROV.UCLDrates, aes(x = Segment, y = Subst.Rate, fill = Segment)) + geom_violin() +
                  theme_bw() + scale_fill_brewer(palette = "Greens") +
                  stat_summary(fun.data = mean_sdl, geom = "pointrange")
RATESviolplot

logRATESviolplot <- ggplot(OROV.UCLDrates2, aes(x = Segment, y = Log.Subst.Rate, fill = Segment)) + geom_violin() +
  theme_bw() + scale_fill_brewer(palette = "Reds")
logRATESviolplot

RATESboxplot <- ggplot(OROV.UCLDrates, aes(x = Segment, y = Subst.Rate, color = Segment)) + geom_boxplot()
RATESboxplot

#Manuscript figure
logRATESboxplot <- ggplot(OROV.UCLDrates2, aes(x = Genotype, y = Log.Subst.Rate, color = Segment, linetype = Partition)) + geom_boxplot() + theme_classic()
#                    + scale_fill_manual(values = alpha(c("red","green","blue"), 0.2))
logRATESboxplot

#ggsave