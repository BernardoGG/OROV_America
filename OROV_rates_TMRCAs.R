### OROV rates and TMRCAs ###
### B Gutierrez ###

setwd("/Users/user/Documents/DPhil_Research/OROV_America")
getwd()

# Packages

library("ggplot2")
library("Hmisc")

### Read data files and create objects with specific parameters

## S segment - relaxed and strict clocks, inferred rates and TMRCAs
OROVS <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_S_Skyline12_nododgies_600m.csv", sep = ";")
OROVSstrict <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_Strict_clocks/OROV_S_Skyline12_strictclock_norate_300m_log.csv", sep = ";")

# Objects
OROVSs.t <- OROVS[order(OROVS$age.root.),]
OROVSs.tscr <- OROVSstrict[order(OROVSstrict$age.root.),]
OROVSs.cr <- OROVS[order(OROVS$ucld.mean),]
OROVSs.scr <- OROVSstrict[order(OROVSstrict$meanRate),]

## M segment - relaxed and strict clocks, inferred rates and TMRCAs
# Full M data set
# Tip calibrated + lineage node calibrated
OROVML1 <- read.csv("OROV_Temporal_analyses/OROV_full_clock_tests/OROV_M_UCLD_nodecal_L1_nd_600m.csv", sep = ";")
OROVML2 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_BEAST_M_nodecal/OROV_M_nodecal_1bill_log.csv", sep = ";")
OROVMstrictL1 <- read.csv("OROV_Temporal_analyses/OROV_full_clock_tests/OROV_M_Strict_nodecal_L2_nd_600m.csv", sep = ";")
OROVMstrictL2 <- read.csv("OROV_Temporal_analyses/OROV_full_clock_tests/OROV_M_Strict_nodecal_L2_nd_600m.csv", sep = ";")

# Objects
OROVMs_L1.t <- OROVML1[order(OROVML1$age.root.),]
OROVMs_L2.t <- OROVML2[order(OROVML2$age.root.),]

OROVMs_L1.tscr <- OROVMstrictL1[order(OROVMstrictL1$age.root.),]
OROVMs_L2.tscr <- OROVMstrictL2[order(OROVMstrictL2$age.root.),]

OROVMs_ucldL1.cr <- OROVML1[order(OROVML1$ucld.mean),]
OROVMs_ucldL2.cr <- OROVML2[order(OROVML2$ucld.mean),]

OROVMs_strictL1.scr <- OROVMstrictL1[order(OROVMstrictL1$meanRate),]
OROVMs_strictL2.scr <- OROVMstrictL2[order(OROVMstrictL2$meanRate),]

# M - Lineage 1 data set
OROVM1 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_M_Gen1_Skyline12_norates_600m_nd.csv", sep = ";")
OROVM1strict <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_Strict_clocks/OROV_M_Gen1_strictclock_Skyline12_nododgies_300m_log.csv", sep = ";")

# Objects
OROVM1s.t <- OROVM1[order(OROVM1$age.root.),]
OROVM1s.cr <- OROVM1[order(OROVM1$ucld.mean),]
OROVM1s.scr <- OROVM1strict[order(OROVM1strict$meanRate),]

# M - Lineage 2 data set
OROVM2 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_M_Gen2_Skyline12_norates_300m_ALL_log.csv", sep = ";")
OROVM2strict <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_Strict_clocks/OROV_M_Gen2_Skyline12_strictclock_300m_log.csv", sep = ";")

# Objects
OROVM2s.t <- OROVM2[order(OROVM$age.root.),]
OROVM2s.cr <- OROVM2[order(OROVM2$ucld.mean),]
OROVM2s.scr <- OROVM2strict[order(OROVM2strict$meanRate),]

## L segment - relaxed and strict clocks, inferred rates and TMRCAs
OROVL <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_L_Skyline12_1.1bill.csv", sep = ";")
OROVLstrict <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_BEAST_clocks/OROV_Strict_clocks/OROV_L_Skyline12_strict_1.1bill.csv", sep = ";")

# Objects
OROVLs.t <- OROVL[order(OROVL$age.root.),]
OROVLs.tscr <- OROVLstrict[order(OROVLstrict$age.root.),]
OROVLs.cr <- OROVL[order(OROVL$ucld.mean),]
OROVLs.scr <- OROVLstrict[order(OROVLstrict$meanRate),]


### Build data frames for plots

## Burn in 10% of lower values per each feature
# S segment
burnS <- round(length(OROVSs.t$age.root.)*0.05)
burnStemps <- round(length(OROVSs.scr$meanRate)*0.05)

OROVS.tmrca <- OROVSs.t$age.root.[(burnS):(length(OROVSs.t$age.root.)-burnS)]
OROVS.stmrca <- OROVSs.tscr$age.root.[(burnStemps):(length(OROVSs.tscr$age.root.)-burnStemps)]
OROVS.ucld <- OROVSs.cr$ucld.mean[(burnS):(length(OROVSs.cr$ucld.mean)-burnS)]
OROVS.strict <- OROVSs.scr$clock.rate[(burnStemps):(length(OROVSs.scr$clock.rate)-burnStemps)]

# M segment - full data set (tip calibrated, tip + node calibrated)
burnM <- round(length(OROVMs_ucldL1.cr$ucld.mean)*0.05)
burnMtemps <- round(length(OROVMs_strictL1.scr$meanRate)*0.05)

OROVM.tmrcaL1 <- OROVMs_L1.t$age.root.[(burnM):(length(OROVMs_L1.t$age.root.)-burnM)]
OROVM.tmrcaL2 <- OROVMs_L2.t$age.root.[(burnM):(length(OROVMs_L2.t$age.root.)-burnM)]
OROVM.stmrcaL1 <- OROVMs_L1.tscr$age.root.[(burnMtemps):(length(OROVMs_L1.tscr$age.root.)-burnMtemps)]
OROVM.stmrcaL2 <- OROVMs_L2.tscr$age.root.[(burnMtemps):(length(OROVMs_L2.tscr$age.root.)-burnMtemps)]
OROVM.ucldL1 <- OROVMs_ucldL1.cr$meanRate[(burnMtemps):(length(OROVMs_ucldL1.cr$meanRate)-burnMtemps)]
OROVM.ucldL2 <- OROVMs_ucldL2.cr$meanRate[(burnM):(length(OROVMs_ucldL1.cr$meanRate)-burnM)]
OROVM.strictL1 <- OROVMs_strictL1.scr$meanRate[(burnMtemps):(length(OROVMs_strictL1.scr$meanRate)-burnMtemps)]
OROVM.strictL2 <- OROVMs_strictL2.scr$meanRate[(burnMtemps):(length(OROVMs_strictL2.scr$meanRate)-burnMtemps)]

# M segment - Lineage 1
burnM1 <- round(length(OROVM1s.t$age.root.)*0.05)

OROVM1.tmrca <- OROVM1s.t$age.root.[(burnM1):(length(OROVM1s.t$age.root.)-burnM1)]
OROVM1.ucld <- OROVM1s.cr$ucld.mean[(burnM1):(length(OROVM1s.cr$ucld.mean)-burnM1)]
OROVM1.strict <- OROVM1s.scr$clock.rate[(burnM1):(length(OROVM1s.scr$clock.rate)-burnM1)]

# M segment - Lineage 1
burnM2 <- round(length(OROVM2s.t$age.root.)*0.05)

OROVM2.tmrca <- OROVM2s.t$age.root.[(burnM2):(length(OROVM2s.t$age.root.)-burnM2)]
OROVM2.ucld <- OROVM2s.cr$ucld.mean[(burnM2):(length(OROVM2s.cr$ucld.mean)-burnM2)]
OROVM2.strict <- OROVM2s.scr$clock.rate[(burnM2):(length(OROVM2s.scr$clock.rate)-burnM2)]

# L segment
burnL <- round(length(OROVLs.t$age.root.)*0.05)
burnLtemps <- round(length(OROVLs.scr$meanRate)*0.05)

OROVL.tmrca <- OROVLs.t$age.root.[(burnL):(length(OROVLs.t$age.root.)-burnL)]
OROVL.stmrca <- OROVLs.tscr$age.root.[(burnLtemps):(length(OROVLs.tscr$age.root.)-burnLtemps)]
OROVL.ucld <- OROVLs.cr$ucld.mean[(burnL):(length(OROVLs.cr$ucld.mean)-burnL)]
OROVL.strict <- OROVLs.scr$clock.rate[(burnL):(length(OROVLs.scr$clock.rate)-burnL)]

## Create data frames for tree TMRCAs and rates (log 10 values)

cumlengthL <- 3*length(OROVL.tmrca)+length(OROVL.stmrca)
cumlengthM <- 3*length(OROVM.tmrcaL1)+length(OROVM.stmrcaL1)
cumlengthS <- 3*length(OROVS.tmrca)+length(OROVS.stmrca)

# TMRCA data frame
OROV.tmrca <- data.frame(Segment = c(rep("L", cumlengthL),
                                         rep("M", cumlengthM),
                                         rep("S", cumlengthS)),
                             Method = c(rep("TempEst", length(OROVL.tmrca)),
                                      rep("TreeTime", length(OROVL.tmrca)),
                                      rep("Strict Clock", length(OROVL.stmrca)),
                                      rep("UCLD", length(OROVL.tmrca)),
                                      rep("TempEst", length(OROVM.tmrcaL1)),
                                      rep("TreeTime", length(OROVM.tmrcaL1)),
                                      rep("Strict Clock", length(OROVM.stmrcaL1)),
                                      rep("UCLD", length(OROVM.tmrcaL1)),
                                      rep("TempEst", length(OROVS.tmrca)),
                                      rep("TreeTime", length(OROVS.tmrca)),
                                      rep("Strict Clock", length(OROVS.stmrca)),
                                      rep("UCLD", length(OROVS.tmrca))),
                             TMRCAs = c(rep(1855.68, length(OROVL.tmrca)),
                                           rep(1797.5, length(OROVL.tmrca)),
                                           OROVL.stmrca,
                                           OROVL.tmrca,
                                           rep(1551.08, length(OROVM.tmrcaL1)),
                                           rep(1490.8, length(OROVM.tmrcaL1)),
                                           OROVM.stmrcaL1,
                                           OROVM.tmrcaL1,                                       
                                           rep(1932.41, length(OROVS.tmrca)),
                                           rep(1926.1, length(OROVS.tmrca)),
                                           OROVS.stmrca,
                                           OROVS.tmrca)
                        )

# Evolutionary rates data frames
OROV.Mrates <- data.frame(Name = c("TempEst", "TreeTime", rep("Strict.Clock", length(OROVM.strict)),
                                   rep("Strict.Clock.L1", length(OROVM.strictL1)), rep("Strict.Clock.L2", length(OROVM.strictL2)),
                                   rep("UCLD", length(OROVM.ucld)), rep("UCLD.L1", length(OROVM.ucldL1)), rep("UCLD.L2", length(OROVM.ucldL2))),
                          Method = c(log10(4.66e-4), log10(4.13e-4),
                                   log10(OROVM.strict), log10(OROVM.strictL1), log10(OROVM.strictL2),
                                   log10(OROVM.ucld), log10(OROVM.ucldL1), log10(OROVM.ucldL2))
                            )

cumlengthM1 <- 2*length(OROVM1.tmrca)+length(OROVM1.stmrca)
cumlengthM2 <- 2*length(OROVM2.tmrca)+length(OROVM2.stmrca)
cumlengthUCLD <- length(OROVL.ucld)+length(OROVM.ucld)+length(OROVM1.ucld)+length(OROVM2.ucld)+length(OROVS.ucld)
cumlengthstrict <- length(OROVL.strict)+length(OROVM.strict)+length(OROVM1.strict)+length(OROVM2.strict)+length(OROVS.strict)

cumlengthLr <- 3*length(OROVL.ucld)+length(OROVL.strict)
cumlengthMr <- 3*length(OROVM.ucldL1)+length(OROVM.strictL1)
cumlengthSr <- 3*length(OROVS.ucld)+length(OROVS.strict)

OROV.Allrates <- data.frame(Segment = c(rep("L", cumlengthLr),
                                        rep("M", cumlengthMr),
                                        rep("S", cumlengthSr)),
                            Method = c(rep("TempEst", length(OROVL.ucld)),
                                     rep("TreeTime", length(OROVL.ucld)),
                                     rep("Strict Clock", length(OROVL.strict)),
                                     rep("UCLD", length(OROVL.ucld)),
                                     rep("TempEst", length(OROVM.ucldL1)),
                                     rep("TreeTime", length(OROVM.ucldL1)),
                                     rep("Strict Clock", length(OROVM.strictL1)),
                                     rep("UCLD", length(OROVM.ucldL1)),
                                     rep("TempEst", length(OROVS.ucld)),
                                     rep("TreeTime", length(OROVS.ucld)),
                                     rep("Strict Clock", length(OROVS.strict)),
                                     rep("UCLD", length(OROVS.ucld))),
                            Log.Rates = c(rep(log10(1.09e-3), length(OROVL.ucld)),
                                          rep(log10(3.22e-4), length(OROVL.ucld)),
                                          log10(OROVL.strict),
                                          log10(OROVL.ucld),
                                          rep(log10(4.66e-4), length(OROVM.ucldL1)),
                                          rep(log10(4.13e-4), length(OROVM.ucldL1)),
                                          log10(OROVM.strictL1),
                                          log10(OROVM.ucldL1),                                       
                                          rep(log10(8.04e-4), length(OROVS.ucld)),
                                          rep(log10(5.13e-4), length(OROVS.ucld)),
                                          log10(OROVS.strict),
                                          log10(OROVS.ucld))
                            )

OROV.UCLDrates <- data.frame(Segment = c(rep("L", length(OROVL.ucld)),
                                    rep("M", length(OROVM.ucldL2)),
                                    rep("M1", length(OROVM1.ucld)),
                                    rep("M2", length(OROVM2.ucld)),
                                    rep("S", length(OROVS.ucld))),
                        Subst.Rate = c(OROVL.ucld,
                                       OROVM.ucldL2,
                                       OROVM1.ucld,
                                       OROVM2.ucld,
                                       OROVS.ucld)
              )

OROV.UCLDrates2 <- data.frame(Genotype = c(rep("L", length(OROVL.ucld)),
                                        rep("M", length(OROVM.ucldL2)),
                                        rep("M - Lineage 1", length(OROVM1.ucld)),
                                        rep("M - Lineage 2", length(OROVM2.ucld)),
                                        rep("S", length(OROVS.ucld))),
                              Segment = c(rep("L", length(OROVL.ucld)),
                                        rep("M", length(OROVM.ucldL2)+length(OROVM1.ucld)+length(OROVM2.ucld)),
                                        rep("S", length(OROVS.ucld))),
                              Partition = c(rep("Full", length(OROVL.ucld)+length(OROVM.ucldL2)),
                                          rep("Split", length(OROVM1.ucld)+length(OROVM2.ucld)),
                                          rep("Full", length(OROVS.ucld))),
                              Log.Subst.Rate = c(log10(OROVL.ucld),
                                  log10(OROVM.ucldL2),
                                  log10(OROVM1.ucld),
                                  log10(OROVM2.ucld),
                                  log10(OROVS.ucld))
              )

### Generate boxplots

## Comparison between TMRCA estimates between segments and methods
OROV.tmrca$Method <- factor(OROV.tmrca$Method, levels = c("TempEst", "TreeTime", "Strict Clock", "UCLD"))
tmrcaRATESboxplot <- ggplot(OROV.tmrca, aes(x = Method, y = TMRCAs, color = Segment)) + geom_boxplot()
tmrcaRATESboxplot

## Relaxed clock UCLD rates between the three viral segments
ucldRATESboxplot <- ggplot(OROV.UCLDrates, aes(x = Segment, y = log10(Subst.Rate), color = Segment)) + geom_boxplot()
ucldRATESboxplot

## Comparison between evolutionary rate estimates between segments and methods
OROV.Allrates$Method <- factor(OROV.Allrates$Method, levels = c("TempEst", "TreeTime", "Strict Clock", "UCLD"))
allRATESboxplot <- ggplot(OROV.Allrates, aes(x = Method, y = Log.Rates, color = Segment)) + geom_boxplot()
allRATESboxplot

## Comparison of rates based on different methods and node calibrations for the M segment
OROV.Mrates$Name <- factor(OROV.Mrates$Name, levels = c("TempEst", "TreeTime", "Strict.Clock", "Strict.Clock.L1", "Strict.Clock.L2",
                                                        "UCLD", "UCLD.L1", "UCLD.L2"))
MRATESboxplot <- ggplot(OROV.Mrates, aes(x = Name, y = Rate, color = Name)) + geom_boxplot()
MRATESboxplot

#ggsave