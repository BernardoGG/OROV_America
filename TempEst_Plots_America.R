##OROV TempEst Plots##
##B Gutierrez##
##January 2019##

#Set working directory and load ggplot2 library
setwd("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Temporal_analyses")
getwd()
library(ggplot2)

#Import TemEst files for specific segment and tree topologies
TTR_L <- read.table ("RAxML_OROV_L_TempEst", header = TRUE)
TTR_Mnorecs <- read.table ("RAxML_OROV_Mnorecs_TempEst", header = TRUE)
TTR_Mgen1.nd <- read.table ("RAxML_OROV_M_Genotype1.nd_TempEst", header = TRUE)
TTR_Mgen2 <- read.table ("RAxML_OROV_M_Genotype2_TempEst", header = TRUE)
TTR_S <- read.table ("RAxML_OROV_S_TempEst", header = TRUE)

#TTR_Sg1 <- read.table ("PhyML_OROV_Sg1_1_TempEst", header = TRUE)
#TTR_Sg2 <- read.table ("PhyML_OROV_Sg2_1_TempEst", header = TRUE)

#Plot individual tip-to-root graphs
#L segment
TTR_L_plot <- data.frame(
  Year = TTR_L$date,
  Divergence = TTR_L$distance
)

TempEstL <- ggplot (TTR_L_plot, aes (x = Year, y = Divergence)) + 
  geom_point(color = 'palegreen3') + ylim(0, 0.13) + xlim(1915, 2016.5) +
  geom_smooth (method='lm', formula= y~x, fullrange = TRUE, se = FALSE, color = 'palegreen4') +
  labs (title = "L segment", x = "Year", y = "Root-to-tip divergence")

TempEstL

#M segment
#M segment, no recombinants data set
TTR_M_plot <- data.frame(
  Year = TTR_Mnorecs$date,
  Divergence = TTR_Mnorecs$distance,
  DivergenceMid = TTR_Mnorecs$distance
)

TempEstM <- ggplot (TTR_M_plot, aes (x = Year, y = DivergenceMid)) + 
  geom_point(color = 'orange1') + ylim(0, 0.24) + xlim(1490, 2016.5) +
  geom_smooth (method='lm', formula = y~x, fullrange = TRUE, se = FALSE, color = 'chocolate') +
  labs (title = "M segment", x = "Year", y = "Root-to-tip divergence")

TempEstM

#M lineage 1 data set
TTR_Mgen1.nd_plot <- data.frame(
  Year = TTR_Mgen1.nd$date,
  Divergence = TTR_Mgen1.nd$distance,
  DivergenceMid = TTR_Mgen1.nd$distance
)

TempEstMgen1.nd <- ggplot (TTR_Mgen1.nd_plot, aes (x = Year, y = DivergenceMid)) + 
  geom_point(color = 'lightsalmon') + ylim(0, 0.06) + xlim(1925, 2016.5) +
  geom_smooth (method='lm', formula = y~x, fullrange = TRUE, se = FALSE, color = 'lightsalmon3') +
  labs (title = "M segment, Lineage 1", x = "Year", y = "Root-to-tip divergence")

TempEstMgen1.nd

#M genotype 2 data set
TTR_Mgen2_plot <- data.frame(
  Year = TTR_Mgen2$date,
  Divergence = TTR_Mgen2$distance,
  DivergenceMid = TTR_Mgen2$distance
)

TempEstMgen2 <- ggplot (TTR_Mgen2_plot, aes (x = Year, y = DivergenceMid)) + 
  geom_point(color = 'bisque3') + ylim(0, 0.06) + xlim(1915, 2016.5) +
  geom_smooth (method='lm', formula = y~x, fullrange = TRUE, se = FALSE, color = 'bisque4') +
  labs (title = "M segment, Clade 2", x = "Year", y = "Root-to-tip divergence")

TempEstMgen2

#S segment
TTR_S_plot <- data.frame(
  Year = TTR_S$date,
  Divergence = TTR_S$distance
)

TempEstS <- ggplot (TTR_S_plot, aes (x = Year, y = Divergence)) + 
  geom_point(color = 'blueviolet') + ylim(0, 0.15) + xlim(1915, 2016.5) +
  geom_smooth (method='lm', formula = y~x, fullrange = TRUE, se = FALSE, color = 'mediumpurple') +
  labs (title = "S segment", x = "Year", y = "Root-to-tip divergence")

TempEstS
