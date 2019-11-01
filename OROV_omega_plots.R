##OROV - Polyprotein dN/dS ratio plotting##
##B Gutierrez##

setwd("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection")
getwd()

#Packages
library("ggplot2")
library("coda")

#PART 1 - Site-wise dN/dS ratio comparisons between viruses
#Read data files
OROV_M <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_MEME/OROV_M_MEME/OROV_M_MEME_1.csv", sep = ",")
BUNV_M <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/BUNV_M_MEME_Genbank.csv", sep = ",")
LACV_M <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/LACV_M_MEME_sites.csv", sep = ",")
SBV_M <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/SBV_M_MEME_sites.csv", sep = ",")

#Crete data frame with dN/dS ratios per site per virus
#Beta.plus / alpha
maxlengthvector <- pmax(OROV_M$Site,BUNV_M$Site,LACV_M$Site,SBV_M$Site)

omega_compare <- data.frame(Site = maxlengthvector,
                            OROV = c((OROV_M$beta..1/OROV_M$alpha),rep(0,(length(maxlengthvector)-length(OROV_M$Site)))),
                            BUNV = c((BUNV_M$X.beta..sup....sup..1/BUNV_M$X.alpha.),rep(0,(length(maxlengthvector)-length(BUNV_M$Site)))),
                            LACV = c((LACV_M$X.beta..sup....sup..1/LACV_M$X.alpha.),rep(0,(length(maxlengthvector)-length(LACV_M$Site)))),
                            SBV = c((SBV_M$X.beta..sup....sup..1/SBV_M$X.alpha.),rep(0,(length(maxlengthvector)-length(SBV_M$Site))))
)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

is.infinite.data.frame <- function(x)
  do.call(cbind, lapply(x, is.infinite))

omega_compare <- replace(omega_compare, is.nan(omega_compare), 0)
omega_compare[is.infinite(omega_compare)] <- 0

#Beta.plus
omega_compare2 <- data.frame(Site = maxlengthvector,
                            OROV = c((OROV_M$beta..1),rep(0,(length(maxlengthvector)-length(OROV_M$Site)))),
                            BUNV = c((BUNV_M$X.beta..sup....sup..1),rep(0,(length(maxlengthvector)-length(BUNV_M$Site)))),
                            LACV = c((LACV_M$X.beta..sup....sup..1),rep(0,(length(maxlengthvector)-length(LACV_M$Site)))),
                            SBV = c((SBV_M$X.beta..sup....sup..1),rep(0,(length(maxlengthvector)-length(SBV_M$Site))))
)

#Plot site position versus dN/dS for all viruses
omega_compare_plot <- ggplot() + geom_point(data = omega_compare, aes(x = Site, y = OROV), color = "forestgreen") + 
                                  geom_point(data = omega_compare, aes(x = Site, y = BUNV), color = "darkolivegreen2") +
                                  geom_point(data = omega_compare, aes(x = Site, y = LACV), color = "dodgerblue3") + 
                                  geom_point(data = omega_compare, aes(x = Site, y = SBV), color = "lightblue") +
                                  theme_classic() + xlab('Site') + ylab('dN/dS')
omega_compare_plot

omega_compare_plot2 <- ggplot() + geom_point(data = omega_compare2, aes(x = Site, y = OROV), color = "forestgreen") + 
                                  geom_point(data = omega_compare2, aes(x = Site, y = BUNV), color = "darkolivegreen2") +
                                  geom_point(data = omega_compare2, aes(x = Site, y = LACV), color = "dodgerblue3") + 
                                  geom_point(data = omega_compare2, aes(x = Site, y = SBV), color = "lightblue") +
                                  theme_classic() + xlab('Site') + ylab('Beta')
omega_compare_plot2


#PART 2 - Mean dN/dS per domain
#Import data from BEAST log files
CoreAll <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_coreregion.csv", sep = ";")
Core1 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_coreregion_Gen1nd.csv", sep = ";")
Core2 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_coreregion_Gen2.csv", sep = ";")
HeadAll <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_head.csv", sep = ";")
Head1 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_Gc_head_Gen1nd.csv", sep = ";")
Head2 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_head_Gen2.csv", sep = ";")
Stalk1All <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_stalk1.csv", sep = ";")
Stalk11 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_stalk1_Gen1nd.csv", sep = ";")
Stalk12 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_stalk1_Gen2.csv", sep = ";")
Stalk2All <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_stalk2.csv", sep = ";")
Stalk21 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_stalk2_Gen1nd.csv", sep = ";")
Stalk22 <- read.csv("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Selection/OROV_polyprotein_analyses/OROV_omega_Gc_stalk2_Gen2.csv", sep = ";")

#Normalise dN/dS ratio from conditional counts to expected counts - total_X/utotal_X
CoreAll_dN <- as.vector(CoreAll$OROV_Gc_coreregion.total_N/CoreAll$OROV_Gc_coreregion.utotal_N)
CoreAll_dS <- as.vector(CoreAll$OROV_Gc_coreregion.total_S/CoreAll$OROV_Gc_coreregion.utotal_S)
Core1_dN <- as.vector(Core1$OROV_Gc_coreregion_nododgies.total_N/Core1$OROV_Gc_coreregion_nododgies.utotal_N)
Core1_dS <- as.vector(Core1$OROV_Gc_coreregion_nododgies.total_S/Core1$OROV_Gc_coreregion_nododgies.utotal_S)
Core2_dN <- as.vector(Core2$OROV_Gc_coreregion_Gen2.total_N/Core2$OROV_Gc_coreregion_Gen2.utotal_N)
Core2_dS <- as.vector(Core2$OROV_Gc_coreregion_Gen2.total_S/Core2$OROV_Gc_coreregion_Gen2.utotal_S)
HeadAll_dN <- as.vector(HeadAll$OROV_Gc_head.total_N/HeadAll$OROV_Gc_head.utotal_N)
HeadAll_dS <- as.vector(HeadAll$OROV_Gc_head.total_S/HeadAll$OROV_Gc_head.utotal_S)
Head1_dN <- as.vector(Head1$OROV_Gc_head_nododgies.total_N/Head1$OROV_Gc_head_nododgies.utotal_N)
Head1_dS <- as.vector(Head1$OROV_Gc_head_nododgies.total_S/Head1$OROV_Gc_head_nododgies.utotal_S)
Head2_dN <- as.vector(Head2$OROV_Gc_head_Gen2.total_N/Head2$OROV_Gc_head_Gen2.utotal_N)
Head2_dS <- as.vector(Head2$OROV_Gc_head_Gen2.total_S/Head2$OROV_Gc_head_Gen2.utotal_S)
Stalk1All_dN <- as.vector(Stalk1All$OROV_Gc_stalk1.total_N/Stalk1All$OROV_Gc_stalk1.utotal_N)
Stalk1All_dS <- as.vector(Stalk1All$OROV_Gc_stalk1.total_S/Stalk1All$OROV_Gc_stalk1.utotal_S)
Stalk11_dN <- as.vector(Stalk11$OROV_Gc_stalk1_nododgies.total_N/Stalk11$OROV_Gc_stalk1_nododgies.utotal_N)
Stalk11_dS <- as.vector(Stalk11$OROV_Gc_stalk1_nododgies.total_S/Stalk11$OROV_Gc_stalk1_nododgies.utotal_S)
Stalk12_dN <- as.vector(Stalk12$OROV_Gc_stalk1_Gen2.total_N/Stalk12$OROV_Gc_stalk1_Gen2.utotal_N)
Stalk12_dS <- as.vector(Stalk12$OROV_Gc_stalk1_Gen2.total_S/Stalk12$OROV_Gc_stalk1_Gen2.utotal_S)
Stalk2All_dN <- as.vector(Stalk2All$OROV_Gc_stalk2.total_N/Stalk2All$OROV_Gc_stalk2.utotal_N)
Stalk2All_dS <- as.vector(Stalk2All$OROV_Gc_stalk2.total_S/Stalk2All$OROV_Gc_stalk2.utotal_S)
Stalk21_dN <- as.vector(Stalk21$OROV_Gc_stalk2_nododgies.total_N/Stalk21$OROV_Gc_stalk2_nododgies.utotal_N)
Stalk21_dS <- as.vector(Stalk21$OROV_Gc_stalk2_nododgies.total_S/Stalk21$OROV_Gc_stalk2_nododgies.utotal_S)
Stalk22_dN <- as.vector(Stalk22$OROV_Gc_stalk2_Gen2.total_N/Stalk22$OROV_Gc_stalk2_Gen2.utotal_N)
Stalk22_dS <- as.vector(Stalk22$OROV_Gc_stalk2_Gen2.total_S/Stalk22$OROV_Gc_stalk2_Gen2.utotal_S)

#Create data frame for dN/dS ratios (MCMC sample)
omega_all <- data.frame(Number = c(1:length(Core1_dN)),
                        Core = CoreAll_dN/CoreAll_dS,
                        CoreG1 = Core1_dN/Core1_dS,
                        CoreG2 = Core2_dN/Core2_dS,
                        Head = HeadAll_dN/HeadAll_dS,
                        HeadG1 = c((Head1_dN/Head1_dS),rep(NA,(length(Core1_dN)-length(Head1_dN)))),
                        HeadG2 = Head2_dN/Head2_dS,
                        Stalk1 = Stalk1All_dN/Stalk1All_dS,
                        Stalk1G1 = Stalk11_dN/Stalk11_dS,
                        Stalk1G2 = Stalk12_dN/Stalk12_dS,
                        Stalk2 = Stalk2All_dN/Stalk2All_dS,
                        Stalk2G1 = Stalk21_dN/Stalk21_dS,
                        Stalk2G2 = Stalk22_dN/Stalk22_dS
)

#Estimate 95% HPD for dN/dS ratios
#Commands for Head Gen 1: (mean(omega_all$HeadG1)) , (HPDinterval(as.mcmc(omega_all$HeadG1))[1,2]) , (HPDinterval(as.mcmc(omega_all$HeadG1))[1,1])
omega_plot_values <- data.frame(ID = c("H","H1","H2","S1","S11","S12","S2","S21","S22","C","C1","C2"),
                               Domain = c(rep("Head",3),rep("Stalk 1",3),rep("Stalk 2",3),rep("Core",3)),
                               Genotype = c("All","Genotype 1","Genotype 2","All","Genotype 1","Genotype 2",
                                            "All","Genotype 1","Genotype 2","All","Genotype 1","Genotype 2"),
                               Mean = c((mean(omega_all$Head)),(mean(omega_all$HeadG1, na.rm = TRUE)),(mean(omega_all$HeadG2)),
                                        (mean(omega_all$Stalk1)),(mean(omega_all$Stalk1G1)),(mean(omega_all$Stalk1G2)),
                                        (mean(omega_all$Stalk2)),(mean(omega_all$Stalk2G1)),(mean(omega_all$Stalk2G2)),
                                        (mean(omega_all$Core)),(mean(omega_all$CoreG1)),(mean(omega_all$CoreG2))),
                               Upper = c((HPDinterval(as.mcmc(omega_all$Head))[1,2]),(HPDinterval(as.mcmc(omega_all$HeadG1),na.rm = TRUE)[1,2]),(HPDinterval(as.mcmc(omega_all$HeadG2))[1,2]),
                                         (HPDinterval(as.mcmc(omega_all$Stalk1))[1,2]),(HPDinterval(as.mcmc(omega_all$Stalk1G1))[1,2]),(HPDinterval(as.mcmc(omega_all$Stalk1G2))[1,2]),
                                         (HPDinterval(as.mcmc(omega_all$Stalk2))[1,2]),(HPDinterval(as.mcmc(omega_all$Stalk2G1))[1,2]),(HPDinterval(as.mcmc(omega_all$Stalk2G2))[1,2]),
                                         (HPDinterval(as.mcmc(omega_all$Core))[1,2]),(HPDinterval(as.mcmc(omega_all$CoreG1))[1,2]),(HPDinterval(as.mcmc(omega_all$CoreG2))[1,2])),
                               Lower = c((HPDinterval(as.mcmc(omega_all$Head))[1,1]),(HPDinterval(as.mcmc(omega_all$HeadG1),na.rm = TRUE)[1,1]),(HPDinterval(as.mcmc(omega_all$HeadG2))[1,1]),
                                         (HPDinterval(as.mcmc(omega_all$Stalk1))[1,1]),(HPDinterval(as.mcmc(omega_all$Stalk1G1))[1,1]),(HPDinterval(as.mcmc(omega_all$Stalk1G2))[1,1]),
                                         (HPDinterval(as.mcmc(omega_all$Stalk2))[1,1]),(HPDinterval(as.mcmc(omega_all$Stalk2G1))[1,1]),(HPDinterval(as.mcmc(omega_all$Stalk2G2))[1,1]),
                                         (HPDinterval(as.mcmc(omega_all$Core))[1,1]),(HPDinterval(as.mcmc(omega_all$CoreG1))[1,1]),(HPDinterval(as.mcmc(omega_all$CoreG2))[1,1]))
)

#Generate plots
omega_plot_values$ID <- factor(omega_plot_values$ID, levels = omega_plot_values$ID)
omegasplot <- ggplot(omega_plot_values, aes(x = ID, y = Mean, ymin = Lower, ymax = Upper, color = Genotype)) +
                    geom_pointrange() + theme_classic() + xlab('Domain') + ylab('dN/dS')
omegasplot

#PART 3 - Site-wise ratios (beta - alpha) inferred through FUBAR
fubar_L <- read.csv("./OROV_FUBAR/OROV_L_FUBAR.csv", sep = ",", header = TRUE)
fubar_M <- read.csv("./OROV_FUBAR/OROV_M_norecombs_nododgies_FUBAR.csv", sep = ",", header = TRUE)
fubar_M1 <- read.csv("./OROV_FUBAR/OROV_M_Lin1_nododgies_FUBAR.csv", sep = ",", header = TRUE)
fubar_M2 <- read.csv("./OROV_FUBAR/OROV_M_Lin2_FUBAR.csv", sep = ",", header = TRUE)
fubar_gc <- read.csv("./OROV_FUBAR/OROV_M_Gc_norecombs_FUBAR.csv", sep = ",", header = TRUE)
fubar_S <- read.csv("./OROV_FUBAR/OROV_S_FUBAR.csv", sep = ",", header = TRUE)

fubar_L_dotplot <- ggplot(fubar_L, aes(x = Site, y = X.beta...alpha.)) +
  geom_bar(stat = "identity", color = "gray85") + theme_classic() + geom_point()  + xlab('') + ylab('Beta-Alpha') +
  geom_hline(yintercept=1, linetype="dashed")

fubar_M_dotplot <- ggplot(fubar_M, aes(x = Site, y = X.beta...alpha.)) +
  geom_bar(stat = "identity", color = "gray85") + theme_classic() + geom_point()  + xlab('') + ylab('Beta-Alpha') +
  geom_hline(yintercept=1, linetype="dashed")

fubar_M1_dotplot <- ggplot(fubar_M1, aes(x = Site, y = X.beta...alpha.)) +
  geom_bar(stat = "identity", color = "gray85") + theme_classic() + geom_point()  + xlab('') + ylab('Beta-Alpha') +
  geom_hline(yintercept=1, linetype="dashed")

fubar_M2_dotplot <- ggplot(fubar_M2, aes(x = Site, y = X.beta...alpha.)) +
  geom_bar(stat = "identity", color = "gray85") + theme_classic() + geom_point()  + xlab('') + ylab('Beta-Alpha') +
  geom_hline(yintercept=1, linetype="dashed")

fubar_gc_dotplot <- ggplot(fubar_gc, aes(x = Site, y = X.beta...alpha.)) +
  geom_bar(stat = "identity", color = "gray85") + theme_classic() + geom_point()  + xlab('') + ylab('Beta-Alpha') +
  geom_hline(yintercept=1, linetype="dashed")

fubar_S_dotplot <- ggplot(fubar_S, aes(x = Site, y = X.beta...alpha.)) +
  geom_bar(stat = "identity", color = "gray85") + theme_classic() + geom_point()  + xlab('') + ylab('Beta-Alpha') +
  geom_hline(yintercept=1, linetype="dashed")

fubar_L_dotplot
fubar_M_dotplot
fubar_M1_dotplot
fubar_M2_dotplot
fubar_gc_dotplot
fubar_S_dotplot
