##OROV SimPlots##
##B Gutierrez##
##January 2019##

#Set working directory and load ggplot2 library
setwd("/Users/user/Documents/DPhil_Research/OROV_America/OROV_Recombination")
getwd()
library('adegenet')
library('ape')
library('spider')
library('ggplot2')

#Import sequences and transform to DNAbin object
L1.seqs <- fasta2DNAbin('OROV_v_IQTV_L_BRA1960.fasta')
L2.seqs <- fasta2DNAbin('OROV_v_IQTV_L_BRA2000.fasta')
L3.seqs <- fasta2DNAbin('OROV_v_IQTV_L_BRA2009.fasta')
L4.seqs <- fasta2DNAbin('OROV_v_IQTV_L_PAN1999.fasta')
L5.seqs <- fasta2DNAbin('OROV_v_IQTV_L_ECU2016.fasta')
L6.seqs <- fasta2DNAbin('OROV_v_IQTV_L_PAN1989.fasta')
L7.seqs <- fasta2DNAbin('OROV_v_IQTV_L_PER2000.fasta')
L8.seqs <- fasta2DNAbin('OROV_v_IQTV_L_PER2008.fasta')

G1a.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G1_BRA1960.fasta')
G1b.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G1_BRA2000.fasta')
G1c.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G1_BRA2009.fasta')
G1d.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G1_PAN1999.fasta')
G2a.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G2_ECU2016.fasta')
G2b.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G2_PAN1989.fasta')
G2c.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G2_PER2000.fasta')
G2d.seqs <- fasta2DNAbin('OROV_v_IQTV_M_G2_PER2008.fasta')

S1.seqs <- fasta2DNAbin('OROV_v_IQTV_S_BRA1960.fasta')
S2.seqs <- fasta2DNAbin('OROV_v_IQTV_S_BRA2000.fasta')
S3.seqs <- fasta2DNAbin('OROV_v_IQTV_S_BRA2009.fasta')
S4.seqs <- fasta2DNAbin('OROV_v_IQTV_S_PAN1999.fasta')
S5.seqs <- fasta2DNAbin('OROV_v_IQTV_S_ECU2016.fasta')
S6.seqs <- fasta2DNAbin('OROV_v_IQTV_S_PAN1989.fasta')
S7.seqs <- fasta2DNAbin('OROV_v_IQTV_S_PER2000.fasta')
S8.seqs <- fasta2DNAbin('OROV_v_IQTV_S_PER2008.fasta')

#Set sliding windows on alignments to estimate genetic distances
L1_wndws <- slidingWindow(L1.seqs, width = round((length(L1.seqs)/2)*0.07), interval = 20)
length(L1_wndws)
L2_wndws <- slidingWindow(L2.seqs, width = round((length(L2.seqs)/2)*0.07), interval = 20)
length(L2_wndws)
L3_wndws <- slidingWindow(L3.seqs, width = round((length(L3.seqs)/2)*0.07), interval = 20)
length(L3_wndws)
L4_wndws <- slidingWindow(L4.seqs, width = round((length(L4.seqs)/2)*0.07), interval = 20)
length(L4_wndws)
L5_wndws <- slidingWindow(L5.seqs, width = round((length(L5.seqs)/2)*0.07), interval = 20)
length(L5_wndws)
L6_wndws <- slidingWindow(L6.seqs, width = round((length(L6.seqs)/2)*0.07), interval = 20)
length(L6_wndws)
L7_wndws <- slidingWindow(L7.seqs, width = round((length(L7.seqs)/2)*0.07), interval = 20)
length(L7_wndws)
L8_wndws <- slidingWindow(L8.seqs, width = round((length(L8.seqs)/2)*0.07), interval = 20)
length(L8_wndws)

G1a_wndws <- slidingWindow(G1a.seqs, width = round((length(G1a.seqs)/2)*0.07), interval = 20)
length(G1a_wndws)
G1b_wndws <- slidingWindow(G1b.seqs, width = round((length(G1b.seqs)/2)*0.07), interval = 20)
length(G1b_wndws)
G1c_wndws <- slidingWindow(G1c.seqs, width = round((length(G1c.seqs)/2)*0.07), interval = 20)
length(G1c_wndws)
G1d_wndws <- slidingWindow(G1d.seqs, width = round((length(G1d.seqs)/2)*0.07), interval = 20)
length(G1d_wndws)
G2a_wndws <- slidingWindow(G2a.seqs, width = round((length(G2a.seqs)/2)*0.07), interval = 20)
length(G2a_wndws)
G2b_wndws <- slidingWindow(G2b.seqs, width = round((length(G2b.seqs)/2)*0.07), interval = 20)
length(G2b_wndws)
G2c_wndws <- slidingWindow(G2c.seqs, width = round((length(G2c.seqs)/2)*0.07), interval = 20)
length(G2c_wndws)
G2d_wndws <- slidingWindow(G2d.seqs, width = round((length(G2d.seqs)/2)*0.07), interval = 20)
length(G2d_wndws)

S1_wndws <- slidingWindow(S1.seqs, width = round((length(S1.seqs)/2)*0.07), interval = 20)
length(S1_wndws)
S2_wndws <- slidingWindow(S2.seqs, width = round((length(S2.seqs)/2)*0.07), interval = 20)
length(S2_wndws)
S3_wndws <- slidingWindow(S3.seqs, width = round((length(S3.seqs)/2)*0.07), interval = 20)
length(S3_wndws)
S4_wndws <- slidingWindow(S4.seqs, width = round((length(S4.seqs)/2)*0.07), interval = 20)
length(S4_wndws)
S5_wndws <- slidingWindow(S5.seqs, width = round((length(S5.seqs)/2)*0.07), interval = 20)
length(S5_wndws)
S6_wndws <- slidingWindow(S6.seqs, width = round((length(S6.seqs)/2)*0.07), interval = 20)
length(S6_wndws)
S7_wndws <- slidingWindow(S7.seqs, width = round((length(S7.seqs)/2)*0.07), interval = 20)
length(S7_wndws)
S8_wndws <- slidingWindow(S8.seqs, width = round((length(S8.seqs)/2)*0.07), interval = 20)
length(S8_wndws)

#Estimate mean genetic distance per window
L1_dist <- sapply(L1_wndws, dist.dna)
x <- length(L1_wndws)
L1_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L1_dist
)
L2_dist <- sapply(L2_wndws, dist.dna)
x <- length(L2_wndws)
L2_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L2_dist
)
L3_dist <- sapply(L3_wndws, dist.dna)
x <- length(L3_wndws)
L3_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L3_dist
)
L4_dist <- sapply(L4_wndws, dist.dna)
x <- length(L4_wndws)
L4_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L4_dist
)
L5_dist <- sapply(L5_wndws, dist.dna)
x <- length(L5_wndws)
L5_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L5_dist
)
L6_dist <- sapply(L6_wndws, dist.dna)
x <- length(L6_wndws)
L6_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L6_dist
)
L7_dist <- sapply(L7_wndws, dist.dna)
x <- length(L7_wndws)
L7_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L7_dist
)
L8_dist <- sapply(L8_wndws, dist.dna)
x <- length(L8_wndws)
L8_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = L8_dist
)

G1a_dist <- sapply(G1a_wndws, dist.dna)
x <- length(G1a_wndws)
G1a_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1a_dist
)
G1b_dist <- sapply(G1b_wndws, dist.dna)
x <- length(G1b_wndws)
G1b_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1b_dist
)
G1c_dist <- sapply(G1c_wndws, dist.dna)
x <- length(G1c_wndws)
G1c_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1c_dist
)
G1d_dist <- sapply(G1d_wndws, dist.dna)
x <- length(G1d_wndws)
G1d_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1d_dist
)
G2a_dist <- sapply(G2a_wndws, dist.dna)
x <- length(G2a_wndws)
G2a_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2a_dist
)
G2b_dist <- sapply(G2b_wndws, dist.dna)
x <- length(G2b_wndws)
G2b_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2b_dist
)
G2c_dist <- sapply(G2c_wndws, dist.dna)
x <- length(G2c_wndws)
G2c_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2c_dist
)
G2d_dist <- sapply(G2d_wndws, dist.dna)
x <- length(G2d_wndws)
G2d_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2d_dist
)

S1_dist <- sapply(S1_wndws, dist.dna)
x <- length(S1_wndws)
S1_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S1_dist
)
S2_dist <- sapply(S2_wndws, dist.dna)
x <- length(S2_wndws)
S2_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S2_dist
)
S3_dist <- sapply(S3_wndws, dist.dna)
x <- length(S3_wndws)
S3_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S3_dist
)
S4_dist <- sapply(S4_wndws, dist.dna)
x <- length(S4_wndws)
S4_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S4_dist
)
S5_dist <- sapply(S5_wndws, dist.dna)
x <- length(S5_wndws)
S5_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S5_dist
)
S6_dist <- sapply(S6_wndws, dist.dna)
x <- length(S6_wndws)
S6_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S6_dist
)
S7_dist <- sapply(S7_wndws, dist.dna)
x <- length(S7_wndws)
S7_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S7_dist
)
S8_dist <- sapply(S8_wndws, dist.dna)
x <- length(S8_wndws)
S8_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = S8_dist
)

#Plot genetic distances for each segment
ggplot() +
  geom_line(data = L1_dist, aes(x = Window.Start, y = Dist), color = "olivedrab") +
  geom_line(data = L2_dist, aes(x = Window.Start, y = Dist), color = "olivedrab4") +
  geom_line(data = L3_dist, aes(x = Window.Start, y = Dist), color = "olivedrab3") +
  geom_line(data = L4_dist, aes(x = Window.Start, y = Dist), color = "olivedrab2") +
  geom_line(data = L5_dist, aes(x = Window.Start, y = Dist), color = "palegreen4") +
  geom_line(data = L6_dist, aes(x = Window.Start, y = Dist), color = "palegreen3") +
  geom_line(data = L7_dist, aes(x = Window.Start, y = Dist), color = "palegreen2") +
  geom_line(data = L8_dist, aes(x = Window.Start, y = Dist), color = "palegreen1") +
  ylim(0, 0.7) +
  xlab("L") +
  ylab("Distance to IQTV")

ggplot() +
  geom_line(data = G1a_dist, aes(x = Window.Start, y = Dist), color = "royalblue4") +
  geom_line(data = G1b_dist, aes(x = Window.Start, y = Dist), color = "royalblue3") +
  geom_line(data = G1c_dist, aes(x = Window.Start, y = Dist), color = "royalblue") +
  geom_line(data = G1d_dist, aes(x = Window.Start, y = Dist), color = "royalblue1") +
  geom_line(data = G2a_dist, aes(x = Window.Start, y = Dist), color = "darkred") +
  geom_line(data = G2b_dist, aes(x = Window.Start, y = Dist), color = "firebrick") +
  geom_line(data = G2c_dist, aes(x = Window.Start, y = Dist), color = "firebrick3") +
  geom_line(data = G2d_dist, aes(x = Window.Start, y = Dist), color = "firebrick1") +
  ylim(0, 0.7) +
  xlab("M") +
  ylab("Distance to IQTV")

ggplot() +
  geom_line(data = S1_dist, aes(x = Window.Start, y = Dist), color = "goldenrod4") +
  geom_line(data = S2_dist, aes(x = Window.Start, y = Dist), color = "goldenrod3") +
  geom_line(data = S3_dist, aes(x = Window.Start, y = Dist), color = "goldenrod2") +
  geom_line(data = S4_dist, aes(x = Window.Start, y = Dist), color = "goldenrod1") +
  geom_line(data = S5_dist, aes(x = Window.Start, y = Dist), color = "lightgoldenrod4") +
  geom_line(data = S6_dist, aes(x = Window.Start, y = Dist), color = "lightgoldenrod3") +
  geom_line(data = S7_dist, aes(x = Window.Start, y = Dist), color = "lightgoldenrod2") +
  geom_line(data = S8_dist, aes(x = Window.Start, y = Dist), color = "lightgoldenrod1") +
  ylim(0, 0.7) +
  xlab("S") +
  ylab("Distance to IQTV")

#Repeat analyses by comparing the M segment with pseudo-outgroup sequences
#within each clade
#G1 outgroup: Trinidad 1955, G2 outgroup: Panama 1989
#Read sequences
G1INa.seqs <- fasta2DNAbin('OROV_M_G1_BRA1960.fasta')
G1INb.seqs <- fasta2DNAbin('OROV_M_G1_BRA2009.fasta')
G1INc.seqs <- fasta2DNAbin('OROV_M_G1_PAN1999.fasta')
G1OUTa.seqs <- fasta2DNAbin('OROV_M_G1_BRA1960out.fasta')
G1OUTb.seqs <- fasta2DNAbin('OROV_M_G1_BRA2009out.fasta')
G1OUTc.seqs <- fasta2DNAbin('OROV_M_G1_PAN1999out.fasta')

G2INa.seqs <- fasta2DNAbin('OROV_M_G2_ECU2016.fasta')
G2INb.seqs <- fasta2DNAbin('OROV_M_G2_PER2008.fasta')
G2INc.seqs <- fasta2DNAbin('OROV_M_G2_PER2000.fasta')
G2OUTa.seqs <- fasta2DNAbin('OROV_M_G2_ECU2016out.fasta')
G2OUTb.seqs <- fasta2DNAbin('OROV_M_G2_PER2008out.fasta')
G2OUTc.seqs <- fasta2DNAbin('OROV_M_G2_PER2000out.fasta')

#Set sliding windows on alignments to estimate genetic distances
G1INa_wndws <- slidingWindow(G1INa.seqs, width = round((length(G1INa.seqs)/2)*0.07), interval = 20)
length(G1INa_wndws)
G1INb_wndws <- slidingWindow(G1INb.seqs, width = round((length(G1INb.seqs)/2)*0.07), interval = 20)
length(G1INb_wndws)
G1INc_wndws <- slidingWindow(G1INc.seqs, width = round((length(G1INc.seqs)/2)*0.07), interval = 20)
length(G1INc_wndws)
G1OUTa_wndws <- slidingWindow(G1OUTa.seqs, width = round((length(G1OUTa.seqs)/2)*0.07), interval = 20)
length(G1OUTa_wndws)
G1OUTb_wndws <- slidingWindow(G1OUTb.seqs, width = round((length(G1OUTb.seqs)/2)*0.07), interval = 20)
length(G1OUTb_wndws)
G1OUTc_wndws <- slidingWindow(G1OUTc.seqs, width = round((length(G1OUTc.seqs)/2)*0.07), interval = 20)
length(G1OUTc_wndws)
G2INa_wndws <- slidingWindow(G2INa.seqs, width = round((length(G2INa.seqs)/2)*0.07), interval = 20)
length(G2INa_wndws)
G2INb_wndws <- slidingWindow(G2INb.seqs, width = round((length(G2INb.seqs)/2)*0.07), interval = 20)
length(G2INb_wndws)
G2INc_wndws <- slidingWindow(G2INc.seqs, width = round((length(G2INc.seqs)/2)*0.07), interval = 20)
length(G2INc_wndws)
G2OUTa_wndws <- slidingWindow(G2OUTa.seqs, width = round((length(G2OUTa.seqs)/2)*0.07), interval = 20)
length(G2OUTa_wndws)
G2OUTb_wndws <- slidingWindow(G2OUTb.seqs, width = round((length(G2OUTb.seqs)/2)*0.07), interval = 20)
length(G2OUTb_wndws)
G2OUTc_wndws <- slidingWindow(G2OUTc.seqs, width = round((length(G2OUTc.seqs)/2)*0.07), interval = 20)
length(G2OUTc_wndws)

#Estimate mean genetic distance per window
G1INa_dist <- sapply(G1INa_wndws, dist.dna)
x <- length(G1INa_wndws)
G1INa_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1INa_dist
)
G1INb_dist <- sapply(G1INb_wndws, dist.dna)
x <- length(G1INb_wndws)
G1INb_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1INb_dist
)
G1INc_dist <- sapply(G1INc_wndws, dist.dna)
x <- length(G1INc_wndws)
G1INc_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1INc_dist
)
G1OUTa_dist <- sapply(G1OUTa_wndws, dist.dna)
x <- length(G1OUTa_wndws)
G1OUTa_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1OUTa_dist
)
G1OUTb_dist <- sapply(G1OUTb_wndws, dist.dna)
x <- length(G1OUTb_wndws)
G1OUTb_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1OUTb_dist
)
G1OUTc_dist <- sapply(G1OUTc_wndws, dist.dna)
x <- length(G1OUTc_wndws)
G1OUTc_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G1OUTc_dist
)
G2INa_dist <- sapply(G2INa_wndws, dist.dna)
x <- length(G2INa_wndws)
G2INa_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2INa_dist
)
G2INb_dist <- sapply(G2INb_wndws, dist.dna)
x <- length(G2INb_wndws)
G2INb_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2INb_dist
)
G2INc_dist <- sapply(G2INc_wndws, dist.dna)
x <- length(G2INc_wndws)
G2INc_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2INc_dist
)
G2OUTa_dist <- sapply(G2OUTa_wndws, dist.dna)
x <- length(G2OUTa_wndws)
G2OUTa_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2OUTa_dist
)
G2OUTb_dist <- sapply(G2OUTb_wndws, dist.dna)
x <- length(G2OUTb_wndws)
G2OUTb_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2OUTb_dist
)
G2OUTc_dist <- sapply(G2OUTc_wndws, dist.dna)
x <- length(G2OUTc_wndws)
G2OUTc_dist <- data.frame(
  Window.Start = c(seq(1, 20*x, by = 20)),
  Dist = G2OUTc_dist
)

#SimPlots
ggplot() +
  geom_line(data = G1INa_dist, aes(x = Window.Start, y = Dist), color = "red") +
  geom_line(data = G1INb_dist, aes(x = Window.Start, y = Dist), color = "red") +
  geom_line(data = G1INc_dist, aes(x = Window.Start, y = Dist), color = "red") +
  geom_line(data = G1OUTa_dist, aes(x = Window.Start, y = Dist), color = "darkred") +
  geom_line(data = G1OUTb_dist, aes(x = Window.Start, y = Dist), color = "darkred") +
  geom_line(data = G1OUTc_dist, aes(x = Window.Start, y = Dist), color = "darkred") +
  geom_line(data = G2INa_dist, aes(x = Window.Start, y = Dist), color = "blue") +
  geom_line(data = G2INb_dist, aes(x = Window.Start, y = Dist), color = "blue") +
  geom_line(data = G2INc_dist, aes(x = Window.Start, y = Dist), color = "blue") +
  geom_line(data = G2OUTa_dist, aes(x = Window.Start, y = Dist), color = "midnightblue") +
  geom_line(data = G2OUTb_dist, aes(x = Window.Start, y = Dist), color = "midnightblue") +
  geom_line(data = G2OUTc_dist, aes(x = Window.Start, y = Dist), color = "midnightblue") +
  ylim(0, 0.35) +
  xlab("M") +
  ylab("Distance to basal group")

#Estimate mean relative genetic distances from each genotype to IQTV outgroup
G1_dist <- apply(cbind(G1a_dist$Dist,G1b_dist$Dist,G1c_dist$Dist),
                 1, mean)
G1_dist <- data.frame(Window = G1a_dist$Window,
                      Dist = G1_dist)
G2_dist <- apply(cbind(G2a_dist$Dist,G2b_dist$Dist,G2c_dist$Dist),
                 1, mean)
G2_dist <- data.frame(Window = G2a_dist$Window,
                      Dist = G2_dist)

Ref1 <- data.frame(
  Window = G1_dist$Window,
  Dist = G1_dist$Dist/(G1_dist$Dist+G2_dist$Dist)
)

Ref2 <- data.frame(
  Window = G2_dist$Window,
  Dist = G2_dist$Dist/(G1_dist$Dist+G2_dist$Dist)
)

ggplot() +
  geom_line(data = Ref1, aes(x = Ref1$Window, y = Ref1$Dist),color = "royalblue4") +
  geom_line(data = Ref2, aes(x = Ref2$Window, y = Ref2$Dist),color = "darkred") +
  ylim(0.4, 0.6) +
  xlab("M") +
  ylab("Relative distance to IQTV")

#Find breakpoints
Break1 <- data.frame(
  Window = Ref1$Window,
  Dist = Ref2$Dist-Ref1$Dist
)
Break2 <- data.frame(
  Window = Break1$Window,
  Dist = Break1$Dist+abs(Break1$Dist)
)

SrcGen <- data.frame(
  Window.Start.Site = Break2$Window,
  Window.End.Site = Break2$Window+round((length(G2a.seqs)/2)*0.07),
  Lineage = ifelse(Break2$Dist > 0, "LineageA", "LineageB")
)

SrcGen
