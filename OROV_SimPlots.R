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

C1.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_BRA1960.fasta')
C2.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_BRA2000.fasta')
C3.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_BRA2009.fasta')
C4.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_PAN1999.fasta')
C5.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_ECU2016.fasta')
C6.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_PAN1989.fasta')
C7.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_PER2000.fasta')
C8.seqs <- fasta2DNAbin('OROV_v_IQTV_Concat_PER2008.fasta')

#Set sliding windows on alignments to estimate genetic distances
L1_wndws <- slidingWindow(L1.seqs, width = round(length(L1.seqs)*0.07), interval = 20)
length(L1_wndws)
L2_wndws <- slidingWindow(L2.seqs, width = round(length(L2.seqs)*0.07), interval = 20)
length(L2_wndws)
L3_wndws <- slidingWindow(L3.seqs, width = round(length(L3.seqs)*0.07), interval = 20)
length(L3_wndws)
L4_wndws <- slidingWindow(L4.seqs, width = round(length(L4.seqs)*0.07), interval = 20)
length(L4_wndws)
L5_wndws <- slidingWindow(L5.seqs, width = round(length(L5.seqs)*0.07), interval = 20)
length(L5_wndws)
L6_wndws <- slidingWindow(L6.seqs, width = round(length(L6.seqs)*0.07), interval = 20)
length(L6_wndws)
L7_wndws <- slidingWindow(L7.seqs, width = round(length(L7.seqs)*0.07), interval = 20)
length(L7_wndws)
L8_wndws <- slidingWindow(L8.seqs, width = round(length(L8.seqs)*0.07), interval = 20)
length(L8_wndws)

G1a_wndws <- slidingWindow(G1a.seqs, width = round(length(G1a.seqs)*0.07), interval = 20)
length(G1a_wndws)
G1b_wndws <- slidingWindow(G1b.seqs, width = round(length(G1b.seqs)*0.07), interval = 20)
length(G1b_wndws)
G1c_wndws <- slidingWindow(G1c.seqs, width = round(length(G1c.seqs)*0.07), interval = 20)
length(G1c_wndws)
G1d_wndws <- slidingWindow(G1d.seqs, width = round(length(G1d.seqs)*0.07), interval = 20)
length(G1d_wndws)
G2a_wndws <- slidingWindow(G2a.seqs, width = round(length(G2a.seqs)*0.07), interval = 20)
length(G2a_wndws)
G2b_wndws <- slidingWindow(G2b.seqs, width = round(length(G2b.seqs)*0.07), interval = 20)
length(G2b_wndws)
G2c_wndws <- slidingWindow(G2c.seqs, width = round(length(G2c.seqs)*0.07), interval = 20)
length(G2c_wndws)
G2d_wndws <- slidingWindow(G2d.seqs, width = round(length(G2d.seqs)*0.07), interval = 20)
length(G2d_wndws)

S1_wndws <- slidingWindow(S1.seqs, width = round(length(S1.seqs)*0.07), interval = 20)
length(S1_wndws)
S2_wndws <- slidingWindow(S2.seqs, width = round(length(S2.seqs)*0.07), interval = 20)
length(S2_wndws)
S3_wndws <- slidingWindow(S3.seqs, width = round(length(S3.seqs)*0.07), interval = 20)
length(S3_wndws)
S4_wndws <- slidingWindow(S4.seqs, width = round(length(S4.seqs)*0.07), interval = 20)
length(S4_wndws)
S5_wndws <- slidingWindow(S5.seqs, width = round(length(S5.seqs)*0.07), interval = 20)
length(S5_wndws)
S6_wndws <- slidingWindow(S6.seqs, width = round(length(S6.seqs)*0.07), interval = 20)
length(S6_wndws)
S7_wndws <- slidingWindow(S7.seqs, width = round(length(S7.seqs)*0.07), interval = 20)
length(S7_wndws)
S8_wndws <- slidingWindow(S8.seqs, width = round(length(S8.seqs)*0.07), interval = 20)
length(S8_wndws)

C1_wndws <- slidingWindow(C1.seqs, width = 600, interval = 20)
length(C1_wndws)
C2_wndws <- slidingWindow(C2.seqs, width = 600, interval = 20)
length(C2_wndws)
C3_wndws <- slidingWindow(C3.seqs, width = 600, interval = 20)
length(C3_wndws)
C4_wndws <- slidingWindow(C4.seqs, width = 600, interval = 20)
length(C4_wndws)
C5_wndws <- slidingWindow(C5.seqs, width = 600, interval = 20)
length(C5_wndws)
C6_wndws <- slidingWindow(C6.seqs, width = 600, interval = 20)
length(C6_wndws)
C7_wndws <- slidingWindow(C7.seqs, width = 600, interval = 20)
length(C7_wndws)
C8_wndws <- slidingWindow(C8.seqs, width = 600, interval = 20)
length(C8_wndws)

#Estimate mean genetic distance per window
L1_dist <- sapply(L1_wndws, dist.dna)
x <- length(L1_wndws)
L1_dist <- data.frame(
  Window = c(1:x),
  Dist = L1_dist
)
L2_dist <- sapply(L2_wndws, dist.dna)
x <- length(L2_wndws)
L2_dist <- data.frame(
  Window = c(1:x),
  Dist = L2_dist
)
L3_dist <- sapply(L3_wndws, dist.dna)
x <- length(L3_wndws)
L3_dist <- data.frame(
  Window = c(1:x),
  Dist = L3_dist
)
L4_dist <- sapply(L4_wndws, dist.dna)
x <- length(L4_wndws)
L4_dist <- data.frame(
  Window = c(1:x),
  Dist = L4_dist
)
L5_dist <- sapply(L5_wndws, dist.dna)
x <- length(L5_wndws)
L5_dist <- data.frame(
  Window = c(1:x),
  Dist = L5_dist
)
L6_dist <- sapply(L6_wndws, dist.dna)
x <- length(L6_wndws)
L6_dist <- data.frame(
  Window = c(1:x),
  Dist = L6_dist
)
L7_dist <- sapply(L7_wndws, dist.dna)
x <- length(L7_wndws)
L7_dist <- data.frame(
  Window = c(1:x),
  Dist = L7_dist
)
L8_dist <- sapply(L8_wndws, dist.dna)
x <- length(L8_wndws)
L8_dist <- data.frame(
  Window = c(1:x),
  Dist = L8_dist
)

G1a_dist <- sapply(G1a_wndws, dist.dna)
x <- length(G1a_wndws)
G1a_dist <- data.frame(
  Window = c(1:x),
  Dist = G1a_dist
)
G1b_dist <- sapply(G1b_wndws, dist.dna)
x <- length(G1b_wndws)
G1b_dist <- data.frame(
  Window = c(1:x),
  Dist = G1b_dist
)
G1c_dist <- sapply(G1c_wndws, dist.dna)
x <- length(G1c_wndws)
G1c_dist <- data.frame(
  Window = c(1:x),
  Dist = G1c_dist
)
G1d_dist <- sapply(G1d_wndws, dist.dna)
x <- length(G1d_wndws)
G1d_dist <- data.frame(
  Window = c(1:x),
  Dist = G1d_dist
)
G2a_dist <- sapply(G2a_wndws, dist.dna)
x <- length(G2a_wndws)
G2a_dist <- data.frame(
  Window = c(1:x),
  Dist = G2a_dist
)
G2b_dist <- sapply(G2b_wndws, dist.dna)
x <- length(G2b_wndws)
G2b_dist <- data.frame(
  Window = c(1:x),
  Dist = G2b_dist
)
G2c_dist <- sapply(G2c_wndws, dist.dna)
x <- length(G2c_wndws)
G2c_dist <- data.frame(
  Window = c(1:x),
  Dist = G2c_dist
)
G2d_dist <- sapply(G2d_wndws, dist.dna)
x <- length(G2d_wndws)
G2d_dist <- data.frame(
  Window = c(1:x),
  Dist = G2d_dist
)

S1_dist <- sapply(S1_wndws, dist.dna)
x <- length(S1_wndws)
S1_dist <- data.frame(
  Window = c(1:x),
  Dist = S1_dist
)
S2_dist <- sapply(S2_wndws, dist.dna)
x <- length(S2_wndws)
S2_dist <- data.frame(
  Window = c(1:x),
  Dist = S2_dist
)
S3_dist <- sapply(S3_wndws, dist.dna)
x <- length(S3_wndws)
S3_dist <- data.frame(
  Window = c(1:x),
  Dist = S3_dist
)
S4_dist <- sapply(S4_wndws, dist.dna)
x <- length(S4_wndws)
S4_dist <- data.frame(
  Window = c(1:x),
  Dist = S4_dist
)
S5_dist <- sapply(S5_wndws, dist.dna)
x <- length(S5_wndws)
S5_dist <- data.frame(
  Window = c(1:x),
  Dist = S5_dist
)
S6_dist <- sapply(S6_wndws, dist.dna)
x <- length(S6_wndws)
S6_dist <- data.frame(
  Window = c(1:x),
  Dist = S6_dist
)
S7_dist <- sapply(S7_wndws, dist.dna)
x <- length(S7_wndws)
S7_dist <- data.frame(
  Window = c(1:x),
  Dist = S7_dist
)
S8_dist <- sapply(S8_wndws, dist.dna)
x <- length(S8_wndws)
S8_dist <- data.frame(
  Window = c(1:x),
  Dist = S8_dist
)

C1_dist <- sapply(C1_wndws, dist.dna)
x <- length(C1_wndws)
C1_dist <- data.frame(
  Window = c(1:x),
  Dist = C1_dist
)
C2_dist <- sapply(C2_wndws, dist.dna)
x <- length(C2_wndws)
C2_dist <- data.frame(
  Window = c(1:x),
  Dist = C2_dist
)
C3_dist <- sapply(C3_wndws, dist.dna)
x <- length(C3_wndws)
C3_dist <- data.frame(
  Window = c(1:x),
  Dist = C3_dist
)
C4_dist <- sapply(C4_wndws, dist.dna)
x <- length(C4_wndws)
C4_dist <- data.frame(
  Window = c(1:x),
  Dist = C4_dist
)
C5_dist <- sapply(C5_wndws, dist.dna)
x <- length(C5_wndws)
C5_dist <- data.frame(
  Window = c(1:x),
  Dist = C5_dist
)
C6_dist <- sapply(C6_wndws, dist.dna)
x <- length(C6_wndws)
C6_dist <- data.frame(
  Window = c(1:x),
  Dist = C6_dist
)
C7_dist <- sapply(C7_wndws, dist.dna)
x <- length(C7_wndws)
C7_dist <- data.frame(
  Window = c(1:x),
  Dist = C7_dist
)
C8_dist <- sapply(C8_wndws, dist.dna)
x <- length(C8_wndws)
C8_dist <- data.frame(
  Window = c(1:x),
  Dist = C8_dist
)

#Plot genetic distances for each segment
ggplot() +
  geom_line(data = L1_dist, aes(x = Window, y = Dist), color = "red") +
  geom_line(data = L2_dist, aes(x = Window, y = Dist), color = "blue") +
  geom_line(data = L3_dist, aes(x = Window, y = Dist), color = "yellow") +
  geom_line(data = L4_dist, aes(x = Window, y = Dist), color = "purple") +
  geom_line(data = L5_dist, aes(x = Window, y = Dist), color = "orange") +
  geom_line(data = L6_dist, aes(x = Window, y = Dist), color = "green") +
  geom_line(data = L7_dist, aes(x = Window, y = Dist), color = "cyan") +
  geom_line(data = L8_dist, aes(x = Window, y = Dist), color = "orchid") +
  ylim(0:1) +
  xlab("L") +
  ylab("Distance to IQTV")

ggplot() +
  geom_line(data = G1a_dist, aes(x = Window, y = Dist), color = "red") +
  geom_line(data = G1b_dist, aes(x = Window, y = Dist), color = "blue") +
  geom_line(data = G1c_dist, aes(x = Window, y = Dist), color = "yellow") +
  geom_line(data = G1d_dist, aes(x = Window, y = Dist), color = "purple") +
  geom_line(data = G2a_dist, aes(x = Window, y = Dist), color = "orange") +
  geom_line(data = G2b_dist, aes(x = Window, y = Dist), color = "green") +
  geom_line(data = G2c_dist, aes(x = Window, y = Dist), color = "cyan") +
  geom_line(data = G2d_dist, aes(x = Window, y = Dist), color = "orchid") +
  ylim(0:1) +
  xlab("M") +
  ylab("Distance to IQTV")

ggplot() +
  geom_line(data = S1_dist, aes(x = Window, y = Dist), color = "red") +
  geom_line(data = S2_dist, aes(x = Window, y = Dist), color = "blue") +
  geom_line(data = S3_dist, aes(x = Window, y = Dist), color = "yellow") +
  geom_line(data = S4_dist, aes(x = Window, y = Dist), color = "purple") +
  geom_line(data = S5_dist, aes(x = Window, y = Dist), color = "orange") +
  geom_line(data = S6_dist, aes(x = Window, y = Dist), color = "green") +
  geom_line(data = S7_dist, aes(x = Window, y = Dist), color = "cyan") +
  geom_line(data = S8_dist, aes(x = Window, y = Dist), color = "orchid") +
  ylim(0:1) +
  xlab("S") +
  ylab("Distance to IQTV")

ggplot() +
  geom_line(data = C1_dist, aes(x = Window, y = Dist), color = "red") +
  geom_line(data = C2_dist, aes(x = Window, y = Dist), color = "blue") +
  geom_line(data = C3_dist, aes(x = Window, y = Dist), color = "yellow") +
  geom_line(data = C4_dist, aes(x = Window, y = Dist), color = "purple") +
  geom_line(data = C5_dist, aes(x = Window, y = Dist), color = "orange") +
  geom_line(data = C6_dist, aes(x = Window, y = Dist), color = "green") +
  geom_line(data = C7_dist, aes(x = Window, y = Dist), color = "cyan") +
  geom_line(data = C8_dist, aes(x = Window, y = Dist), color = "orchid") +
  ylim(0:1) +
  xlab("Concatenated genome (ORFs)") +
  ylab("Distance to IQTV")