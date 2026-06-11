#### 数据导入_At部分####
setwd("R_downstrem")
level3 <- read.csv("level-3.csv", row.names = 1)
n1 <- ncol(level3);n1
level3 <- level3[,1:(n1-2)]
level3.t <- as.data.frame(t(level3))
level3.t$sum <- apply(level3.t, 1, sum)
level3.a <- read.csv("level-3.csv", header = F, row.names = 1)
#level3.a <- level3.a[,-1]
n2 <- ncol(level3.a);n2
name <- unlist(level3.a[1,1:(n2-2)])
level3.t$name <- name 
n3 <- ncol(level3.t)
level3.t <- level3.t[,c(n3,n3-1,1:(n3-2))]
library(tidyr)
level3.t <- separate(level3.t, name, into = c("K", "P", "C"), sep = ";")
level3.t$K <- gsub(".*__", "", level3.t$K )
level3.t$P <- gsub(".*__", "", level3.t$P)
level3.t$C <- gsub(".*__", "", level3.t$C)
library(ggplot2)
p <- aggregate(level3.t[,-1:-3], by = list(level3.t$P), sum)
p <- p[order(-p$sum),]
names(p)[1] <- "phylum"
c <- aggregate(level3.t[,-1:-3], by = list(level3.t$C), sum)
c <- c[order(-c$sum),]
names(c)[1] <- "class"
proteobacteria <- subset(level3.t,P == "Proteobacteria")
proteobacteria <- proteobacteria[order(-proteobacteria$sum),]
proteobacteria <- proteobacteria[,-1:-2]
names(proteobacteria)[1] <- "phylum"
p_no_proteobacteria <- subset(p, p$phylum != "Proteobacteria")
match(names(proteobacteria),names(p_no_proteobacteria))
phylum <- rbind(proteobacteria, p_no_proteobacteria)
phylum <- phylum[order(-phylum$sum),]

dim(phylum)
phylum_top10 <- phylum[1:9,]
row.names(phylum_top10) <- phylum_top10$phylum
phylum_top10 <- phylum_top10[,-1:-2]
n4 <- nrow(phylum)
phylum_top10[10,] <- apply(phylum[10:n4,-1:-2], 2, sum)
row.names(phylum_top10)[10] <- "Others"

phylum_top20 <- phylum[1:28,]
row.names(phylum_top20) <- c(1:28)
phylum_top20 <- phylum_top20[-c(14,18,19,20,21,23,24,26,28),]
row.names(phylum_top20) <- phylum_top20$phylum
phylum_top20 <- phylum_top20[,-1:-2]
n4 <- nrow(phylum)
phylum_top20[20,] <- apply(phylum[29:n4,-1:-2], 2, sum)
row.names(phylum_top20)[20] <- "Others"

#write.csv(phylum_top10, "phylum_top10.csv")
##################突变体部分的数据#############################
metadata <- read.csv("At_20240729_16S_metadata_m.csv", header = T,row.names = 1)
metadata$Group <- paste(metadata$Genotype, metadata$Niche, sep = "-")
#phylum_top10 <- read.csv("phylum_top10.csv", row.names = 1)
phylum_top10$tax <- row.names(phylum_top10)

library(reshape2)
phylum_top10.melt <- as.data.frame(melt(phylum_top10, id.vars = "tax"))

phylum_top10.melt.m <- merge(phylum_top10.melt, metadata, 
                             by.x = "variable", by.y = "row.names", all.y  = T)
phylum_top10.melt.m$Niche <- factor(phylum_top10.melt.m$Niche, 
                                    levels = c("Rhizosphere", "Root endo"))
View(phylum_top10.melt.m)


color <- c("#EE766D", "#7CAE00", "#00BFC4", "#C77CFF",  "#00BE67",
           "#00A9FF", "#FF61CC","#FF8808","#FF4500","#CD9600")


bp1 <- ggplot(phylum_top10.melt.m,
              aes(x = variable, y = value, fill = tax )) +
  geom_bar(stat = "identity", position = "fill", width = 1)+
  labs(x = (""),
       y = ("percentage (%)"),
       title = ("Phylum distribution")) + 
  scale_fill_manual(values = color) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid( ~ Niche + Genotype, scales = "free_x", 
              space = "free", switch = "x") +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, #使标题居中
                                  size = 16,color="black", face = "bold"), 
        legend.position = "right",
        legend.title = element_text(size = 20,color = "black", face = "bold"),
        legend.text = element_text(size = 20,color = "black"),
        axis.text = element_text(size = 30,color = "black"),
        axis.text.x = element_text(size= 22,angle = 45,vjust = 0.25, hjust = 0),
        axis.title = element_text(size = 30,color="black", face= "plain"));bp1
########################用mean去做###############
#phylum_top10.melt.m.ex <- phylum_top10.melt.m.ex[c(-61:-70,-191:-200,-261:-270,-461:-470),]
phylum_top10.melt.m.ex <- phylum_top10.melt.m 
df <- aggregate(phylum_top10.melt.m.ex[,1:5],
                   by=list(
                     Genotype=phylum_top10.melt.m.ex$Genotype,Niche=phylum_top10.melt.m.ex$Niche,tax=phylum_top10.melt.m.ex$tax),FUN=mean)
df <- df[,c(1:3,6)]
df <- df[order(df$Genotype),]
df <- df[order(df$Niche),]
df$Genotype <- gsub("Col0","Col-0",df$Genotype)
#df$Genotype <- gsub("Dm","p455cam7",df$Genotype)
#df$Genotype <- gsub("Cam7","cam7",df$Genotype)
#把genotype改成factor，下面是修改的代码
#lvls <- unique(df$Genotype)
aa <- c("Col-0","p180","p455","p444","p328")
df$Genotype <- factor(df$Genotype,levels = aa )
#df <- df[,c(1:3,6)]
bp <- ggplot(df,
                aes(x = Genotype, y = value, fill = tax )) +
  geom_bar(stat = "identity", position = "fill", width = 0.75)+
  labs(x = (""),
       y = ("percentage (%)"),
       title = ("Phylum distribution")) + 
  scale_fill_manual(values = color) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid( ~ Niche, scales = "free_x", 
              space = "free", switch = "x") +
  theme_bw() + 
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, #使标题居中
                                  size = 30,color="black",face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 20,color = "black", face = "bold"),
        legend.text = element_text(size = 20,color = "black"),
        axis.text = element_text(size = 30,color = "black"),
        axis.text.x = element_text(size= 22,angle = 45,vjust = 0.25, hjust = 0),
        axis.title = element_text(size = 30,color="black", face= "plain"));bp
