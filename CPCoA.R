getwd()
setwd("PCoA")
At_aha_metadata <- read.csv(file = "At_20240729_16S_metadata_aha.csv", header = T, row.names = 1)
At_feature_table <- read.csv(file = "At_dada2_table_with_Bacteria_and_Archaea_no_Chloroplast_from_biom.csv", header = T)
row.names(At_feature_table) <- At_feature_table[,1]
At_feature_table <- At_feature_table[,-1]
###只用RS
At_feature_table_aha <- At_feature_table[,c(65:73,85:95,106:119,131:139)]
At_feature_table_aha <- as.data.frame(lapply(At_feature_table_aha,as.numeric))
At_feature_table_aha$sum <- apply(At_feature_table_aha,1,sum)
At_feature_table_aha <- At_feature_table_aha[order(-At_feature_table_aha$sum),]
#colnames(At_feature_table)[1] <- 'Feature_ID'
At_otu_aha <- rep(c("OTU"),53688)
At_number_aha <- c(1:53688)
At_id_aha <- paste(At_otu_aha,At_number_aha,sep = "")
At_feature_table_aha <- cbind(At_id_aha,At_feature_table_aha)
colnames(At_feature_table_aha)[1] <- 'OTU_ID'
row.names(At_feature_table_aha) <- At_feature_table_aha[,1]
#At_feature_table_aha <- At_feature_table_aha[,-1]
#At_feature_table <- At_feature_table[,-2]
At_count_matrix_aha <- as.data.frame(t(At_feature_table_aha[,2:45]))
At_metadata_rhizosphere_aha <- subset(At_aha_metadata, Niche == "Rhizosphere")
At_feature_table_rhizosphere_aha <- t(At_feature_table_aha[,2:45])[row.names(At_count_matrix_aha) %in% row.names(At_metadata_rhizosphere_aha),]
library(vegan)
capscale.feature_table_At_rhizosphere_aha <- capscale(At_feature_table_rhizosphere_aha ~ Genotype, data = At_metadata_rhizosphere_aha, add=F, dist="bray")
R2_capscale.feature_table_At_rhizosphere_aha <- RsquareAdj(capscale.feature_table_At_rhizosphere_aha)$r.squared
adjR2_capscale.feature_table_At_rhizosphere_aha <- RsquareAdj(capscale.feature_table_At_rhizosphere_aha)$adj.r.squa
capscale.feature_table_At_rhizosphere_aha.anova <- anova.cca(capscale.feature_table_At_rhizosphere_aha, permutations = 10000)
capscale.feature_table_At_rhizosphere_aha.p.val <- capscale.feature_table_At_rhizosphere_aha.anova[1, 4]
wa.capscale.feature_table_At_rhizosphere_aha <- capscale.feature_table_At_rhizosphere_aha$CCA$wa
cap.prop.feature_table_At_rhizosphere_aha <- capscale.feature_table_At_rhizosphere_aha$CCA$eig/sum(capscale.feature_table_At_rhizosphere_aha$CCA$eig)
metadata_feature_table_At_rhizosphere_aha <- cbind(At_metadata_rhizosphere_aha, wa.capscale.feature_table_At_rhizosphere_aha)
metadata_feature_table_At_rhizosphere_aha$Genotype <- factor(metadata_feature_table_At_rhizosphere_aha$Genotype,levels=c("Col0","aha1_1","aha1_2","p455"))
At_col_Genotype_aha <- c("#BC7FF7","#56B4E9", "#009E73", "#F0E442")
metadata_feature_table_At_rhizosphere_aha$Genotype <- gsub("Col0","Col-0",metadata_feature_table_At_rhizosphere_aha$Genotype)
metadata_feature_table_At_rhizosphere_aha$Genotype <- gsub("p455","raps1",metadata_feature_table_At_rhizosphere_aha$Genotype)
library(ggplot2)
p_rarefied_At_aha = ggplot(metadata_feature_table_At_rhizosphere_aha, aes(CAP1, CAP2, color=Genotype)) +
  geom_point(cex=10) +
  scale_color_manual(values=At_col_Genotype_aha)  +
  theme_bw() + 
  labs(x=paste("PCoA 1 (", format(100 * cap.prop.feature_table_At_rhizosphere_aha[1], digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * cap.prop.feature_table_At_rhizosphere_aha[2], digits=4), "%)", sep="")) +
  ggtitle(paste(format(100 * R2_capscale.feature_table_At_rhizosphere_aha, digits=3)," % of variance; p=",format(capscale.feature_table_At_rhizosphere_aha.p.val,digits=2),sep=""))+
  theme(legend.position= "right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text = element_text(size = 10, colour = "black"))
p_rarefied_At_aha = p_rarefied_At_aha + stat_ellipse(level=0.68);p_rarefied_At_aha
