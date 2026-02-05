library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
# ----------------------------------------------------------------------------------------- #
#                                             figure2. A                                    #
# ----------------------------------------------------------------------------------------- #
gene_transcription_error <- gene_transcription_error[gene_transcription_error$Error.rate!=0,]
gene_transcription_error$type <- "Mistranscription"
gene_translation_error$type <- "Mistranslation"
persite_data <- data.frame(type=c("Mistranscription","Mistranslation"),
                           mean=c(mean(gene_transcription_error$Error.rate),mean(gene_translation_error$Error_Rate)),
                           mean_log=c(mean(gene_transcription_error$log),mean(gene_translation_error$log)),
                           sd=c(sd(gene_transcription_error$Error.rate),sd(gene_translation_error$Error_Rate)),
                           sd_log = c(sd(gene_transcription_error$log),sd(gene_translation_error$log)),
                           se = c(sd(gene_transcription_error$Error.rate)/sqrt(length(gene_transcription_error$Error.rate)),sd(gene_translation_error$Error_Rate)/sqrt(length(gene_translation_error$Error_Rate))),
                           se_log = c(sd(gene_transcription_error$log)/sqrt(length(gene_transcription_error$log)),sd(gene_translation_error$log)/sqrt(length(gene_translation_error$log))))
test <- wilcox.test(gene_transcription_error$Error.rate,gene_translation_error$Error_Rate)
test$p.value
pdf("persite.pdf");
ggplot(persite_data) +
  geom_bar(aes(x = type, y = mean_log), stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("gray15", "gray20")) +
  geom_errorbar(aes(x = type, ymin = mean_log - sd_log, ymax = mean_log + sd_log), width = 0.4, colour = "black", alpha = 0.9, linewidth = 1.3) +
  labs(x = "", y = "", title = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 23,face = "bold"),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 23),
        plot.title = element_text(size = 15, hjust = 0.5, vjust = 2.5),
        panel.border = element_blank(),
        axis.line = element_line(size = 1)) +
  scale_y_log10()  
dev.off()

# ----------------------------------------------------------------------------------------- #
#                                             figure2. B                                    #
# ----------------------------------------------------------------------------------------- #
pergene_TR <- pergene_TR[pergene_TR$Error_rate!=0,]
pergene_TR$err1 <- 1-(1-pergene_TR$Error_rate)^pergene_TR$Lnt
pergene_TL$err2 <- 1-(1-pergene_TL$Error_rate)^pergene_TL$Lp

pergene_data <- data.frame(type=c("Mistranscription","Mistranslation"),
                           mean=c(mean(pergene_TR$err1),mean(pergene_TL$err2)),
                           mean_log=c(mean(pergene_TR$log),mean(pergene_TL$log)),
                           sd=c(sd(pergene_TR$err1),sd(pergene_TL$err2)),
                           sd_log = c(sd(pergene_TR$log),sd(pergene_TL$log)),
                           se = c(sd(pergene_TR$err1)/sqrt(length(pergene_TR$err1)),sd(pergene_TL$err2)/sqrt(length(pergene_TL$err2))),
                           se_log = c(sd(pergene_TR$log)/sqrt(length(pergene_TR$log)),sd(pergene_TL$log)/sqrt(length(pergene_TL$log))))
wilcox.test(pergene_TR$err1,pergene_TL$err2)

pdf("pergene.pdf");
ggplot(pergene_data) +
  geom_bar(aes(x = type, y = mean_log), stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("gray15", "gray20")) +
  geom_errorbar(aes(x = type, ymin = mean_log - sd_log, ymax = mean_log + sd_log), width = 0.4, colour = "black", alpha = 0.9, linewidth = 1.3) +
  labs(x = "", y = "", title = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 23,face = "bold"),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 23),
        plot.title = element_text(size = 15, hjust = 0.5, vjust = 2.5),
        panel.border = element_blank(),
        axis.line = element_line(size = 1)) +
  scale_y_log10()  
dev.off()

# ----------------------------------------------------------------------------------------- #
#                                             figure2. C                                    #
# ----------------------------------------------------------------------------------------- #
pergene_TL$SYMBOL <- 0
for(i in seq(1,nrow(pergene_TL)))
{
  if(pergene_TL[i,1] %in% ensembl$protein)
  {
    loc <- which(ensembl$protein == pergene_TL[i,1])
    pergene_TL[i,6] <- ensembl[loc,2]
  }
}
gene <- intersect(pergene_TL$SYMBOL,pergene_TR$gene)
spec_genedata <- data.frame(gene=gene,err1=0,err2=0,frac_TR_TL=0)
for(i in seq(1,length(gene)))
{
  loc_1 <- which(pergene_TR$gene == gene[[i]])
  spec_genedata[i,2] <- pergene_TR[loc_1,4]
  loc_2 <- which(pergene_TL$SYMBOL == gene[[i]])
  spec_genedata[i,3] <- pergene_TL[loc_2,4]
  spec_genedata[i,4] <- as.numeric(spec_genedata[i,2]/spec_genedata[i,3])
}
spec_genedata$type <- ifelse(spec_genedata$frac_TR_TL > 1, "Transcription Error rate > Translation Error rate", "Transcription Error rate < Translation Error rate")

spec_genedata$frac_TR_TL_adjusted <- ifelse(spec_genedata$frac_TR_TL < 1, 
                                            -abs(spec_genedata$frac_TR_TL - 1), 
                                            abs(spec_genedata$frac_TR_TL - 1))

spec_genedata$bar_color <- ifelse(spec_genedata$frac_TR_TL_adjusted > 0, "black", "gray60")

pdf("special_pergene.pdf");
ggplot(spec_genedata, aes(x = reorder(gene, frac_TR_TL_adjusted), y = frac_TR_TL_adjusted, fill = bar_color)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  labs(y = "Ratio between per-gene rates of") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5), 
                     labels = function(x) ifelse(x < 0, 1 - abs(x), 1 + abs(x))) + 
  labs(x = "Gene", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 2.5)
  ) +
  scale_fill_identity()
dev.off()
