library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
# ----------------------------------------------------------------------------------------- #
#                                             figure4. B                                    #
# ----------------------------------------------------------------------------------------- #
fn_double <- fn_fitness[fn_fitness$n_mut<=2,]
fn_select <- data.frame(geno = fn_double$geno,D7_rg1=fn_double$D7_rep1_rg,D7_rg2=fn_double$D7_rep2_rg,D7_rg3=fn_double$D7_rep3_rg,type=fn_double$n_mut)
fn_select <- fn_select[fn_select$type!=0,]
fn_select$AGS_D7_rg <- rowMeans(fn_select[, c("D7_rg1", "D7_rg2", "D7_rg3")], na.rm = TRUE)  
C0T7 <- data.frame(gene=fn_select$geno,COT7_fn=(fn_select$AGS_D7_rg),tag=fn_select$type,gene1=NA,gene2=NA,fitness1=NA,fitness2=NA,fitness12=NA)
single_fitness <- C0T7[C0T7[,3]==1,2]
double_fitness <- C0T7[C0T7[,3]==2,2]

for(i in seq(1,nrow(C0T7)))
{
  if(C0T7[i,3]==1)
  {C0T7[i,4] <- C0T7[i,1]}
  if(C0T7[i,3]==2)
  {
    a <- unlist(strsplit(C0T7[i,1],split = " "))
    C0T7[i,4] <- a[[1]]
    C0T7[i,5] <- a[[2]]
  }
}

for(i in seq(1,nrow(C0T7)))
{
  if(C0T7[i,3]==1)
  {C0T7[i,6] <- C0T7[i,2]}
  if(C0T7[i,3]==2)
  {
    C0T7[i,8] <- C0T7[i,2]
    a <- which(C0T7[,1]==C0T7[i,4])
    if(length(a)>0)
    { C0T7[i,6] <- C0T7[a,2]}
    b <- which(C0T7[,1]==C0T7[i,5])
    if(length(b)>0)
    {C0T7[i,7] <- C0T7[b,2]}
  }
}
fitness_double <- C0T7[C0T7[,3]==2,]
fitness_double <- na.omit(fitness_double)
fitness_double$s1 <- (fitness_double$fitness1-1)
fitness_double$s2 <- (fitness_double$fitness2-1)
fitness_double$s12 <- (fitness_double$fitness12-1)
fitness_double$e_s <- fitness_double$s12-(fitness_double$s1+fitness_double$s2)
tr<- 0.8;tl<- 0.3
fitness_double$tre <- fitness_double$s1*tr*(1-tl)
fitness_double$tle <- fitness_double$s2*tl*(1-tr)
fitness_double$tr_tl <- fitness_double$tre+fitness_double$tle+fitness_double$s12*tr*tl

fitness_double$u1_s <- (fitness_double$s1*tr)+(fitness_double$s2*tl)
fitness_double$u2_s <- (fitness_double$tr_tl)


fitness_double$u2_u1_s <- fitness_double$u2_s-fitness_double$u1_s
summary(fitness_double$u2_u1_s);
summary(fitness_double$e_s)

fitness_double$relative_s <- fitness_double$u2_u1_s/(1+fitness_double$u1_s)
count_fitness <- fitness_double %>%
  group_by(gene1) %>%
  dplyr::summarise(count_gene = n(),mean_mutant_s=mean(e_s),mean_error_s=mean(u2_u1_s),mean_relative_s=mean(relative_s))
summary(count_fitness$mean_mutant_s);
summary(count_fitness$mean_error_s)

special_gene <- fitness_double[fitness_double$gene1=="C315T",]
special_gene$label <- paste(special_gene$gene1,special_gene$gene2,sep = "+")
pdf("C315T.pdf", width = 10, height = 8)
ggplot(special_gene, aes(x =s1_s2 , y = s12)) +
  geom_point(aes(
    shape = ifelse(s12 < s1_s2, 25,17),
    fill = ifelse(s12 < s1_s2,  "gray30", "gray60"),
    color = ifelse(s12 < s1_s2, "gray30", "gray60")
  ), size = 10, stroke = 1, show.legend = FALSE) +  
   geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_shape_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  theme_bw() +
  labs(x = "Fitness effect expected by additive model", y = "Observed fitness effect of\n double mutants", title = "") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    axis.title.x = element_text(size = 25,face = "bold",margin = margin(t = 20)),
    axis.title.y = element_text(size = 25, face = "bold",margin = margin(r = 20)),
    legend.position = "none",  
    plot.margin = unit(c(1, 3, 1, 1), "lines")  
  ) +
  scale_x_continuous(limits = c(-0.0085,0))+
  scale_y_continuous(limits = c(-0.0085,0))

dev.off()

# ----------------------------------------------------------------------------------------- #
#                                             figure4. C                                    #
# ----------------------------------------------------------------------------------------- #
fitness_list <- list.files(pattern = "\\.RData$")
genename <- c("AGS","AGY","AUS","AUU","AUY","TGS","TGY","TUS","TUU","TUY")
gene_fitness <- list()
for(i in seq(1,length(fitness_list)))
{
  tt <- load(fitness_list[[i]])
  gene_fitness[[i]] <- count_fitness[ sample(1:length(count_fitness[,1]),100),]
}
names(gene_fitness) <- genename
fitness <- matrix(0,nrow = 100,ncol = 10)
colnames(fitness) <- genename
fitness <- as.data.frame(fitness)
for(i in seq(1,10))
{
  tt <- gene_fitness[[i]]
  fitness[,i] <- tt[,4]-tt[,3]
}
fitness_long <- melt(fitness, variable.name = "Gene", value.name = "Fitness")

ggplot(fitness_long, aes(x = factor(Gene, levels = c("TGS","AGS","TUS","AUS","TGY","AGY","TUY","AUY","TUU","AUU")), y = Fitness)) +
  geom_boxplot() +  
  geom_jitter(width = 0.2, alpha = 0.5, color = "gray30") +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",linewidth = 0.8) +  
  theme_minimal() +
  labs(title = "", x = "", y = "")+
  theme(
    axis.title.x = element_text(size = 16),                
    axis.text.x = element_blank(),               
    axis.text.y = element_text(size = 22)                
  )
