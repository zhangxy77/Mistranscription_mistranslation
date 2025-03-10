library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
# ----------------------------------------------------------------------------------------- #
#                                             figure5                                       #
# ----------------------------------------------------------------------------------------- #
num_generations <- 10000;
num_genes <- 500; ## only consider the top 100 genes. Other lower expressed genes are of no consequence in terms of translation error
min_expr <- 100000;
powerlaw_param <- 2;
sum(rpar(num_genes,min_expr,powerlaw_param)); ## this should be at around the 1e8 level (total number of proteins in Ecoli at about 1e6, yeast at about 1e8, mammalian at about 1e10, according to Bionumbers)
pop_size <- 1000;
Lnt <- 1000; ## nucleotide length
Lp <- 300; ## protein length
e <- 100;##epistasis
shape_est_TR <- 5.86
scale_est_TR <- 0.173
shape_est_TL <- 5.44
scale_est_TL <- 0.197

ancPop <- data.frame(
  ind = rep(1:pop_size,each=num_genes),
  beta1=rep(runif(num_genes, 10^-5, 10^-5),pop_size),
  beta2=rep(runif(num_genes, 10^-4, 10^-4),pop_size),
  N = rep(rpar(num_genes,min_expr,powerlaw_param),pop_size) ## using power-law distribution for expression
);
pop <- ancPop;

dfParamSearch <- expand.grid(
  mut_rate = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01),
  param_c = c(1e-5,1e-6,1e-7,1e-8)
);



dfParamSearch <- dfParamSearch %>%
  mutate(filename = paste0("01.sim_",mut_rate,"_",param_c,".RData"));
allPop <- list();
allTrackFitness <- list();
allTrackCor <- list();
for(thisParam in c(1:nrow(dfParamSearch))) {
  cat(paste0("\rLoading ",thisParam,"/",nrow(dfParamSearch)));
  load(dfParamSearch[thisParam,"filename"]);
  allPop[[filename]] <- pop;
  allTrackFitness[[filename]] <- trackFitness;
  allTrackCor[[filename]] <- trackCor;}

## summarize them
allCor <- data.frame();
tt_test <- matrix(0,nrow = 2,ncol = 2)
colnames(tt_test) <- c("β2 < medianvalue","β2 > medianvalue")
rownames(tt_test) <- c("β1 < medianvalue","β1 > medianvalue")
tt_test <- as.data.frame(tt_test)
for(thisParam in c(1:nrow(dfParamSearch))) {
  myF <- dfParamSearch[thisParam,"filename"];
  pearson_list <- list();pearson_plist <- list();
  spearman_list <- list();spearman_plist <- list();
  for(j in seq(1,num_genes*1000,num_genes))
  {
    k <- j+num_genes-1
    pearson.obj <- cor.test(allPop[[myF]]$beta1[j:k],allPop[[myF]]$beta2[j:k],method="pearson");
    spearman.obj <- cor.test(allPop[[myF]]$beta1[j:k],allPop[[myF]]$beta2[j:k],method="spearman");
    pearson_list <- append(pearson_list,pearson.obj$estimate)
    pearson_plist <- append(pearson_plist,pearson.obj$p.value)
    spearman_list <- append(spearman_list,spearman.obj$estimate)
    spearman_plist <- append(spearman_plist,spearman.obj$p.value)
  }
  usingCor <- rev(unlist(pearson_list));#取100个cell的cor
  tt_test[1,1] <- table(allPop[[myF]]$beta1[1:num_genes]<median(allPop[[myF]]$beta1[1:num_genes])&allPop[[myF]]$beta2[1:num_genes]<median(allPop[[myF]]$beta2[1:num_genes]))[[2]];
  tt_test[1,2] <- table(allPop[[myF]]$beta1[1:num_genes]<median(allPop[[myF]]$beta1[1:num_genes])&allPop[[myF]]$beta2[1:num_genes]>median(allPop[[myF]]$beta2[1:num_genes]))[[2]];
  tt_test[2,1] <- table(allPop[[myF]]$beta1[1:num_genes]>median(allPop[[myF]]$beta1[1:num_genes])&allPop[[myF]]$beta2[1:num_genes]<median(allPop[[myF]]$beta2[1:num_genes]))[[2]];
  tt_test[2,2] <- table(allPop[[myF]]$beta1[1:num_genes]>median(allPop[[myF]]$beta1[1:num_genes])&allPop[[myF]]$beta2[1:num_genes]>median(allPop[[myF]]$beta2[1:num_genes]))[[2]];
  tt.obj <- chisq.test(tt_test,correct = T)
  OR_result <- oddsratio(as.matrix(tt_test))
  allCor <- allCor %>% 
    rbind.fill(
      data.frame(dfParamSearch[thisParam,]) %>%
        mutate(r = mean(unlist(pearson_list)),
               r.p = mean(unlist(pearson_plist)),
               r.p.adjust = mean(p.adjust(unlist(pearson_plist), method = "bonferroni")),
               rho = mean(unlist(spearman_list)),
               rho.p =mean(unlist(spearman_plist)),
               rho.p.adjust = mean(p.adjust(unlist(spearman_plist), method = "bonferroni")),
               fracCorNeg = mean(usingCor<0,na.rm=T), ## sample from the last 3000 generations
               fracSigCorNeg = mean(rev(unlist(pearson_list))< -0.088,na.rm=T),
               chisq_p = -log(tt.obj$p.value),
               OR = OR_result$measure[[2]],
               OR_p = OR_result$p.value[[4]],
               log_OR_p <- -log(OR_result$p.value[[4]]) 
        )
    )
}

allCor$log_OR_p <- -log(allCor$OR_p)
allCor <- allCor %>% mutate(r_p_text = case_when(
  r.p < 0.05 ~ "*"))
allCor <- allCor %>% mutate(OR_p_text = case_when(
  OR_p < 0.05 ~ "*"))
allCor <- allCor %>% mutate(Frac_p_text = case_when(
  fracSigCorNeg > 0.5 ~ "*"))

save(allCor,allTrackCor,allTrackFitness,dfParamSearch,file="epistasis_allCor.RData")

plotname <- "mix_heatmap.pdf"
pdf(plotname)  
allCor %>%
  ggplot(aes(fill=fracCorNeg, x=as.factor(mut_rate), y=as.factor(param_c))) +
  geom_tile() +
  geom_text(aes(label=Frac_p_text), col ="black", size = 7) +
  labs(title = "", x = "Per-generation probability of error rate changes", y = "Toxicity coefficient of erroneous molecules", 
       fill = "Fraction of samples with negative \nmistranscription-mistranslation correlation") +
  scale_fill_distiller(type="div", palette = "RdYlBu", limits = c(0, 1)) +
  theme(
    axis.text = element_text(size = 13, face = "bold"), 
    axis.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 18, hjust = 0.5, vjust = 1.5, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title.align = 0.5, 
    legend.key.width = unit(2, "cm"), 
    legend.key.height = unit(0.5, "cm"),
    legend.margin = margin(t = 0.05, b = 0),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 20)
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.8, "npc")))

allCor %>%
  ggplot(aes(fill=r, x=as.factor(mut_rate),y=as.factor(param_c))) +
  geom_tile() +
  geom_text(aes(label=r_p_text),col ="black",size = 7)+
  labs(title = "", x = "Per-generation probability of error rate changes", y = "Toxicity coefficient of erroneous molecules", 
       fill = "Averaged mistranscription-mistranslation correlation") +
  scale_fill_distiller(type="div",palette = "RdYlBu", limits=c(-0.13, 0.13))+
  theme(axis.text = element_text(size = 13, face = "bold"), 
        axis.title = element_text(size = 21, face = "bold"),
        legend.text = element_text(size = 13,face = "bold"),
        legend.title  = element_text(size = 18,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,vjust = 1.5, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal", 
        legend.justification = "center",
        legend.title.align = 0.5,
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.5, "cm"), 
        legend.margin = margin(t = 5, b = 0), 
        plot.margin = margin(t = 20, r = 20, b = 40, l = 20) 
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.8, "npc")))  
dev.off();

## pick one specific parameter set and plot
thisParam <- 26
load(dfParamSearch[thisParam,"filename"])
myF <- dfParamSearch[thisParam,"filename"];
filename <- myF
plotname <- paste(filename,".pdf",sep = "")
cor.test(pop$beta1[1:num_genes],pop$beta2[1:num_genes],method="spearman");
cor.test(pop$beta1[1:num_genes],pop$beta2[1:num_genes],method="pearson");
pop %>% 
  mutate(geneId = rep(1:num_genes,pop_size) ) %>%
  group_by(geneId) %>%
  dplyr::summarise(beta1 = mean(beta1),beta2=mean(beta2),N=mean(N)) %>%
  filter(beta1 > 1e-8, beta2 > 1e-8) %>%
  cor.test(~ beta1 + beta2, data=., method="spearman")
cor.test(pop$N[1:num_genes],pop$beta1[1:num_genes],method="spearman");
cor.test(pop$N[1:num_genes],pop$beta2[1:num_genes],method="spearman");
exp_pop <- pop[1:num_genes,]
exp_pop <- exp_pop[order(exp_pop$N,decreasing = T),]
exp_level <- exp_pop$N[length(exp_pop$N)*0.25]
exp_pop$exp_type <- ifelse(exp_pop$N>= exp_level ,"high","low")
exp_pop$exp_type <- ifelse(exp_pop$N>= mean(exp_pop$N) ,"high","low")
exp_pop_high <- exp_pop[exp_pop$exp_type == "high",]
exp_pop_low <- exp_pop[exp_pop$exp_type == "low",]
pop$N_CPM <- (pop$N/sum(pop$N))*1e6

pdf(plotname);
pop %>% 
  mutate(geneId = rep(1:num_genes,pop_size) ) %>%
  group_by(geneId) %>%
  dplyr::summarise(beta1 = mean(beta1),beta2=mean(beta2),N=mean(N)) %>%
  ggplot(aes(x=beta1,y=beta2,color=log(N))) +
  geom_point() +
  labs(title = myF)+
  scale_color_distiller(type="div",palette=1) +
  scale_x_log10() + scale_y_log10()+
  theme(axis.text= element_text(size = 23, face = "bold"),
        axis.title= element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 25, face = "bold"))

pop$log_N <- log(pop$N)
plot(pop$beta1[1:num_genes],pop$beta2[1:num_genes], cex.lab = 2, cex.axis = 1.5, cex.main = 2,font.lab = 2,font.axis = 2,font.main = 2);
plot(log(pop$N_CPM[1:num_genes], base = 10),pop$beta1[1:num_genes],xlab = "Expression", ylab = "Mistranscription rate", cex.lab = 2, cex.axis = 1.5, cex.main = 2,font.lab = 2,font.axis = 2,font.main = 2);
plot(log(pop$N_CPM[1:num_genes], base = 10),pop$beta2[1:num_genes], xlab = "Expression", ylab = "Mistranslation rate",cex.lab = 2, cex.axis = 1.5, cex.main = 2,font.lab = 2,font.axis = 2,font.main = 2);
# expression
ggplot(pop[1:num_genes, ], aes(x = N_CPM, y = beta1)) +
  geom_point(shape = 1,size=2)  +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(
    labels = c(expression(italic(0)),   
               expression(0.5),   
               expression(1.0),   
               expression(1.5),
               expression(2.0)),   
    expand = c(0, 0),   
    breaks = c(0, 0.000005, 0.000010, 0.000015, 0.000020),   
    limits = c(-0.000001, 0.000025)
  ) +
  labs(x = "Expression (CPM)", y = "Mistranscription rate") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 23, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank()  
  )
ggplot(pop[1:num_genes, ], aes(x = N_CPM, y = beta2)) +
  geom_point(shape = 1,size=2) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(labels = c(expression(italic(0)),   
                                expression(2.0),   
                                expression(4.0),   
                                expression(6.0)),   
                     expand = c(0, 0),   
                     breaks = c(0, 0.0002, 0.0004, 0.0006),   
                     limits = c(-0.00001, 0.00065)) +
  labs(x = "Expression (CPM)", y = "Mistranslation rate") +
  theme_bw() + 
  theme(
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 23,face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank()
  )

plot(1:length(trackFitness),unlist(trackFitness),type="l", xlab = "Generation", ylab = "Fitness",cex.lab = 2, cex.axis = 1.5, cex.main = 2,font.lab = 2,font.axis = 2,font.main = 2);
plot(1:length(trackCor),unlist(trackCor),type="l", xlab = "Generation", ylab = "Correlation",cex.lab = 2, cex.axis = 1.5, cex.main = 2,font.lab = 2,font.axis = 2,font.main = 2);
#expression
percent_thresholds <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1.0)
result_total <- data.frame(Percent = percent_thresholds, RealCorrelation = NA, RandomCorrelation = NA,sd_real= NA,sd_random=NA)
results <- data.frame(Percent = percent_thresholds, RealCorrelation = NA, RandomCorrelation = NA)
real_random_1 <- data.frame(real=rep(0,1000),random=rep(0,1000))
real_random_5 <- data.frame(real=rep(0,1000),random=rep(0,1000))
real_random_10 <- data.frame(real=rep(0,1000),random=rep(0,1000))
real_random_20 <- data.frame(real=rep(0,1000),random=rep(0,1000))
real_random_50 <- data.frame(real=rep(0,1000),random=rep(0,1000))
real_random_100 <- data.frame(real=rep(0,1000),random=rep(0,1000))
t <- 1
for(i in seq(1,length(pop[,1]),500))
{
  end <- i+num_genes-1 
  exp_pop <- pop[i:end,]
  exp_pop <- exp_pop[order(exp_pop$N,decreasing = T),]
 
  for (threshold in percent_thresholds) {
    exp_level <- quantile(exp_pop$N, 1 - threshold)
    
    high_exp_genes <- exp_pop[exp_pop$N >= exp_level, ]
    
    if (nrow(high_exp_genes) > 1) {
      results$RealCorrelation[results$Percent == threshold] <- cor(high_exp_genes$beta1, high_exp_genes$beta2)
    } else {
      results$RealCorrelation[results$Percent == threshold] <- NA
    }
    
    random_genes <- exp_pop[sample(nrow(exp_pop), size = ceiling(threshold * nrow(exp_pop))), ]
    
    if (nrow(random_genes) > 1) {
      results$RandomCorrelation[results$Percent == threshold] <- cor(random_genes$beta1, random_genes$beta2)
    } else {
      results$RandomCorrelation[results$Percent == threshold] <- NA
    }
  }
  real_random_1[t,1] <- results[1,2];real_random_1[t,2] <- results[1,3]
  real_random_5[t,1] <- results[2,2];real_random_5[t,2] <- results[2,3]
  real_random_10[t,1] <- results[3,2];real_random_10[t,2] <- results[3,3]
  real_random_20[t,1] <- results[4,2];real_random_20[t,2] <- results[4,3]
  real_random_50[t,1] <- results[5,2];real_random_50[t,2] <- results[5,3]
  real_random_100[t,1] <- results[6,2];real_random_100[t,2] <- results[6,3]
  t <- t+1
}
result_total[1,2] <- mean(real_random_1[,1]);result_total[1,3] <- mean(real_random_1[,2]);result_total[1,4] <- sd(real_random_1[,1]);result_total[1,5] <- sd(real_random_1[,2])
result_total[2,2] <- mean(real_random_5[,1]);result_total[2,3] <- mean(real_random_5[,2]);result_total[2,4] <- sd(real_random_5[,1]);result_total[2,5] <- sd(real_random_5[,2])
result_total[3,2] <- mean(real_random_10[,1]);result_total[3,3] <- mean(real_random_10[,2]);result_total[3,4] <- sd(real_random_10[,1]);result_total[3,5] <- sd(real_random_10[,2])
result_total[4,2] <- mean(real_random_20[,1]);result_total[4,3] <- mean(real_random_20[,2]);result_total[4,4] <- sd(real_random_20[,1]);result_total[4,5] <- sd(real_random_20[,2])
result_total[5,2] <- mean(real_random_50[,1]);result_total[5,3] <- mean(real_random_50[,2]);result_total[5,4] <- sd(real_random_50[,1]);result_total[5,5] <- sd(real_random_50[,2])
result_total[6,2] <- mean(real_random_100[,1]);result_total[6,3] <- mean(real_random_100[,2]);result_total[6,4] <- sd(real_random_100[,1]);result_total[6,5] <- sd(real_random_100[,2])

result_long <- melt(result_total, id.vars = "Percent", measure.vars = c("RealCorrelation", "RandomCorrelation"))



ggplot(result_long, aes(x = factor(Percent, levels = unique(result_total$Percent)), 
                        y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = value - ifelse(variable == "RealCorrelation", 
                                          result_total$sd_real,
                                          result_total$sd_random),
                    ymax = value + ifelse(variable == "RealCorrelation", 
                                          result_total$sd_real,
                                          result_total$sd_random)),
                position = position_dodge(0.7), width = 0.2) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = c("gray30", "gray60"), 
                    labels = c("real", "random"),
                    name = "") +
  scale_y_continuous(limits = c(-0.7, 0.6), 
                     breaks = seq(-0.7, 0.6, by = 0.1),
                     labels = label_number(accuracy = 0.1)) + 
  scale_x_discrete(labels = c("0.01" = "1%", "0.05" = "5%", "0.1" = "10%", "0.2" = "20%", "0.5" = "50%", "1" = "100%")) +
  labs(x = "High expression Genes", y = "Correlation", title = "") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 20, face = "bold"), 
        axis.title = element_text(size = 25, face = "bold"), 
        axis.text = element_text(size = 23, face = "bold"),
        plot.title = element_text(size = 25, face = "bold"))

dev.off();

beta_plot<- pop %>% 
  mutate(geneId = rep(1:num_genes,pop_size) ) %>%
  group_by(geneId) %>%
  dplyr::summarise(beta1 = mean(beta1),beta2=mean(beta2),N=mean(N)) 
beta_p1 <- ggplot(beta_plot,aes(y = reorder(geneId, beta1),x=beta1))+ 
  geom_point(size = 1)+  # Increase the size of the points 
  scale_x_continuous(labels = c(expression(italic(0)), 
                                expression(0.5), 
                                expression(1.0), 
                                expression(1.5),
                                expression(2.0),
                                expression(2.5)),
                     expand = c(0,0), 
                     breaks = c(0,0.000005,0.00001,0.000015,0.00002,0.000025), 
                     limits = c(0,0.000025)) +
  labs(x = "", y = "", title = "")+
  theme(axis.title.x = element_text(size=30, margin = margin(r = 20)), 
        axis.text.x = element_text(size = 23,face = "bold"), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=30, margin = margin(r = 20)), # Increase the distance between y-axis title and labels
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_line(color = "grey93"))+
  coord_cartesian(clip = "off")mean

beta_p2 <- ggplot(beta_plot,aes(y = reorder(geneId, beta2),x=beta2))+ 
  geom_point(size = 1)+  # Increase the size of the points 
  scale_x_continuous(labels = c(expression(italic(0)), 
                                expression(2), 
                                expression(4), 
                                expression(6),
                                expression(8)), 
                     expand = c(0,0), 
                     breaks = c(0,0.0002,0.0004,0.0006,0.0008), 
                     limits = c(0,0.0008))+
  labs(x = "", y = "", title = "")+
  theme(axis.title.x = element_text(size=30, margin = margin(r = 20)), 
        axis.text.x = element_text(size = 23, face = "bold"), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=30, margin = margin(r = 20)), # Increase the distance between y-axis title and labels
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_line(color = "grey93"))+
  coord_cartesian(clip = "off")

pdf("beta1.pdf")
ggMarginal(beta_p1, type="histogram", fill = "grey48", xparams = list(bins=40), yparams = list(bins=40))
dev.off()
pdf("beta2.pdf")
ggMarginal(beta_p2, type="histogram", fill = "grey48", xparams = list(bins=40), yparams = list(bins=40))
dev.off()
