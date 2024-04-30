suppressMessages(library(ggplot2))

setwd("/project2/xinhe/xsun/neuron_simulation/2.torus")

load("enrichment.rdata")

se <- (enrichment$high - enrichment$low) / (2*1.96)
z <- enrichment$estimate/se
p <- exp(-0.717*2 - 0.416*z^2)
enrichment$p <- p
enrichment$lp <- -log10(p)
#save(enrichment, file = "enrichment_p.rdata")

class_snp <- levels(as.factor(enrichment$snp))

for (i in 1:length(class_snp)){
  
  enrichment_tmp <- enrichment[enrichment$snp == class_snp[i],]
  
  rownames(enrichment_tmp) <- enrichment_tmp$term
  enrichment_tmp <- enrichment_tmp[c("schizophrenia","bipolar","autism","depression","adhd","neuroticism","intelligence","insomnia","alzheimer2","Alzheimer_hg19_hm","bmi","T2D","Parkinson's disease"),]
  enrichment_tmp$term <- c("Schizophrenia","Bipolar","Autism","Depression","ADHD","Neuroticism","Intelligence","Insomnia","Alzheimer","Alzheimer_newer","BMI","T2D","Parkinson's disease")
  
  enrichment_tmp$order <- seq(1,nrow(enrichment_tmp))
  enrichment_tmp <- enrichment_tmp[order(enrichment_tmp$order,decreasing = T),]
  
  colnames(enrichment_tmp)[1] <- "Traits"
  enrichment_tmp$Traits <- as.factor(enrichment_tmp$Traits)
  
  if (min(enrichment_tmp$low) < 0) {
    xlim_l <- floor(min(enrichment_tmp$low)) -1
  }else {
    xlim_l <- ceiling(min(enrichment_tmp$low)) + 1
  }  
  if (max(enrichment_tmp$high) < 0) {
    xlim_u <- floor(max(enrichment_tmp$high)) -1
  }else {
    xlim_u <- ceiling(max(enrichment_tmp$high)) +1
  }

  
  interval <- ggplot(enrichment_tmp, aes(y=seq(1,nrow(enrichment_tmp)), x=estimate)) + theme_bw(base_line_size =0.3) +
    geom_errorbar(aes(xmin=low, xmax=high), width=0.1, color = "black",alpha=0.6) +
    geom_point( color = "brown",aes(size = lp)) + 
    #aes(size = lp) +
    ylab("Traits") +  xlab(bquote(.(class_snp[i]) ~ "enrichment (log"[2] ~ ")")) + 
    xlim(xlim_l,xlim_u)  + 
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12, color = "black")) +
    geom_vline(aes(xintercept = 0),linetype="longdash") +
    scale_y_continuous(breaks = seq(1,nrow(enrichment_tmp)),labels = enrichment_tmp$Traits) +
    labs(size = expression(-log[10]("p-value"))) 
  
  file_save <- paste0(class_snp[i],".pdf")
  ggsave(filename = file_save, plot = interval, dpi = 600,width = 5.5, height = 3, limitsize = FALSE )
  
  
}