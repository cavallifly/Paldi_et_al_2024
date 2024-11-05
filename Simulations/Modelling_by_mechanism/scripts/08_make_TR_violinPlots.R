options(scipen=9999)
library(dplyr)
library(ggplot2)
library(ggpubr)
library("ggrepel")
library(scales)
library("RColorBrewer")
pi=3.14159265359

data <- read.table("AvgTransRatio_perChrom.tab", header=T)
print(head(data))
#quit()
print(unique(data$sample))

# Effect of the number of copies
levels_Allexp <- c("DMSO","DMSO_A","DMSO_B",
	           "TSA" , "TSA_A" , "TSA_B" ,
		   "TSA24hREC", "TSA24hREC_A", "TSA24hREC_B")

levels_exp <- c("DMSO","DMSO_A","DMSO_B",
                "TSA" , "TSA_A" , "TSA_B")
#levels_exp <- c("TSA" , "TSA_A" , "TSA_B")

levels_sim <- unique(data$sample[-grep("TSA",data$sample)])
levels_sim <- unique(levels_sim[-grep("DMSO",levels_sim)])

levels <- c(levels_exp,levels_sim)
colors <- c('blue','white','white',
            'red','white','white',
# 	    'orange',
	    rep(c(rgb(1,0,0,0.50),"white","white"),as.integer(length(levels_sim)/3))
	    )

print(levels)
print(length(levels))
print(length(colors))

data$sample<- factor(data$sample, levels = levels)
data <- data[!is.na(data$sample),]
print(data$sample)

mC <- c("DMSO","DMSO_A","DMSO_B", "TSA", "TSA_A", "TSA_B", "TSA24hREC", "TSA24hREC_A", "TSA24hREC_B")
for(sample in mC)
{
    #print(sample)
    result <- c(sample,as.numeric(quantile(data[data$sample==sample,]$avg_trans_ratio, c(.25, .50, .75), na.rm = TRUE)))
    #    print(result)
    if (!exists("quartiles"))
    {
	print(paste0("Adding ",sample))
	quartiles <- data.frame()    
        quartiles <- as.data.frame(t(result))
    } else {
	print(paste0("Adding ",sample))
        quartiles <- rbind(quartiles,as.data.frame(t(result)))
    }
}
colnames(quartiles) <- c("sample","first","second","third")
print(quartiles)

for(sample in unique(data$sample))
{
    print(sample)
    result <- c(sample,as.numeric(quantile(data[data$sample==sample,]$avg_trans_ratio, c(.25, .50, .75), na.rm = TRUE)))
    print(result)
}
#quit()

# Plot the ashapeVol per replica for some values of timestep
outFile <- paste0("AvgTransRatio_perChrom.pdf")
if(!file.exists(outFile))
{
    print(paste0("outFile ",outFile))

    p <- ggplot(data, aes(x=factor(sample), y=avg_trans_ratio, fill=factor(sample))) +
             annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(quartiles[c(1),]$first), ymax = as.numeric(quartiles[c(1),]$third), fill = rgb(0,0,1,0.50), color=NA) + # DMSO
	     annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(quartiles[c(2),]$first), ymax = as.numeric(quartiles[c(2),]$third), fill = rgb(0,0,1,0.80), color=NA) + # DMSO_A
	     annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(quartiles[c(3),]$first), ymax = as.numeric(quartiles[c(3),]$third), fill = rgb(0,0,1,0.20), color=NA) + # DMSO_B
             #annotate(geom = "rect", ymin = as.numeric(quartiles[c(4),]$first), ymax = as.numeric(quartiles[c(4),]$third), xmin = -Inf, xmax = Inf, fill = rgb(1,0,0,0.50), color=NA) + # TSA
             #annotate(geom = "rect", ymin = as.numeric(quartiles[c(5),]$first), ymax = as.numeric(quartiles[c(5),]$third), xmin = -Inf, xmax = Inf, fill = rgb(1,0,0,0.20), color=NA) + # TSA_B
             #annotate(geom = "rect", ymin = as.numeric(quartiles[c(6),]$first), ymax = as.numeric(quartiles[c(6),]$third), xmin = -Inf, xmax = Inf, fill = rgb(1,0,0,0.80), color=NA) + # TSA_A
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     ylim(0.15,.75) +
	     #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
	     #labels = trans_format("log10", math_format(10^.x))) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 #axis.title = element_text(size=12),	   
		 axis.text = element_text(face="bold",size=12),
		 #axis.text.y = element_text(face="bold",size=18),
		 axis.text.x = element_text(angle=60, hjust=1),
		 #axis.ticks.x = element_blank(),
		 #axis.ticks.y = element_blank(),
           legend.position = "none") +
	   labs(x="",y="Avg. trans-ratio",fill="") +     
	   guides(fill="none")

    #pdf(outFile,width=28)
    pdf(outFile)    
    print(p)
    dev.off()

}

quit()
