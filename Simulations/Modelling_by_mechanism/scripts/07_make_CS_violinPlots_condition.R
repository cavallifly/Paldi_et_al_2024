options(scipen=9999)
library(dplyr)
library(ggplot2)
library(ggpubr)
library("ggrepel")
library(scales)
library("RColorBrewer")

args = commandArgs(trailingOnly=TRUE)
pi=3.14159265359

data <- read.table("CS_perChrom.tab", header=T)
print(head(data))
print(unique(data$sample))

condition <- args[1]

levels_sim <- unique(data$sample[-grep("TSA",data$sample)])
levels_sim <- unique( levels_sim[-grep("DMSO",levels_sim)])
levels_exp <- c(paste0(condition,"_A"),paste0(condition,"_B"))

levels <- c(levels_exp,levels_sim)

data$sample<- factor(data$sample, levels = levels)

data <- data[!is.na(data$sample),]

print(length(levels))

colors <- c('red','blue',
	    rep(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),as.integer(length(levels_sim)))
	    )


# Compute stats on the dataset at a fixed timestep
#stats <- data %>%
#    group_by(sample) %>%
#    dplyr::reframe(n=n(), quantile_df(CS), mean=mean(CS), sd=sd(CS), median=median(CS), mad=mad(CS))
#stats <- unique(stats)
#print(head(stats))

#for(sample in unique(data$sample))

for(sample in levels_exp)
{
    print(sample)
    result <- c(sample,as.numeric(quantile(data[data$sample==sample,]$CS, c(.25, .50, .75), na.rm = TRUE)))
    print(result)
    if (!exists("quartiles"))
    {
	print(sample)
	quartiles <- data.frame()    
        quartiles <- as.data.frame(t(result))
    } else {
	print(paste0("Adding ",sample))
        quartiles <- rbind(quartiles,as.data.frame(t(result)))
    }
}
#for(sample in levels)
#{
#	#print(sample)
#    result <- as.data.frame(t(c(sample,as.numeric(quantile(data[data$sample==sample,]$CS, c(.25, .50, .75), na.rm = TRUE)))))
#    colnames(result) <- c("sample","first","second","third")
#    print(result)
#}
colnames(quartiles) <- c("sample","first","second","third")
print(quartiles)
print(as.numeric(quartiles[c(1),]$first))

# Plot the ashapeVol per replica for some values of timestep
outFile <- paste0("CS_perChrom.pdf")
if(!file.exists(outFile))
{
    print(paste0("outFile ",outFile))

    p <- ggplot(data, aes(x=factor(sample), y=CS, fill=factor(sample))) +
             annotate(geom = "rect", ymin = as.numeric(quartiles[c(1),]$first), ymax = as.numeric(quartiles[c(1),]$third), xmin = -Inf, xmax = Inf, fill = rgb(1,0,0,0.5), color=NA) +
             annotate(geom = "rect", ymin = as.numeric(quartiles[c(2),]$first), ymax = as.numeric(quartiles[c(2),]$third), xmin = -Inf, xmax = Inf, fill = rgb(0,0,1,0.5), color=NA) +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     ylim(0.80,1.75) +
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
	   labs(x="",y="Compartment strength",fill="") +     
	   guides(fill="none")

    #pdf(outFile,width=28)
    pdf(outFile, height=10)    
    print(p)
    dev.off()

}

quit()
