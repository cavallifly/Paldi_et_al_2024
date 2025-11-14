#cat <(awk '{print $0,"A"}' A_compartments_cis_sTSA.bed) <(awk '{print $0,"B"}' B_compartments_cis_sTSA.bed) | sort -k 1,1 -k 2,2n -k 3,3n > compartments_cis_sTSA.bed 
#cat <(awk '{print $0,"A"}' A_compartments_cis_sDMSO.bed) <(awk '{print $0,"B"}' B_compartments_cis_sDMSO.bed) | sort -k 1,1 -k 2,2n -k 3,3n > compartments_cis_sDMSO.bed 


options(scipen=9999)
library(dplyr)
library(ggplot2)
library(ggpubr)
library("ggrepel")
library(scales)
library("RColorBrewer")
pi=3.14159265359

allLengths <- data.frame()

#levels <- c(list.files(pattern="comp"),list.files(pattern="comp"))
#levels <- gsub(".bed","",levels)
#print(levels)
#quit()

for(file in c("compartments_cis_sTSA.bed","compartments_cis_sDMSO.bed"))
{
    #print(file)

    name <- gsub(".bed","",file)
    name <- gsub("compartments_cis_","",name)

    compartments = read.table(file)
    colnames(compartments) <- c("chrom","start","end","compartment")

    compartments$nextStart <- c(compartments$start[2:length(compartments$start)],compartments$end[length(compartments$end)])   
    compartments$gap   <- compartments$nextStart-compartments$end
    compartments[compartments$gap<0,]$gap <- 0

    # Identify the gaps between compartments
    compartmentsTmp <- data.frame()
    for(r in 1:nrow(compartments))
    {
        row <- compartments[r,]
	print(row)
	line <- data.frame(row$chrom,row$start,row$end,row$compartment)
        names(line) <- c("chrom","start","end","compartment")
	if(is.null(compartmentsTmp))
	{
	    compartmentsTmp <- line
	} else {
 	    compartmentsTmp <- rbind(compartmentsTmp,line)
	}

	if(row$gap > 0)
	{
	    newline <- data.frame(row$chrom,row$end,row$nextStart,"GAP")
	    names(newline) <- c("chrom","start","end","compartment")
	    compartmentsTmp <- rbind(compartmentsTmp,newline)
	}

    }
    colnames(compartmentsTmp) <- c("chrom","start","end","compartment")
    compartments <- compartmentsTmp
    print(compartments)

    compartments$length <- compartments$end-compartments$start

    tmpDf <- data.frame(paste0(compartments$compartment,"_",name),compartments$length/1e3,compartments$compartment)
    if(is.null(allLengths))
    {
        allLengths <- tmpDf       
    } else {
        allLengths <- rbind(allLengths,tmpDf)
    }

    print(nrow(allLengths))
}
colnames(allLengths) <- c("sample","compartmentLength","compartment")
#print(head(allLengths))

data <- allLengths
print(head(data))

# Computing the statistics
stats <- data %>%
    group_by(sample) %>%
    dplyr::reframe(n=n(), sample=unique(sample), mean=mean(compartmentLength), sd=sd(compartmentLength), median=median(compartmentLength), mad=mad(compartmentLength))
print(head(stats))

# Making the plot now
outFile <- paste0("compartmentLength.pdf")

#print(unique(data$sample))
#quit()
#levels = unique(data$sample)
levels <- c("A_sDMSO", "B_sDMSO", "GAP_sDMSO",
            "A_sTSA" , "B_sTSA",  "GAP_sTSA") 
print(levels)
print(paste0("Number of replicates ",length(levels)))

data$sample<- factor(data$sample, levels = levels)

colors <- c(rgb(0  , 0, 255, max = 255, alpha = 225, names = "blue"),
       	    rgb(0  , 0, 255, max = 255, alpha = 125, names = "blue50"),
	    "gray",
       	    rgb(255, 0,   0, max = 255, alpha = 125, names = "red"),
       	    rgb(255, 0,   0, max = 255, alpha = 125, names = "red50"),
	    "gray"
	    )	    
names(colors) <- as.character(levels)
print(colors)

if(!file.exists(outFile))
{
    print(paste0("outFile ",outFile))

    p <- ggplot(data, mapping=aes(x=sample, y=compartmentLength, fill=sample)) +
     	     geom_violin(position = position_dodge(1), trim=T, scale = "width") +
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +
	     #geom_jitter(width=0.1, size = 0.1) +
	     scale_color_manual(values=colors) +
	     scale_fill_manual(values=colors) +
	     #ylim(0.0,.80) +
	     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
	     labels = trans_format("log10", math_format(10^.x))) +
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
	   labs(x="",y="Compartment length (kb)",fill="") +
   	   geom_point(data=stats,mapping=aes(x=sample,y=mean), color="magenta", group=1) +
   	   geom_text_repel(data=stats,mapping=aes(x=sample,label=as.integer(mean),  y=mean  ), color="magenta", group=1) +
   	   geom_text_repel(data=stats,mapping=aes(x=sample,label=as.integer(median),y=median), color="black", group=1) +	   
	   geom_text(data=stats,mapping=aes(x=as.factor(sample),y=5,label=n), color="black", group=1) +
	   guides(fill="none")

    pdf(outFile)
    print(p)
    dev.off()

}
