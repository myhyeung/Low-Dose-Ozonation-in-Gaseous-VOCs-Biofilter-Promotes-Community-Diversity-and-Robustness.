
#alpha diversity plotting, uses package 'Phyloseq" for analysis and plotting, OTU and taxonomy was generated in QIIME2 and imported for downstream analysis
plot_richness(data,x="group",measures=c("Observed", "Shannon","Simpson")) + geom_boxplot()  + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1), legend.text = element_text(size=16),panel.border = element_rect(fill=NA,size = 1),panel.background = element_blank(),axis.text = element_text(size=14))


#beta analysis
or<-ordinate(data, formula = ~ Complex.carbs. + Carboxylic + Aminos + Amines + Polymers,method = "CCA")
or<-ordinate(data, formula = ~ Actinobacteria + Proteobacteria + Bacteroidetes + TM7 + Firmicutes,method = "CCA")
p0<-plot_ordination(data,or,shape = "group") + 
  theme(panel.border = element_rect(size=1,fill = NA), legend.text = element_text(size=16),panel.background = element_blank(),axis.text = element_blank(),axis.title = element_text(size = 20),legend.title = element_blank()) + 
  stat_ellipse(type = "norm",level = 0.6) + scale_shape_manual(values=c(4, 19))
p0 = p0 + geom_point(size=4)

arrowmat = vegan::scores(or, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map = aes(xend = 1*CCA1, yend = 1*CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 0.8 * CCA1, y = 0.8* CCA2, shape = NULL, color = NULL, 
                label = labels)
arrowhead = arrow(length = unit(0.5,"lines"))
p0 + geom_segment(arrow_map, size = 1, data = arrowdf, color = "red", 
                  arrow = arrowhead) 
geom_text(label_map, size = 4, data = arrowdf)

#Aldex filtering of PICRUSt2 pathways
path<-read.delim('picrust2.tsv',sep = '\t',header = T,row.names = 1)
path<-path[-1:-24]
path<-round(path)
groupn<-c(rep('0v',2),rep('100v',3))
groupn<-c(rep('ozone',12),rep('control',12))
pathpmaclean<-round(pathpma[rowSums(pathpma)>1500,])
x.all <- aldex(pathpmaclean, groupn, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)
ggplot2::ggplot(a, ggplot2::aes(x=effect,y=reorder(X1,effect))) + ggplot2::geom_bar(fill="white",color='black',position = 'dodge',stat="identity") 
ggplot2::theme(panel.background = element_blank())