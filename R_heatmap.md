Orhtologous genes

```R
install.packages("ComplexHeatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs-heatmap.pdf",width=7,height=7)
orthologs_heatmap <- heatmap.2(orthologs_matrix, margins = c(6, 12), dendogram="column", Colv=FALSE, cexCol=0.8, cexRow=0.8 , trace="none", tracecol="Gray", col=scale)
dev.off()
```


#Orthologous light genes
```R
install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

library("gplots")
library("Cairo")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_light_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "snow4", "gray36"), space = "rgb")(100)
#pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs_light_heatmap2.pdf",width=6,height=7)
#svg(file="Users/lopeze/Desktop/Bioinformatics/orthologs_light_heatmap2.svg", width=10, height=8)
image = heatmap.2(orthologs_matrix, margins = c(5, 17), dendrogram="none", Colv=FALSE, Rowv=FALSE, cexCol=0.6, cexRow=0.6 , trace="none", col=scale, cellnote=orthologs_matrix, notecol="white", notecex=0.4,
as.matrix(order_by27_T[rowsToDraw,]), colsep=c(5,8,9), rowsep=c(3, 14), sepwidth=c(0.1,0.1), sepcolor=c("white"))
plot(image)
dev.off()
ggsave("/Users/lopeze/Desktop/Bioinformatics/orthologs_light_heatmap2.svg",image, width=6)
```

#Orthologous light genes using superheat
```R
library("superheat")
library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_light_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs_light_heatmap.pdf",width=7,height=7)
superheat(X = orthologs_matrix, # heatmap matrix
          # change the size of the labels
          pretty.order.rows = FALSE,
          pretty.order.cols = FALSE,
          left.label.size = 0.4,
          bottom.label.size = 0.1,
          left.label.text.size = 3,
          bottom.label.text.size = 2,
          left.label.col = "white",
          bottom.label.col = "white",
          bottom.label.text.angle = 90,
          left.label.text.alignment = "left",
          grid.hline.size = 0.1,
          grid.vline.size = 0.1,
          scale = FALSE)
dev.off()
```

Alternative script: /Users/lopeze/Desktop/Bioinformatics

library("gplots") orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_genes.txt", col.names=1, row.names=1, check.names=FALSE) orthologs_matrix<-data.matrix(orthologs) scale <- colorRampPalette(c("white", "forestgreen", "green"), space = "rgb")(100) pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs-heatmap.pdf",width=7,height=7) orthologs_heatmap <- heatmap.2(orthologs_matrix, geom_tile(aes(fill = scale)), margins = c(5, 12), Colv=FALSE, dendrogram="row", cexCol=0.8, cexRow=0.8, trace="none", tracecol="Gray", col=scale) dev.off()

Alternative script using Pheatmap:

library("pheatmap")
library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_light_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs_light_heatmap2.pdf",width=7,height=7)
#mylwid = c(1.5,4,0.5) mylhei = c(1.5,4,1) mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
orthologs_heatmap <- pheatmap(orthologs_matrix, lmat=mylmat, lwid=mylwid, lhei=mylhei, margins = c(6, 12), dendogram="column", Colv=FALSE, cexCol=0.8, cexRow=0.8 , trace="none", tracecol="Gray", col=scale) dev.off()




==========
Promoter motif heatmap
==========

library("pheatmap")
library("gplots")
promoter<-read.table("/Users/lopeze/Desktop/Bioinformatics/Promoter_motif/Promotermotif_dunlap.txt", col.names=,1, row.names=1, check.names=FALSE)
promoter_matrix<-data.matrix(promoter)
scale <- colorRampPalette(c("white", "lightgreen", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/promoter-heatmap2.pdf",width=7,height=7)
promoter_heatmap <- heatmap.2(promoter_matrix,
cellnote=as.matrix(promoter), #put the values in the cell
notecol="black", #color of the values in the cell
margins = c(5, 12),
dendogram="none",
Rowv=FALSE, #implies dendrogram is not computed and reordered based on row means
Colv=FALSE, #implies dendrogram is not computed and reordered based on column means
cexCol=0.8, #specify the size of the column label
cexRow=0.8,
lhei=c(1,8),
lwid=c(2,2),#modifiy the size of the cells
trace="none", tracecol="Gray", col=scale)
dev.off()


library("pheatmap")
library("gplots")
promoter<-read.table("/Users/lopeze/Desktop/Bioinformatics/Promoter_motif/Promotermotif.txt", col.names=,1, row.names=1, check.names=FALSE)
promoter_matrix<-data.matrix(promoter)
scale <- colorRampPalette(c("white", "yellow", "orange", "red"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/Promoter_motif/promoter-heatmap_DH.pdf",width=6,height=7)
promoter-heatmap <- heatmap.2(promoter_matrix,
cellnote=as.matrix(promoter_matrix), #put the values in the cell
notecol="black", #color of the values in the cell
margins = c(5, 12),
dendogram="none",
Rowv=FALSE, #implies dendrogram is not computed and reordered based on row means
Colv=FALSE, #implies dendrogram is not computed and reordered based on column means
cexCol=0.8, #specify the size of the column label
cexRow=0.8,
lhei=c(1,8),
lwid=c(2,2),#modifiy the size of the cells
trace="none", tracecol="Gray", col=scale)
dev.off()

--------------
STATISTICS
--------------

#E1-E2 experiment

ANOVA 2 Factors

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/E1-E2/E1_data.csv")
attach(data)
boxplot(Diameter~Strain*Conditions)
hist(Diameter)
shapiro.test(data$Diameter) #p-value has to be over 0.5

anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model))

library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(aov.model, pairwise~Conditions|Strain, adjust="tukey")
groups<-cld(means, alpha= .05)
```

#D8-D10 V.spp experiment

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/D8-D10-D14V.spp/D8-D10_data.csv")
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/D8-D10-D14V.spp/all_data.csv")
attach(data)

#RMEL for unbalanced designs
library(lme4)
library(lmertest)
l1 <- lmer(Diameter ~ Experiment + Strain * Conditions + (1|Experiment:Strain:Conditions))
summary(l1)
anova(l1)


library(lsmeans)
means<-(lsmeans(l1, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#ANOVA for balanced sesign
anv.model<-aov(Diameter~Experiment+Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

#Test if data is Normally distributed
res=resid(anv.model)
hist(res)
qqnorm(res)
shapiro.test(res) #p-value has to be over 0.5

#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)


library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#ggplot
library("ggplot2")
library(Cairo)
svg("/Users/lopeze/Desktop/Statistics_R/Rings/D8-D10-D14V.spp/Vdspp_all.svg", width=10, height=8)
image=ggplot(data, aes(x=Conditions, fill=Conditions, y=Diameter)) +
stat_boxplot(geom="errorbar", size=0.5) +
     geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
     labs(y=("Colony diameter(mm)")) +
     facet_grid(~Strain, scale="free") +
     scale_fill_manual(values=c(DD="#999999", LL="#3366CC", Lr="#99CCFF", RD="#660066"))+
     #scale_x_discrete(limits=c("DD","LL","LD","RD"))+
     ylim(30,60)+
     ylab("Colony diameter (mm)")+
     coord_fixed(ratio = 0.2)+
     theme_bw() +
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=10),
           axis.text.x = element_text(colour="black", size=10 ,angle=45, hjust=1),
           axis.text.y = element_text(colour="black", size=10),
           axis.title.x = element_text(size=14),
           axis.title.y = element_text(size=14),
           strip.text = element_text(size=9),
           strip.background = element_rect(colour = "white"))
   plot(image)
   dev.off()
```

#Light differences in WT strains

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Light/WT/all_data.csv")
attach(data)

#RMEL for unbalanced designs
library(lme4)
l1 <- lmer(Diameter ~ Experiment + Strain * Conditions + (1|Experiment:Strain:Conditions))
summary(l1)
anova(l1)


library(lsmeans)
means<-(lsmeans(l1, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#ANOVA
anv.model<-aov(Diameter~Experiment+Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

#Test if data is Normally distributed
res=resid(anv.model)
hist(res)
qqnorm(res)
shapiro.test(res) #p-value has to be over 0.5

#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)


library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#ggplot
library("ggplot2")
library(Cairo)
svg("/Users/lopeze/Desktop/Statistics_R/Rings/D8-D10-D14V.spp/Vdspp_all.svg", width=10, height=8)
ggplot(data, aes(x=Conditions, fill=Conditions, y=Diameter)) +
stat_boxplot(geom="errorbar", size=0.5) +
     geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
     labs(y=("Colony diameter(mm)")) +
     facet_grid(~Strain, scale="free") +
     scale_fill_manual(values=c(DD="#999999", LL="#3366CC", LD="#99CCFF", RD="#660066"))+
     #scale_x_discrete(limits=c("DD","LL","LD","RD"))+
     ylim(30,60)+
     ylab("Colony diameter (mm)")+
     coord_fixed(ratio = 0.2)+
     theme_bw() +
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=10),
           axis.text.x = element_text(colour="black", size=10 ,angle=45, hjust=1),
           axis.text.y = element_text(colour="black", size=10),
           axis.title.x = element_text(size=14),
           axis.title.y = element_text(size=14),
           strip.text = element_text(size=9),
           strip.background = element_rect(colour = "white"))
   plot(image)
   dev.off()
```



#D3 Wc1

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/AFrq/D16_data/D16_53_data.csv")
attach(data)
boxplot(Diameter~Strain*Conditions, las = 2, cex.axis=0.6)


hist(Diameter)
shapiro.test(data$Diameter) #p-value has to be over 0.5

#If the data is not normally distributed:
datalog<-log(data$Diameter)
shapiro.test(datalog$Diameter)
hist(data.log)

#ANOVA
anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H
```

#D15 D16 and D17 medium

```R   
#data<-read.csv("/Users/lopeze/Desktop/Statistics_R/D15-D16_medium/D15_data/D15_med_data.csv")
#data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/D15-D16_medium/D16_data/D16_med_data.csv")
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/D15-D16_medium/D16_D17.csv")
attach(data)
boxplot(Diameter~Strain*Conditions, las = 2, cex.axis=0.6)


#ANOVA
options(max.print=100000)

anv.model<-aov(Diameter~Experiment+Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#Test for normality
res=resid(aov)
hist(res)
qqnorm(res)
shapiro.test(res) #p-value has to be over 0.5
#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)

#Boxplot

library("ggplot2")
library(Cairo)
svg("/Users/lopeze/Desktop/Statistics_R/Rings/D15-D16_medium/D16D17_medium.svg", width=10, height=8)
data$Strain_new= factor(data$Strain, levels=c("WT_12008", "WT_12253","Frq_12008", "Frq_12253", "Wc1_12253","Wc2_1", "Wc2_10"))
image=ggplot(data, aes(x=Conditions, fill=Conditions, y=Diameter)) +
stat_boxplot(geom="errorbar", size=0.5) +
     geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
     #labs(y=("Colony diameter(mm)")) +
     #facet_grid(~Conditions, scale="free") +
     facet_wrap(~Strain_new, scale="free") +
     #scale_fill_manual(values=c(WT_12008="#CC0099", WT_12253="#009900",Wc1_12253="#006600", Frq_12008="#0099CC", Frq_12253="#CC0000", Wc2_1="#FF9900", Wc2_10="#FFCC66"))+
     scale_fill_manual(values=c(PLYA="#CC3333", DOX="#009900",BMM="#FFCC33", MM="#CC0099"))+
     #scale_x_discrete(limits=c("WT_12008", "WT_12253","Frq_12008", "Frq_12253", "Wc1_12253","Wc2_1", "Wc2_10"))+
     scale_x_discrete(limits=c("PLYA","DOX","BMM", "MM"))+
     ylim(25,70)+
     ylim(25,70)+
     ylab("Colony diameter (mm)")+
     theme_bw() +
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=12),
           axis.text.x = element_text(colour="black", size=10 ,angle=45, hjust=1),
           axis.text.y = element_text(colour="black", size=10),
           axis.title.x = element_text(size=14),
           axis.title.y = element_text(size=14),
           strip.text = element_text(size=12),
           strip.background = element_rect(colour = "white"))
   plot(image)
   dev.off()

```
#Pathogenicity test E1, E2

```R   
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Pathogenicity/E1/E1_12253.csv")
attach(data)
boxplot(Score~Strain*Time, las = 2, cex.axis=0.6)

#ANOVA
options(max.print=1000000)

anv.model<-aov(Score~Strain*Time)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Strain", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Strain|Time, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#Standard Error
data_summary <- function(data, Score, Strain){
     require(plyr)
     summary_func <- function(x, col){
         c(mean = mean(x[[col]], na.rm=TRUE),
          se = sd(x[[col]]/sqrt(10), na.rm=TRUE))
     }
     data_sum<-ddply(data, Strain, .fun=summary_func, Score)
     data_sum <- rename(data_sum, c("mean" = Score))
     return(data_sum)
 }

df <- data_summary(data, Score="Score",
                  Strain=c("Strain", "Time"))

df$Strain=as.factor(df$Strain)
head(df)

#Ggplot curve

library(ggplot2)
library(Cairo)
ggplot(df, aes(x=Time, color=Strain, y=Score, shape=Strain)) + geom_point(size=5)
m4 <- lm(Score ~ Strain * Time, data=data)
data$y.hat <- predict(m4)

library(ggplot2)
library(Cairo)
svg("/Users/lopeze/Desktop/Statistics_R/Pathogenicity/E1/E1_12253.svg", width=10, height=8)
image = ggplot(df, aes(x=Time, color=Strain, shape=Strain, y=Score)) +
#scale_color_brewer(palette="Paired")+
scale_color_manual(values=c(WT_12253="#009900",Frq_12253="#CC0000",Wc1_12253="#FFCC00",Mock="#999999"))+
#scale_color_manual(values=c(WT_12008="#0099CC",Frq_12008="#CC0099",Wc2_12008="#FF9900",Mock="#999999"))+
geom_errorbar(aes(ymin=Score - se, ymax=Score + se), width=.1, size=0.6) +
geom_point(size=1) +
geom_line(aes(x=as.integer(Time)), size=0.6)+
ylab("Disease score") +
xlab("Days after inoculation")+
coord_fixed(ratio = 0.5)+
#ylim(0,10)+
scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
theme_bw() +
      theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.text.x = element_text(colour="black",size=12),
        axis.text.y = element_text(colour="black",size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
plot(image)
dev.off()
#AUDCP

sym<-read.table("/Users/lopeze/Desktop/Statistics_R/Pathogenicity/E3/E3_12253_audpc.txt", col.names=,1)
attach(sym)
sym[,1]=as.factor(sym[,1]) #Make leaf a factor  
attach(sym)

evaluation <- sym [ ,c (2,3,4,5,6)]
days<-c(0,7,14,21,28)
sym[,7]=audpc(evaluation,days)
sym
#write.table(sym, "/Users/lopeze/Desktop/Statistics_R/Pathogenicity/E2/audcp_E2_53.txt")

#ANOVA from AUDCP
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Pathogenicity/E3/audcp_E3_08.csv")
attach(data)
options(max.print=1000000)

anv.model<-aov(Value~Strain)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Strain", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Strain|Time, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups



```


#RING experiment

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/AFrq/Frq08_fix.csv")
attach(data)

data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/AFrq/Frq53_fix.csv")
attach(data)

data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/AWc1/Wc1_53_fix.csv")
attach(data)

data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/AWc2/Wc2_08_all.csv")
attach(data)

#ANOVA
anv.model<-aov(Diameter~Conditions*Strain)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

# Check normality with hist, qqplots
hist(Diameter)
shapiro.test(data$Diameter)
qqnorm(Diameter)
qqline(Diameter)

# If data isnt normal Run a Box-Cox proceedure to obtain optimal transformation
library("MASS")
boxcox(anv.model)
# Produces a plot of likelihood of the parameter lambda against values of lambda
# from -2 to 2
# Dotted vertical lines indicate the ideal value of lambda
# Refine range of lambda eg from 0 to 0.5 in increments of 0.1
boxcox(anv.model, lambda = seq(0, 2, 0.5))

# Plot boxcoxs
par(mfrow=c(1,1))
boxcox(anv.model)
boxcox(anv.model, lambda = seq(0, 0.5, 0.1))


# Create a variable with transformed values
Dm1<-(Diameter^0.05) #More normally distributed
Dm2<-(Diameter^1.27)
Dm3<-(Diameter^1.95)

# Check normality with hist, qqplots
par(mfrow=c(3,3))
hist(Dm1)
qqnorm(Dm1)
qqline(Dm1)
hist(Dm2)
qqnorm(Dm2)
qqline(Dm2)
hist(Dm3)
qqnorm(Dm3)
qqline(Dm3)

# Check normalitly with shapiro-wilks on all transformed data eg Ex1 Ex2
shapiro.test(Dm2)

# Run 2-way Anova with normalised data Ex2 was the best transformation
anv.model<-aov(Dm2~Experiment*Strain)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

anv.model<-aov(Dm1~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))


library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups


#Boxplot Afrq
library("ggplot2")
ggplot(data, aes(x=Strain, fill=Strain, y=Diameter)) +
stat_boxplot(geom="errorbar", size=0.5) +
     geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
     #labs(y=("Colony diameter(mm)")) +
     facet_grid(~Conditions, scale="free") +
     #scale_fill_manual(values=c(WT_12008="#0099CC",Frq_12008="#CC0099"))+
     #scale_fill_manual(values=c(WT_12253="#009900",Frq_12253="#CC0000"))+
     scale_fill_manual(values=c(WT_12253="#009900",Wc1_12253="#FFCC00"))+
     #scale_fill_manual(values=c(WT_12008="#CC0099",Wc2_1="#FF9900", Wc2_10="#FFCC66"))+
     #scale_fill_brewer(palette="Paired") +
     #scale_x_discrete(limits=c("WT_12008","Wc2_1","Wc2_10"))+
     #scale_x_discrete(limits=c("WT_12008","Frq_12008"))+
     #scale_x_discrete(limits=c("WT_12253","Frq_12253"))+
     scale_x_discrete(limits=c("WT_12253","Wc1_12253"))+
     ylim(20,60)+
     ylab("Colony diameter (mm)")+
     theme_bw() +
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=12),
           axis.text.x = element_text(colour="black", size=12),
           axis.text.y = element_text(colour="black", size=12),
           axis.title.x = element_text(size=14),
           axis.title.y = element_text(size=14),
           strip.text = element_text(size=12),
           strip.background = element_rect(colour = "white"))



           library("ggplot2")
           ggplot(data, aes(x=Strain, fill=Strain, y=Diameter)) +
           stat_boxplot(geom="bar", size=0.5) +
                geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
                #labs(y=("Colony diameter(mm)")) +
                facet_grid(~Conditions, scale="free") +
                scale_fill_brewer(palette="Paired") +
                #scale_fill_brewer(palette="Greens") +
                #scale_fill_brewer(palette="Paired") +
                ylim(20,60)+
                ylab("Colony diameter (mm)")+
                theme_bw() +
                theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      text = element_text(size=12),
                      axis.text.x = element_text(size=11, angle=45, hjust=1),
                      axis.text.y = element_text(size=11),
                      axis.title.x = element_text(size=15),
                      axis.title.y = element_text(size=15),
                      strip.text = element_text(face="bold", size=10),
                      strip.background = element_rect(colour="black",size=0.5))

#LIGHT PULSE experiment
```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Light_Pulse/LP0102.csv")
attach(data)

data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Light_Pulse/LP_Nc.csv")
attach(data)

#ANOVA
anv.model<-aov(Expression~Conditions*Strain)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

# Check normality with hist, qqplots
par(mfrow=c(4,2))
hist(Expression)
shapiro.test(data$Expression)
qqnorm(Expression)
qqline(Expression)

# If data isnt normal Run a Box-Cox proceedure to obtain optimal transformation
library("MASS")
boxcox(anv.model)
# Produces a plot of likelihood of the parameter lambda against values of lambda
# from -2 to 2
# Dotted vertical lines indicate the ideal value of lambda
# Refine range of lambda eg from 0 to 0.5 in increments of 0.1
boxcox(anv.model, lambda = seq(-2, -0.5, 0.5))

# Plot boxcoxs
par(mfrow=c(1,1))
boxcox(anv.model)
boxcox(anv.model, lambda = seq(0, 0.5, 0.1))

# Add data to original data set (Not done)
lamEx1<-cbind(anv.data, anv.data$Expression^0.17)
lamEx2<-cbind(anv.data, anv.data$Expression^0.26)
lamEx2<-cbind(anv.data, anv.data$Expression^0.26)

# Create a variable with transformed values
Ex1<-(Expression^-1.42)
Ex2<-(Expression^-0.62)
Ex3<-(Expression^-0.52)


# Check normality with hist, qqplots
par(mfrow=c(2,2))
hist(Ex1)
qqnorm(Ex1)
qqline(Ex1)
hist(Ex2)
qqnorm(Ex2)
qqline(Ex2)
hist(Ex3)
qqnorm(Ex3)
qqline(Ex3)

# Check normalitly with shapiro-wilks on all transformed data eg Ex1 Ex2
shapiro.test(Ex1)

# Run 2-way Anova with normalised data Ex2 was the best transformation
anv.model<-aov(Ex1~Experiment*Strain)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

anv.model<-aov(Ex1~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))


library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups




library("ggplot2")
library("Rmisc")
#calculate means
mn <- summarySE(data, measurevar="Expression", groupvars=c("Strain","Conditions"))
mn


ggplot(mn, aes(x=Conditions, fill=Conditions, y=Expression)) +
     geom_bar(stat="identity", width=0.5) +
     facet_grid(~Strain) +
     geom_bar(aes(ymin=Expression-se, ymax=Expression+se), width=.1) +
     theme_bw() +
     #scale_fill_brewer(palette="YIOrRd") +
     scale_fill_brewer(palette="Paired") +
     ylim(0,12)+
     ylab("Normalized relative expression (ddCt)")+
     coord_fixed(ratio = 0.2)+
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           #panel.border = element_rect(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=14),
           strip.text = element_text(face="bold", size=12),
           axis.text.x = element_text(colour="black", size=12),
           axis.text.y = element_text(colour="black", size=12),
           strip.background = element_rect(colour = "white"))



#TEMPERATURE PULSE experiment
 ```R
 data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Temp_Pulse/TP_2820.csv")
 attach(data)

 #ANOVA
 anv.model<-aov(Expression~Conditions*Gene)
 summary(anv.model)
 print(posthoc <- TukeyHSD(anv.model ))

 # Check normality with hist, qqplots
 par(mfrow=c(4,2))
 hist(Expression)
 shapiro.test(data$Expression)
 qqnorm(Expression)
 qqline(Expression)

 # If data isnt normal Run a Box-Cox proceedure to obtain optimal transformation
 library("MASS")
 boxcox(anv.model)
 # Produces a plot of likelihood of the parameter lambda against values of lambda
 # from -2 to 2
 # Dotted vertical lines indicate the ideal value of lambda
 # Refine range of lambda eg from 0 to 0.5 in increments of 0.1
 boxcox(anv.model, lambda = seq(-2, -0.5, 0.5))

 # Plot boxcoxs
 par(mfrow=c(1,1))
 boxcox(anv.model)
 boxcox(anv.model, lambda = seq(0, 0.5, 0.1))

 # Add data to original data set (Not done)
 lamEx1<-cbind(anv.data, anv.data$Expression^0.17)
 lamEx2<-cbind(anv.data, anv.data$Expression^0.26)
 lamEx2<-cbind(anv.data, anv.data$Expression^0.26)

 # Create a variable with transformed values
 Ex1<-(Expression^-1.42)
 Ex2<-(Expression^-0.52)
 Ex3<-(Expression^-0.52)


 # Check normality with hist, qqplots
 par(mfrow=c(2,2))
 hist(Ex1)
 qqnorm(Ex1)
 qqline(Ex1)
 hist(Ex2)
 qqnorm(Ex2)
 qqline(Ex2)
 hist(Ex3)
 qqnorm(Ex3)
 qqline(Ex3)

 # Check normalitly with shapiro-wilks on all transformed data eg Ex1 Ex2
 shapiro.test(Ex1)

 # Run 2-way Anova with normalised data Ex2 was the best transformatio

 anv.model<-aov(Ex2~Gene*Conditions)
 summary(anv.model)
 print(posthoc <- TukeyHSD(anv.model ))


 library(lsmeans)
 means<-(lsmeans(anv.model, pairwise~Conditions|Gene, adjust="tukey"))
 groups<-cld(means, alpha= .05)
 groups




 library("ggplot2")
 library("Rmisc")
 #calculate means
 mn <- summarySE(data, measurevar="Expression", groupvars=c("Gene","Conditions"))
 mn

mn$Conditions= factor(mn$Conditions, levels=c("28°C","20°C"), labels=c("28°C","20°C"))
mn$Gene_f = factor(mn$Gene, levels=c('frq','wc-1','wc-2','ccg-16'))
 ggplot(mn, aes(x=Conditions, fill=Conditions, y=Expression)) +
      geom_bar(stat="identity", width=0.5) +
      facet_grid(~Gene_f, switch="x") +
      geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=.1) +
      theme_bw() +
      #scale_fill_manual(breaks=c("28°C","20°C"), values=c("maroon4","lightskyblue1")) +
      scale_fill_manual(breaks=c("20°C","28°C"), values=c("lightskyblue1","maroon4")) +
      ylim(0,2)+
      ylab("Normalized relative expression (ddCt)")+
      #coord_fixed(ratio = 1)+
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size=12),
            strip.text = element_text(face="bold", size=10),
            #strip.background = element_rect(colour="black",size=1),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10))
```

#D5-D9 Menadione experiment

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Rings/D5-D9_Men/D9_53_data.csv")
attach(data)
boxplot(Diameter~Strain*Conditions, las = 2, cex.axis=0.6)

#ANOVA
anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Strain", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups

#Test for normality
res=resid(aov)
hist(res)
qqnorm(res)
shapiro.test(res) #p-value has to be over 0.5
#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)

#Standard Error
data_summary <- function(data, Diameter, Strain){
     require(plyr)
     summary_func <- function(x, col){
         c(mean = mean(x[[col]], na.rm=TRUE),
          se = sd(x[[col]]/sqrt(6), na.rm=TRUE))
     }
     data_sum<-ddply(data, Strain, .fun=summary_func, Diameter)
     data_sum <- rename(data_sum, c("mean" = Diameter))
     return(data_sum)
 }

df <- data_summary(data, Diameter="Diameter",
                  Strain=c("Strain", "Conditions"))

df$Strain=as.factor(df$Strain)
head(df)


#Plot results Bar Plot
library("ggplot2")
library("Rmisc")
library(Cairo)

#calculate means
#mn <- summarySE(data, measurevar="Diameter", groupvars=c("Strain","Conditions"))
#mn
svg("/Users/lopeze/Desktop/Statistics_R/Rings/D5-D9_Men/D5_D9_53_boarplot_ME.svg", width=10, height=8)
image=ggplot(df, aes(x=Strain, fill=Strain, y=Diameter)) +
     geom_bar(stat="identity", width=0.5) +
     facet_grid(~Conditions) +
     scale_x_discrete(limits=c("WT_12008","Frq_12008"))+
     #scale_x_discrete(limits=c("WT_12253","Wc1_12253"))+
     geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), width=.1) +
     scale_fill_manual(values=c(WT_12008="#CC0099", Frq_12008="#0099CC"))+
     #scale_fill_manual(values=c(WT_12253="#009900",Wc1_12253="#336666"))+
     theme_bw() +
     ylim(0,50)+
     ylab("Colony diameter (mm)")+
     coord_fixed(ratio = 0.2)+
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=12),
           #strip.text = element_text(face="bold", size=10),
           strip.background = element_rect(colour = "white"),
           axis.text.x = element_text(size=10, angle=45, hjust=1),
           axis.text.y = element_text(size=10))
         plot(image)
         dev.off()

   #Boxplot

   library("ggplot2")
   library(Cairo)
   svg("/Users/lopeze/Desktop/Statistics_R/Rings/D5-D9_Men/D9_53_boxplot_ME.svg", width=10, height=8)
   image=ggplot(data, aes(x=Strain, fill=Strain, y=Diameter)) +
   stat_boxplot(geom="errorbar", size=0.5) +
        geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
        facet_grid(~Conditions) +
        #scale_fill_manual(values=c(WT_12008="#CC0099", Frq_12008="#0099CC"))+
        scale_fill_manual(values=c(WT_12253="#009900",Wc1_12253="#336666"))+
        #scale_x_discrete(limits=c("WT_12008","Frq_12008"))+
        scale_x_discrete(limits=c("WT_12253","Wc1_12253"))+
        ylim(15,60)+
        ylab("Colony diameter (mm)")+
        coord_fixed(ratio = 0.2)+
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              text = element_text(size=12),
              axis.text.x = element_text(colour="black", size=10 ,angle=45, hjust=1),
              axis.text.y = element_text(colour="black", size=10),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              strip.text = element_text(size=11),
              strip.background = element_rect(colour = "white"))
      plot(image)
      dev.off()
