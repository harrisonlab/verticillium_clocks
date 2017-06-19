Michelle script:

library(gplots)
virulence<-read.table("/Users/hulinm/Documents/effectors_reg_heatmap.txt", col.names=,1, row.names=1, check.names=FALSE)
virulence_matrix<-data.matrix(virulence)
scale <- colorRampPalette(c("white", "yellow", "forestgreen"), space = "rgb")(100)
pdf(file = "/Users/hulinm/Documents/e-heatmap.pdf",width=7,height=7)
virulence_heatmap <- heatmap.2(virulence_matrix, margins = c(10, 18), Rowv=NA, cexCol=0.3, cexRow=0.4 ,lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ), trace="none", hline=NULL, vline=NULL, tracecol="Gray", col=scale)
dev.off()


Script used:

install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs-heatmap.pdf",width=7,height=7)
orthologs_heatmap <- heatmap.2(orthologs_matrix, margins = c(6, 12), dendogram="column", Colv=FALSE, cexCol=0.8, cexRow=0.8 , trace="none", tracecol="Gray", col=scale)
dev.off()


Alternative script:
/Users/lopeze/Desktop/Bioinformatics

library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_genes.txt", col.names=1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "green"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs-heatmap.pdf",width=7,height=7)
orthologs_heatmap <- heatmap.2(orthologs_matrix, geom_tile(aes(fill = scale)), margins = c(5, 12), Colv=FALSE, dendrogram="row", cexCol=0.8, cexRow=0.8, trace="none", tracecol="Gray", col=scale)
dev.off()

Alternative script using Pheatmap:

library("pheatmap")
library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs-heatmap.pdf",width=7,height=7)
mylwid = c(1.5,4,0.5)
mylhei = c(1.5,4,1)
mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
orthologs_heatmap <- pheatmap(orthologs_matrix, lmat=mylmat, lwid=mylwid, lhei=mylhei, margins = c(6, 12), dendogram="column", Colv=FALSE, cexCol=0.8, cexRow=0.8 , trace="none", tracecol="Gray", col=scale)
dev.off()


==========
Promoter motif heatmap
==========

library("pheatmap")
library("gplots")
promoter<-read.table("/Users/lopeze/Desktop/Bioinformatics/Promotermotif.txt", col.names=,1, row.names=1, check.names=FALSE)
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
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/Promoter_motif/promoter-heatmap_DH.pdf",width=7,height=7)
promoter-heatmap <- heatmap.2(promoter_matrix,
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
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/D8-D10_V.spp/D8-D10_data.csv")
attach(data)
boxplot(Diameter~Strain*Conditions)
hist(Diameter)
shapiro.test(data$Diameter) #p-value has to be over 0.5

#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)

#ANOVA
anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(aov.model, pairwise~Conditions|Strain, adjust="tukey")
groups<-cld(means, alpha= .05)
```

#D5-D9 Menadione experiment

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/D5-D9_Men/D9_53_data.csv")
attach(data)
boxplot(Diameter~Strain*Conditions, las = 2, cex.axis=0.6)


hist(Diameter)
shapiro.test(data$Diameter) #p-value has to be over 0.5

#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)

#ANOVA
anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(aov.model, pairwise~Conditions|Strain, adjust="tukey")
groups<-cld(means, alpha= .05)
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

#D15 medium

```R   
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/D15-D16_medium/D15_data/D15_med_data.csv")
attach(data)
boxplot(Diameter~Strain*Conditions, las = 2, cex.axis=0.6)


hist(Diameter)
shapiro.test(data$Diameter) #p-value has to be over 0.5

#If the data is not normally distributed:
data.log<-log(data$Diameter)
shapiro.test(data.log$Diameter)
hist(data.log)

#ANOVA
options(max.print=1000000)

anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

```
#Pathogenicity test E1

```R   
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Pathogenicity/E1/E1_53_data.csv")
attach(data)

#ANOVA
options(max.print=1000000)

anv.model<-aov(Score~Strain*Time)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Strain", group=TRUE)
H

#Standard error
df <- with(data , aggregate(Score, list(Strain=Strain, Time=Time), mean))
df$se <- with(data , aggregate(Score, list(Strain=Strain, Time=Time),
              function(x) sd(x)/sqrt(10)))

#-------------#

data_summary <- function(data, Score, Strain){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]]/sqrt(10), na.rm=TRUE))
  }
  data_sum<-ddply(data, Strain, .fun=summary_func,
                  Score)
  data_sum <- rename(data_sum, c("mean" = Score))
 return(data_sum)
}

df <- data_summary(data, Score="Score",
                    Strain=c("Strain", "Time"))

df$Strain=as.factor(df$Strain)
head(df)


#Standard deviation
data_summary <- function(data, Score, Strain){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, Strain, .fun=summary_func,
                  Score)
  data_sum <- rename(data_sum, c("mean" = Score))
 return(data_sum)
}

df2 <- data_summary(data, Score="Score",
                    Strain=c("Strain", "Time"))

df2$Strain=as.factor(df2$Strain)
head(df2)

#Line plot with error bars
library("ggplot2")
m4 <- lm(Score ~ Strain * Time, data=data)
dats$y.hat <- predict(m4)
ggplot(df, aes(x=Time, color=Strain, y=Score)) +
    geom_point(size=1) +
    geom_line(aes(y=Score, x=as.integer(Time))) +
    labs(x="Days after inoculation", y=("Disease score")) +
    geom_errorbar(aes(ymin=Score - df$se, ymax=Score + df$se), width=.2,
                 position=position_dodge(0)) +
    ylim(0,10)+
    scale_fill_manual(values=c("blue", "red", "green", "grey")) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=12),
        aspect.ratio=1/1,
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))


#Area under the curve
days<-c(15,23,30,37,44,51,58)

AGNd08[,13]=audpc(AGNd08[,6:12],days)
AGNd08[,14]=audpc(AGNd08[,6:12],days,“relative”)
colnames(AGNd08)[13]=c(‘audpc’)
colnames(AGNd08)[14]=c(‘audpc_r’)

```
