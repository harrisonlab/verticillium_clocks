Michelle script:

library(gplots)
virulence<-read.table("/Users/hulinm/Documents/effectors_reg_heatmap.txt", col.names=,1, row.names=1, check.names=FALSE)
virulence_matrix<-data.matrix(virulence)
scale <- colorRampPalette(c("white", "yellow", "forestgreen"), space = "rgb")(100)
pdf(file = "/Users/hulinm/Documents/e-heatmap.pdf",width=7,height=7)
virulence_heatmap <- heatmap.2(virulence_matrix, margins = c(10, 18), Rowv=NA, cexCol=0.3, cexRow=0.4 ,lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25 ), trace="none", hline=NULL, vline=NULL, tracecol="Gray", col=scale)
dev.off()


Script used:
```R
install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

library("gplots")
orthologs<-read.table("/Users/lopeze/Desktop/Bioinformatics/Ortholog_genes.txt", col.names=,1, row.names=1, check.names=FALSE)
orthologs_matrix<-data.matrix(orthologs)
scale <- colorRampPalette(c("white", "forestgreen", "darkgreen"), space = "rgb")(100)
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/orthologs-heatmap.pdf",width=7,height=7)
orthologs_heatmap <- heatmap.2(orthologs_matrix, margins = c(6, 12), dendogram="column", Colv=FALSE, cexCol=0.8, cexRow=0.8 , trace="none", tracecol="Gray", col=scale)
dev.off()
```

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
boxplot(Score~Strain*Time, las = 2, cex.axis=0.6)

#ANOVA
options(max.print=1000000)

anv.model<-aov(Score~Strain*Time)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Strain", group=TRUE)
H

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
ggplot(df, aes(x=Time, color=Strain, y=Score, shape=Strain)) + geom_point(size=5)

m4 <- lm(Score ~ Strain * Time, data=data)
data$y.hat <- predict(m4)

ggplot(df, aes(x=Time, color=Strain, shape=Strain, y=Score)) +
geom_errorbar(aes(ymin=Score - se, ymax=Score + se), width=.1, size=0.6) +
geom_point(size=1) +
geom_line(aes(x=as.integer(Time)), size=0.6)+
ylab("Disease score") +
xlab("Days after inoculation")+
coord_fixed(ratio = 0.5)+
ylim(0,10)+
scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
theme_bw() +
      theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))




```
#RING experiment

```R
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Light/Interspecies/Intersp.csv")
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


#Boxplot
library("ggplot2")
ggplot(data, aes(x=Conditions, fill=Conditions, y=Diameter)) +
stat_boxplot(geom="errorbar", size=0.5) +
     geom_boxplot(outlier.shape=16, outlier.size=1, fatten=1) +
     #labs(y=("Colony diameter(mm)")) +
     facet_grid(~Strain, scale="free") +
     scale_fill_brewer(palette="Paired") +
     #scale_fill_manual(breaks=c("WT_12008","ΔFrq_12008"), values=c("royalblue","steelblue1")) +
     #scale_fill_manual(breaks=c("WT_12253","ΔWc1_12253"), values=c("firebrick","antiquewhite")) +
     #scale_fill_manual(breaks=c("WT_12008","ΔWc2_12008"), values=c("royalblue","azure3")) +
     ylim(20,55)+
     ylab("Colony diameter (mm)")+
     theme_bw() +
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           text = element_text(size=12),
           axis.text.x = element_text(size=11),
           axis.text.y = element_text(size=11),
           axis.title.x = element_text(size=15),
           axis.title.y = element_text(size=15),
           strip.text = element_text(face="bold", size=10),
           strip.background = element_rect(colour="black",size=0.5))



           library("ggplot2")
           ggplot(data, aes(x=Strain, fill=Strain, y=Diameter)) +
           stat_boxplot(geom="errorbar", size=0.5) +
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
data<-read.csv("/Users/lopeze/Desktop/Statistics_R/Light_Pulse/LP_0102.csv")
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
     geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=.1) +
     theme_bw() +
     #scale_fill_brewer(palette="YIOrRd") +
     scale_fill_brewer(palette="Paired") +
     ylim(0,12)+
     ylab("Normalized relative expression (ddCt)")+
     coord_fixed(ratio = 1)+
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_rect(),
           panel.background = element_blank(),
           text = element_text(size=12),
           strip.text = element_text(face="bold", size=10),
           strip.background = element_rect(colour="black",size=1),
           axis.text.x = element_text(size=10),
           axis.text.y = element_text(size=10))


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


hist(Diameter)
shapiro.test(data$Diameter) #p-value has to be over 0.5


#ANOVA
anv.model<-aov(Diameter~Strain*Conditions)
summary(anv.model)
print(posthoc <- TukeyHSD(anv.model ))

library(agricolae)
H<-HSD.test(anv.model, "Conditions", group=TRUE)
H

library(lsmeans)
means<-(lsmeans(anv.model, pairwise~Conditions|Strain, adjust="tukey"))
groups<-cld(means, alpha= .05)
groups


#Plot results Bar Plot
library("ggplot2")
library("Rmisc")
#calculate means
mn <- summarySE(data, measurevar="Diameter", groupvars=c("Strain","Conditions"))
mn

#mn$Conditions= factor(mn$Conditions, levels=c("28°C","20°C"), labels=c("28°C","20°C"))
#mn$Gene_f = factor(mn$Gene, levels=c('frq','wc-1','wc-2','ccg-16'))

ggplot(mn, aes(x=Strain, fill=Strain, y=Diameter)) +
     geom_bar(stat="identity", width=0.5) +
     facet_grid(~Conditions) +
     geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), width=.1) +
     theme_bw() +
     #scale_fill_manual(breaks=c("28°C","20°C"), values=c("maroon4","lightskyblue1")) +
     #scale_fill_manual(breaks=c("20°C","28°C"), values=c("lightskyblue1","maroon4")) +
     ylim(0,60)+
     ylab("Colony diameter (mm)")+
     #coord_fixed(ratio = 1)+
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
