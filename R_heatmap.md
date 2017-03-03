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
pdf(file = "/Users/lopeze/Desktop/Bioinformatics/promoter-heatmap.pdf",width=7,height=7)
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
