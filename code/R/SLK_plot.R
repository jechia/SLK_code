library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)

setwd("D:/2021/SLK")


quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

trim_ENSG<-function(x){
  ENSG.l<-Reduce(c,lapply(unique(x),
                          function(x){
                            strsplit(x,split=".",fixed = T)[[1]][1]
                          }))
}

sample_id<-c("SLK_L_1","SLK_L_2",
             "SLK_S_1","SLK_S_2",
             "SLK_WT_1","SLK_WT_2",
             "SLK_KO_1","SLK_KO_2")
files<-paste0(sample_id,".genes.results")
names(files)<-c("SLK_L_1","SLK_L_2",
                "SLK_S_1","SLK_S_2",
                "SLK_WT_1","SLK_WT_2",
                "SLK_KO_1","SLK_KO_2")
txi.rsem<- tximport(files, type = "rsem", txIn = F, txOut = F)
sample_condition <- data.frame(condition=factor(c(rep("SLK_L",2),
                                                  rep("SLK_S",2),
                                                  rep("SLK_WT",2),
                                                  rep("SLK_KO",2)),
                                                levels=c("SLK_WT","SLK_S",
                                                         "SLK_L","SLK_KO")))

txi.rsem$length[txi.rsem$length==0]=1
dds <- DESeqDataSetFromTximport(txi.rsem, 
                                colData=sample_condition,
                                design=~condition)

keep <- rowMeans(counts(dds))>50
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds<- DESeq(dds)

new_count.m<-counts(dds,normalized=T)

deg.pca<-prcomp(t(new_count.m))
plot.df <- cbind(deg.pca$x[,1:2], rep(c("L","S","WT","KO"),each=2)) %>% 
  as.data.frame()
plot.df$PC1 <- as.numeric(plot.df$PC1) / (deg.pca$sdev[1] * sqrt(ncol(new_count.m)))
plot.df$PC2 <- as.numeric(plot.df$PC2) / (deg.pca$sdev[2] * sqrt(ncol(new_count.m)))
plot.df$class <- as.factor(plot.df$V3)
ve <- deg.pca$sdev^2 / sum(deg.pca$sdev^2)
ve <- ve[c(1, 2)]
labs <- paste0(PC, " (", round(ve * 100, 2), "%)")
p<-ggplot(plot.df, aes(PC1, PC2, colour = class)) +
  geom_point(size = 1.8)+
  scale_color_manual(name="Groups", values=c("#cccccc",
                                             "#ff9966", 
                                             "#a8d7ff",
                                             "#797979"))+
  xlab(labs[1])+
  ylab(labs[2])+
  theme_classic()
ggsave(p,file="PCA_7.7.pdf",device="pdf",height = 2,width = 3)


cluster.m<-new_count.m

row.idx<-head(order(deg.pca[["rotation"]][,1],decreasing = T),500)


mat_breaks <- quantile_breaks(log(cluster.m+1), n = 100)

grad_cols<-colorRampPalette(c("#d42517","#f28e16", "#fffbfa"))
#5*2.76
pheatmap(t(log(cluster.m+1)),breaks=mat_breaks,
         show_colnames = F,
         color=rev(grad_cols(100)))


L_res <- results(dds, contrast=c("condition","SLK_L","SLK_S"))
WT_res <- results(dds, contrast=c("condition","SLK_WT","SLK_S"))

L_res.m <- as.data.frame(L_res)
WT_res.m <- as.data.frame(WT_res)

plot.df<-merge(L_res.m,WT_res.m,by="row.names")

plot.df$status <- "NS"
plot.df$status[plot.df$padj.x<0.05&plot.df$padj.y<0.05&plot.df$log2FoldChange.x>0&plot.df$log2FoldChange.y>0] <- "up"
plot.df$status[plot.df$padj.x<0.05&plot.df$padj.y<0.05&plot.df$log2FoldChange.x<0&plot.df$log2FoldChange.y<0] <- "down"

p<-ggplot(plot.df,aes(x=log2FoldChange.x,
                   y=log2FoldChange.y))+
  geom_point(aes(color=status))+
  scale_color_manual(labels=c("down", "NS", "up"), 
                     values=c("#4e8bba", "grey", "#db524f"))+
  stat_cor(aes(label = ..r.label..),method = "pearson", label.x = -10, label.y = 4.5)+
  theme_classic()+
  labs(x="L vs S",y="WT vs S")
ggsave(p,file="L_WT_cor_7.7.pdf",device="pdf",height = 2,width = 3)


L_ENSG.l<-trim_ENSG(rownames(subset(L_res.m,padj<0.05)))
ego <- clusterProfiler::enrichGO(gene          = unique(L_ENSG.l),
                                   keyType = "ENSEMBL",
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
  
L_GO.df<-ego@result
L_GO.df$ratio<-L_GO.df$Count/2979

GO.l<-c("Wnt signaling pathway","cell-substrate adhesion",
        "negative regulation of growth",
        "regulation of cell growth",
        "actin filament organization",
        "epithelial cell development",
        "cellular response to epidermal growth factor stimulu",
        "regulation of cell morphogenesis",
        "regulation of epithelial cell proliferation",
        "mesenchyme development")
plot_go.df<-subset(L_GO.df,Description%in%GO.l)
plot_go.df<-plot_go.df[order(plot_go.df$ratio,decreasing = F),]
plot_go.df$Description<-factor(plot_go.df$Description,levels=plot_go.df$Description)
p<-ggplot(plot_go.df,
       aes(x=ratio,y=Description,color=p.adjust,size=Count))+geom_point()+
  scale_colour_gradientn(colours=c("#ac4f41","#db524f","#efefef","#4e8bba"),
                         breaks=c(0.000001,0.0001,0.05,1),
                         limits=c(0.000001,1),trans="log2")+
  scale_size(range = c(3,5))+
  theme_classic()+
  scale_x_continuous(limits = c(0.018,0.04),breaks = c(0.02, 0.03, 0.04))+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.text.y = element_text(size=12, color="black"))+
  labs(x="Tissues",y="")
ggsave(p,file="GO_L_S_7.17.pdf",device="pdf",height = 3.5,width = 6)

L_gene.df<-data.frame(ENSG=rownames(L_res.m),FC=L_res.m$log2FoldChange)
ENSG.l<-Reduce(c,lapply(unique(L_gene.df$ENSG),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
L_gene.df$ENSG<-ENSG.l
gene.l<-L_gene.df$FC
names(gene.l)<-L_gene.df$ENSG
gene.l<-gene.l[order(gene.l,decreasing=T)]
L_gene.df<-L_gene.df[order(L_gene.df$FC,decreasing=T),]

write.table(unique(L_gene.df),"SLK_gsea_L_S.rnk",quote=F,col.names = F,row.names=F,sep="\t")


L_res <- results(dds, contrast=c("condition","SLK_L","SLK_WT"))
S_res <- results(dds, contrast=c("condition","SLK_S","SLK_WT"))
KO_res <- results(dds, contrast=c("condition","SLK_KO","SLK_WT"))

L_res.m <- as.data.frame(L_res)
S_res.m <- as.data.frame(S_res)
KO_res.m <- as.data.frame(KO_res)



heatmap.l<-Reduce(intersect,list(rownames(subset(L_res.m,padj<0.05)),
                                 rownames(subset(S_res.m,padj<0.05)),
                                 rownames(subset(KO_res.m,padj<0.05))))

L_heatmap.m<-subset(L_res.m,rownames(L_res.m)%in%heatmap.l)[,c(2,6)]
S_heatmap.m<-subset(S_res.m,rownames(S_res.m)%in%heatmap.l)[,c(2,6)]
KO_heatmap.m<-subset(KO_res.m,rownames(KO_res.m)%in%heatmap.l)[,c(2,6)]
tmp1<-merge(L_heatmap.m,S_heatmap.m,by.x="row.names",by.y="row.names")
heatmap_lfc.df<-merge(tmp1,KO_heatmap.m,by.x="Row.names",by.y="row.names")
heatmap_lfc.m<-heatmap_lfc.df[,c(2,4,6)]
rownames(heatmap_lfc.m)<-heatmap_lfc.df[,1]

colnames(heatmap_lfc.m)<-c("SLK_L","SLK_S","SLK_KO")

heatmap_lfc.m[heatmap_lfc.m>2]=2
heatmap_lfc.m[heatmap_lfc.m < -2]= -2
mat_breaks <- quantile_breaks(heatmap_lfc.m[,1], n = 100)
pheatmap(heatmap_lfc.m,breaks = mat_breaks,
         cluster_cols = T,show_rownames = F)

L_ENSG.l<-trim_ENSG(rownames(subset(L_res.m,padj<0.05)))

S_ENSG.l<-trim_ENSG(rownames(subset(S_res.m,padj<0.05)))

KO_ENSG.l<-trim_ENSG(rownames(subset(KO_res.m,padj<0.05)))

total.l<-list(L_ENSG.l,
              S_ENSG.l,
              KO_ENSG.l)

p<-venn.diagram(
  x = total.l,
  category.names = c("L vs WT" , "S vs WT","KO vs WT"),
  filename=NULL,
  output = F ,
  height = 350, 
  width = 350 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff','#fde725ff'),
  fill = c("#440154ff",'#21908dff','#fde725ff'),
  alpha=0.3,
  cex = 0.4,
  fontfamily = "sans",
  cat.cex = 0.45,
  cat.default.pos = "outer",
  cat.pos = c(-10,10,180),
  cat.dist = c(0.05, 0.05,0.05),
  cat.fontfamily = "sans",
  cat.col = c('black', 'black','black')
)
ggsave(p,file="venn_3_7.9.pdf",device="pdf",height = 1,width = 1)

ego <- clusterProfiler::enrichGO(gene          = unique(KO_ENSG.l),
                                 keyType = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)

KO_GO.df<-ego@result
KO_GO.df$ratio<-KO_GO.df$Count/4893

GO.l<-c("Wnt signaling pathway","cytoplasmic translation",
        "chromosome segregation",
        "regulation of nuclear division",
        "DNA replication",
        "establishment of organelle localization",
        "RNA localization",
        "regulation of cell morphogenesis",
        "mRNA transport",
        "ribonucleoprotein complex assembly")
plot_go.df<-subset(KO_GO.df,Description%in%GO.l)
plot_go.df<-plot_go.df[order(plot_go.df$ratio,decreasing = F),]
plot_go.df$Description<-factor(plot_go.df$Description,levels=plot_go.df$Description)
plot_go.df$p.adjust[plot_go.df$p.adjust<10^-10]<-10^-10
p<-ggplot(plot_go.df,
          aes(x=ratio,y=Description,color=p.adjust,size=Count))+geom_point()+
  scale_colour_gradientn(colours=c("#ac4f41","#db524f","#efefef","#4e8bba"),
                         breaks=c(10^-10,0.0001,0.05,1),
                         limits=c(10^-10,1),trans="log2")+
  scale_size(range = c(3,5))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.text.y = element_text(size=12, color="black"))+
  labs(x="Tissues",y="")
ggsave(p,file="GO_KO_WT_7.9.pdf",device="pdf",height = 3.5,width = 6)
