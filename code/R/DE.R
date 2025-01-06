library(NOISeq)
library(tximport)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

setwd("D:/2021/SLK")

sample_name<-c("HCT116-WT","HCT116-L","HCT116-S","HCT116-KO")
files<-paste0(sample_name,".genes.results")
names(files)<-sample_name
SLK_rsem <- tximport(files, type = "rsem", txIn = F,txOut = F)
rsem_count.m<-SLK_rsem$counts

myfactors = data.frame(condition= c("WT", "L", "S", "KO"))

mylengths = rowMax(SLK_rsem[["length"]])
names(mylengths)<-rownames(SLK_rsem[["length"]])

myTMM = tmm(rsem_count.m, long = mylengths, lc = 0)
idx<-rowMeans(myTMM)>50

mydata <- readData(data = myTMM[idx,], factors = myfactors,length=mylengths[idx])

L_results <- noiseq(mydata, factor = "condition", conditions=c("L","WT"),
                    k = NULL, norm = "n", pnr = 0.2,
                    nss = 5, v = 0.02, lc = 1, replicates = "no")
L_deg = subset(degenes(L_results, q = 0.8, M = NULL),abs(M)>0.5)

#DE.plot(L_results, q = 0.8, graphic = "MD")

S_results <- noiseq(mydata, factor = "condition", conditions=c("S","WT"),
                    k = NULL, norm = "n", pnr = 0.2,
                    nss = 5, v = 0.02, lc = 1, replicates = "no")
S_deg = subset(degenes(S_results, q = 0.8, M = NULL),abs(M)>0.5)


L_S_results <- noiseq(mydata, factor = "condition", conditions=c("L","S"),
                    k = NULL, norm = "n", pnr = 0.2,
                    nss = 5, v = 0.02, lc = 1, replicates = "no")
L_S_deg = subset(degenes(L_S_results, q = 0.8, M = NULL),abs(M)>0.5)

KO_results <- noiseq(mydata, factor = "condition", conditions=c("KO","WT"),
                     k = NULL, norm = "n", pnr = 0.2,
                     nss = 5, v = 0.02, lc = 1, replicates = "no")
KO_deg = subset(degenes(KO_results, q = 0.8, M = NULL),abs(M)>0.5)

library(VennDiagram)
total.l<-list(rownames(L_deg),
              rownames(S_deg),
              rownames(L_S_deg),
              rownames(KO_deg))

p<-venn.diagram(
  x = total.l,
  category.names = c("L vs WT" , "S vs WT","L vs S","KO vs WT"),
  filename=NULL,
  output = F ,
  height = 350, 
  width = 350 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff','#fde725ff','#f28e16'),
  fill = c("#440154ff",'#21908dff','#fde725ff','#f28e16'),
  alpha=0.3,
  cex = 0.4,
  fontfamily = "sans",
  cat.cex = 0.45,
  cat.default.pos = "outer",
  cat.pos = c(-10,10,-20,20),
  cat.dist = c(0.2, 0.2,0.12,0.11),
  cat.fontfamily = "sans",
  cat.col = c('black', 'black','black','black')
)
ggsave(p,file="venn_4_3.27.pdf",device="pdf",height = 1.5,width = 1.5)


total.l<-list(rownames(L_deg),
              rownames(S_deg),
              rownames(KO_deg))

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
ggsave(p,file="venn_3_4.27.pdf",device="pdf",height = 1.5,width = 1.5)

pie.df<-data.frame(value=c(486,383,156,152,76,416,261,138),
                   group=c("other","L","KO & L","KO & L & S",
                           "S & L","S","S & KO","KO"))
pie.df<-pie.df %>% mutate(freq=value/sum(value))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
ggplot(pie.df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(label = value),
            position = position_stack(vjust=0.5),
            size=5)
pie.df$group<-factor(pie.df$group,levels=c("KO & L & S","KO & L","S & L",
                                           "S & KO","KO","S","L","other"))
p<-ggplot(pie.df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity",position="stack")+
  blank_theme +
  scale_fill_manual(values=c("#8dd3c7",
                             "#ffffb3",
                             "#bebada", 
                             "#fb8072",
                             "#80b1d3",
                             "#fdb462",
                             "#b3de69",
                             "#fccde5"))+
  theme(axis.text.x=element_blank()) +
  geom_text(aes(label = value),
            position = position_stack(vjust=0.5),
            size=5)
ggsave(p,file="barplot_4.27.pdf",device="pdf",height = 3,width = 2.5)
            #paste0(format(freq*100,digits=2), "%")

inter_gene.l<-Reduce(union,list(intersect(rownames(L_deg),
                                          rownames(L_S_deg)),
                                intersect(rownames(S_deg),
                                          rownames(L_S_deg)),
                                intersect(rownames(KO_deg),
                                          rownames(L_S_deg))))

#inter_gene.l<-Reduce(intersect,list(rownames(L_deg),
#                     rownames(S_deg),
#                     rownames(L_S_deg)))

sel.idx<-rownames(myTMM)%in%inter_gene.l
library(pheatmap)
library(biomaRt)

sel_ENSG.df<-data.frame(ENSG=c("ENSG00000124216","ENSG00000122691",
                               "ENSG00000039068","ENSG00000170558",
                               "ENSG00000147065","ENSG00000104067",
                               "ENSG00000026025","ENSG00000115414",
                               "ENSG00000163347","ENSG00000197822"),
                        name=c("SNAI1","TWIST1","CDH1","CDH2","MSN",
                               "TJP1","VIM","FN1","CLDN1","OCLN"))

sel_heat.df<-as.data.frame(myTMM[,1:4])
sel_heat.df$ENSG<-Reduce(c,lapply(rownames(sel_heat.df),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
sel_merge.df<-unique(merge(sel_heat.df,sel_ENSG.df,by.x="ENSG",by.y="ENSG"))

sel_plot_heat.df<-sel_merge.df[,2:5]
rownames(sel_plot_heat.df)<-sel_merge.df[,6]
sel_mat_breaks <- quantile_breaks(log(sel_plot_heat.df$`HCT116-S`+1), n = 100)

write.table(sel_plot_heat.df,"small_heatmap.tsv",sep="\t",row.names = T,
            col.names = T,quote=F)
colnames(sel_plot_heat.df)<-c("SLK-WT","SLK-L","SLK-S","SLK-KO")
pheatmap(log(sel_plot_heat.df+1),
         breaks = sel_mat_breaks,
         na_col = "lightgrey",color=rev(grad_cols(100)),
         border_color = NA)

heatmap_ENSG.l<-Reduce(c,lapply(unique(inter_gene.l),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl",host = "https://useast.ensembl.org")
sel_gs <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),
                          filters = 'ensembl_gene_id', values =heatmap_ENSG.l , mart = ensembl)

Wnt.data<-AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0016055", columns="ENSEMBL")
epi.data <-AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0002064", columns="ENSEMBL")
junction.data <-AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0034329", columns="ENSEMBL")
mesen.data<-AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0060485", columns="ENSEMBL")
  

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
heat.df<-as.data.frame(myTMM[sel.idx,1:4])
heat.df$ENSG<-Reduce(c,lapply(rownames(heat.df),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
merge.df<-unique(merge(heat.df,sel_gs,by.x="ENSG",by.y="ensembl_gene_id"))
path_gene.l<-Reduce(c,list(Wnt.data$ENSEMBL,
                           epi.data$ENSEMBL,
                           junction.data$ENSEMBL,
                           mesen.data$ENSEMBL))
idx<-which(merge.df$ENSG%in%path_gene.l)
plot_heat.df<-merge.df[idx,2:5]
rownames(plot_heat.df)<-merge.df[idx,6]

write.table(plot_heat.df,"large_heatmap.tsv",sep="\t",row.names = T,
            col.names = T,quote=F)
mat_breaks <- quantile_breaks(log(plot_heat.df$`HCT116-S`+1), n = 100)


annotation_col = data.frame(
  Epi = as.character(as.numeric(merge.df[idx,1]%in%epi.data$ENSEMBL)),
  Mesen = as.character(as.numeric(merge.df[idx,1]%in%mesen.data$ENSEMBL)),
  Junc = as.character(as.numeric(merge.df[idx,1]%in%junction.data$ENSEMBL)),
  Wnt = as.character(as.numeric(merge.df[idx,1]%in%Wnt.data$ENSEMBL))
)
rownames(annotation_col) = rownames(plot_heat.df)

#pheatmap(plot_heat.df,breaks = mat_breaks)

ann_colors = list(
  Epi = c("1"="#61649f", "0"="#f8f8fb"),
  Junc = c("1"="#ee8055", "0"="#fef5f2"),
  Mesen = c("1"="#f9d367", "0"="#fefcf3"),
  Wnt = c("1" = "#add5a2", "0"="#fbfdfb")
)

grad_cols<-colorRampPalette(c("#d42517","#f28e16", "#fffbfa"),interpolate = c("linear"))

sort.idx<-rep(4,nrow(plot_heat.df))
sort.idx[merge.df[idx,1]%in%mesen.data$ENSEMBL]<-3
sort.idx[merge.df[idx,1]%in%junction.data$ENSEMBL]<-2
sort.idx[merge.df[idx,1]%in%Wnt.data$ENSEMBL]<-1
sort.df<-data.frame(idx=sort.idx,count=rowMeans(plot_heat.df))

with(sort.df, order(idx,count))
library(pheatmap)
pheatmap(log(plot_heat.df[with(sort.df, order(idx,count)),]+1),
         breaks = mat_breaks,
         cluster_rows = F,
                 na_col = "lightgrey",color=rev(grad_cols(100)),
                 annotation_row = annotation_col,
                 annotation_colors = ann_colors,
         border_color = NA)


getGOdf<-function(ENSG.l,type){
  ego <- clusterProfiler::enrichGO(gene          = unique(ENSG.l),
                                   keyType = "ENSEMBL",
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
  
  result.m<-ego@result
  GO.df<-data.frame(d=result.m[["Description"]],
                    c=result.m[["Count"]],
                    p=result.m[["p.adjust"]],
                    type=rep(type,length(result.m[["p.adjust"]])))
  return(GO.df)
}

L_S_ENSG.l<-Reduce(c,lapply(unique(rownames(L_S_deg)),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))

L_S_GO.df<-getGOdf(L_S_ENSG.l,"L Vs S")

L_ENSG.l<-Reduce(c,lapply(unique(rownames(L_deg)),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
L_GO.df<-getGOdf(L_ENSG.l,"L Vs WT")

S_ENSG.l<-Reduce(c,lapply(unique(rownames(S_deg)),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
S_GO.df<-getGOdf(S_ENSG.l,"S Vs WT")

KO_ENSG.l<-Reduce(c,lapply(unique(rownames(KO_deg)),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
KO_GO.df<-getGOdf(KO_ENSG.l,"KO Vs WT")

total_GO.df<-rbind(L_S_GO.df,L_GO.df,S_GO.df,KO_GO.df)

new_total_GO.df<-rbind(L_GO.df,S_GO.df,KO_GO.df)

GO_count<-total_GO.df%>%group_by(d)%>%summarise(n=n(),mcount=mean(c),count=max(c))
sig_GO_count<-total_GO.df%>%filter(p<0.05)%>%group_by(d)%>%summarise(n=n(),mcount=mean(c),count=max(c),mp=max(p))



GO.l<-c("Wnt signaling pathway","tissue migration",
        "epithelial cell migration",
        "regulation of epithelial cell proliferation",
        "mesenchyme development",
        "epithelial cell development",
        "cell junction assembly",
        "mesenchymal cell differentiation",
        "embryonic skeletal system development",
        "regulation of cell growth")
sel_GO.df<-subset(total_GO.df,d%in%GO.l)
plot.df<-sel_GO.df
#6.4*4
library(ggplot2)
ggplot(plot.df,aes(x=type,y=d,color=p,size=c))+geom_point()+
  scale_colour_gradient2(high="blue",
                         mid="#efefef",
                         low="red",
                         midpoint = 0.05,trans="log2")+
  scale_size(range = c(4,10))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.text.y = element_text(size=12, color="black"))+
  labs(x="Tissues",y="")

make_quantile_trans <- function(x, format = scales::label_number()) {
  name <- paste0("quantiles_of_", deparse1(substitute(x)))
  xs <- sort(x)
  N <- length(xs)
  transform <- function(x) findInterval(x, xs)/N # find the last element that is smaller
  inverse <- function(q) xs[1+floor(q*(N-1))]
  
  scales::trans_new(
    name = name,
    transform = transform,
    inverse = inverse,
    breaks =  function(x, n = 5) inverse(scales::extended_breaks()(transform(x), n)),
    minor_breaks = function(x, n = 5) inverse(scales::regular_minor_breaks()(transform(x), n)),
    format = format,
    domain = xs[c(1, N)]
  )
}

p<-ggplot(subset(plot.df,type%in%c("KO Vs WT","L Vs WT","S Vs WT")),aes(x=type,y=d,color=p,size=c))+geom_point()+
  scale_colour_gradientn(colours=c("#ac4f41","#db524f","#efefef","#4e8bba"),
                         breaks=c(0.000001,0.0001,0.05,1),
                         limits=c(min(plot.df$p),1),trans="log2")+
  scale_size(limits=c(min(plot.df$c),max(plot.df$c)),
             range = c(3,8.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.text.y = element_text(size=12, color="black"))+
  labs(x="Tissues",y="")
ggsave(p,file="GO_3_4.27.pdf",device="pdf",height = 4,width = 5.5)

p<-ggplot(subset(plot.df,type%in%c("L Vs S")),aes(x=type,y=d,color=p,size=c))+geom_point()+
  scale_colour_gradientn(colours=c("#ac4f41","#db524f","#efefef","#4e8bba"),
                                                              breaks=c(0.000001,0.0001,0.05,1),
                                                              limits=c(min(plot.df$p),1),trans="log2")+
  scale_size(limits=c(min(plot.df$c),max(plot.df$c)),
             range = c(3,8.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.text.y = element_blank())+
  labs(x="Tissues",y="")
ggsave(p,file="GO_1_4.27.pdf",device="pdf",height = 4,width = 1.8)



ego <- clusterProfiler::enrichGO(gene          = unique(L_S_ENSG.l),
                                 keyType = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)

go.df<-ego@result
go.df$ratio<-go.df$Count/1784

GO.l<-c("Wnt signaling pathway","tissue migration",
        "negative regulation of growth",
        "regulation of cell growth",
        "actin filament organization",
        "epithelial cell development",
        "cell junction assembly",
        "regulation of cell morphogenesis",
        "regulation of epithelial cell proliferation",
        "mesenchyme development")
plot_go.df<-subset(go.df,Description%in%GO.l)
plot_go.df<-plot_go.df[order(plot_go.df$ratio,decreasing = F),]
plot_go.df$Description<-factor(plot_go.df$Description,levels=plot_go.df$Description)
ggplot(plot_go.df,
       aes(x=ratio,y=Description,color=p.adjust,size=Count))+geom_point()+
  scale_colour_gradientn(colours=c("#ac4f41","#db524f","#efefef","#4e8bba"),
                         breaks=c(0.000001,0.0001,0.05,1),
                         limits=c(0.000001,1),trans="log2")+
  scale_size(range = c(3.5,6))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"), 
        axis.text.y = element_text(size=12, color="black"))+
  labs(x="Tissues",y="")


gene.l<-L_S_deg$M
names(gene.l)<-L_S_ENSG.l
cnetplot(ego,color.params = list(foldChange = gene.l),
         cex.params = list(gene_label = 0.65),
         showCategory = c("Wnt signaling pathway",
                          "cell junction assembly",
                           "regulation of epithelial cell proliferation"))


L_S.m<-na.omit(L_S_results@results[[1]])
#idx<-rowMax(cbind(L_S.m$L_mean,L_S.m$S_mean))>100
#filtered.m<-L_S.m[idx,]

filtered.m<-L_S.m
L_S_gene.df<-data.frame(ENSG=rownames(filtered.m),FC=filtered.m[,3])
ENSG.l<-Reduce(c,lapply(unique(L_S_gene.df$ENSG),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
L_S_gene.df$ENSG<-ENSG.l
gene.l<-L_S_gene.df$FC
names(gene.l)<-L_S_gene.df$ENSG
gene.l<-gene.l[order(gene.l,decreasing=T)]
L_S_gene.df<-L_S_gene.df[order(L_S_gene.df$FC,decreasing=T),]

write.table(L_S_gene.df,"HCT116-L_S.rnk",quote=F,col.names = F,row.names=F,sep="\t")

L.m<-na.omit(L_results@results[[1]])
filtered.m<-L.m
L_gene.df<-data.frame(ENSG=rownames(filtered.m),FC=filtered.m[,3])
ENSG.l<-Reduce(c,lapply(unique(L_gene.df$ENSG),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
L_gene.df$ENSG<-ENSG.l
gene.l<-L_gene.df$FC
names(gene.l)<-L_gene.df$ENSG
gene.l<-gene.l[order(gene.l,decreasing=T)]
L_gene.df<-L_gene.df[order(L_gene.df$FC,decreasing=T),]
write.table(L_S_gene.df,"HCT116-L.rnk",quote=F,col.names = F,row.names=F,sep="\t")

S.m<-na.omit(S_results@results[[1]])
filtered.m<-S.m
S_gene.df<-data.frame(ENSG=rownames(filtered.m),FC=filtered.m[,3])
ENSG.l<-Reduce(c,lapply(unique(S_gene.df$ENSG),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
S_gene.df$ENSG<-ENSG.l
gene.l<-S_gene.df$FC
names(gene.l)<-S_gene.df$ENSG
gene.l<-gene.l[order(gene.l,decreasing=T)]
S_gene.df<-S_gene.df[order(S_gene.df$FC,decreasing=T),]
write.table(S_gene.df,"HCT116-S.rnk",quote=F,col.names = F,row.names=F,sep="\t")


KO.m<-na.omit(KO_results@results[[1]])
filtered.m<-KO.m
KO_gene.df<-data.frame(ENSG=rownames(filtered.m),FC=filtered.m[,3])
ENSG.l<-Reduce(c,lapply(unique(KO_gene.df$ENSG),function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
KO_gene.df$ENSG<-ENSG.l
gene.l<-KO_gene.df$FC
names(gene.l)<-KO_gene.df$ENSG
gene.l<-gene.l[order(gene.l,decreasing=T)]
KO_gene.df<-KO_gene.df[order(KO_gene.df$FC,decreasing=T),]
write.table(KO_gene.df,"HCT116-KO.rnk",quote=F,col.names = F,row.names=F,sep="\t")


em2 <- enricher(names(gene.l), TERM2GENE = m_t2g)

clusterProfiler::dotplot(em2, showCategory=10)

cnetplot(em2,showCategory = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                              "HALLMARK_TGF_BETA_SIGNALING"))


library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, human_ensembl_gene)

sel_t2g <- subset(m_t2g,gs_name=="WU_CELL_MIGRATION")
em2 <- GSEA(gene.l, TERM2GENE = sel_t2g,pvalueCutoff = 1)

enrichplot::gseaplot2(em2, geneSetID = 1, title = em2$Description[1])

anno <- em2[1, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

p1 <- enrichplot::gseaplot2(em2, geneSetID = 1,
               title = em2$Description[1])+
  annotate("text", 10000, em2[1, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)


gene.l<-read.delim("gene.list",header=F,sep=" ")
ENSG.l<-Reduce(c,lapply(gene.l$V1,function(x){
  strsplit(x,split=".",fixed = T)[[1]][1]
}))
ego <- clusterProfiler::enrichGO(gene          = unique(ENSG.l),
                                 keyType = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
tcga_result.m<-ego@result
heatplot(ego, showCategory = c("actin filament organization",
                               "cell junction assembly",
                               "cell-substrate adhesion"))
fc.l<-gene.l$V2
names(fc.l)<-ENSG.l
cnetplot(ego,color.params = list(foldChange = fc.l),
         showCategory = c("actin filament organization",
                              "cell junction assembly",
                              "cell-substrate adhesion"))