
tcga.l<-read.delim("sel_TCGA_normal.list",
                          header=F,sep="\t")

total_eve.m<-Reduce(rbind,apply(tcga.l,1,function(x){
  tcga_fn<-paste(x[1],"_sig_eve_tcga.txt",
                 sep="",collapse = "")
  print(tcga_fn)
  tcga_eve.m<-read.delim(tcga_fn,sep=" ",header=F)
  re.m<-tcga_eve.m[abs(tcga_eve.m[,4])>0.05,]
  return(re.m)
  }))

eve_name.m<-read.delim("ideas_ase_info.tsv",header=T)

library(dplyr)
eve_count<-total_eve.m%>%group_by(V1)%>%summarise(count=n(),pos=sum(V4>0),
                                                  neg=sum(V4<0))

eve_count<-total_eve.m%>%filter(V3<0.05)%>%group_by(V1)%>%summarise(count=n(),pos=sum(V4>0),
                                                  neg=sum(V4<0))

eve_count=eve_count[order(eve_count$count, decreasing = T),]

sel_eve<-subset(eve_count,count>=10)$V1
plot.df<-subset(total_eve.m,V1%in%sel_eve)

merge.df<-merge(plot.df,eve_name.m,by.x="V1",by.y="ase_id")
out.df<-merge.df%>%group_by(gene_id)%>%summarise(psi=median(V4))
write.table(out.df,"gene.list",quote=F,row.names = F,col.names = F)
#merge.df[merge.df$V4>0.2,]$V4<-0.2
#merge.df[merge.df$V4 < -0.2,]$V4<- -0.2

counts_i_eve=data.frame(table(eve_count$count))

counts_i_eve$fraction=(counts_i_eve$Freq)/sum(counts_i_eve$Freq)

#4.5*3
ggplot(data=counts_i_eve, aes(x=Var1,y=fraction))+
  geom_bar(stat = "identity", color="black", fill="grey80", size=.2)+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=14),
        axis.title = element_text(color="black", size=18))+
  ylab("Fraction\nof Sig. Events")+
  xlab("No. of Cancer Types\nShowing Significance")

total_psi.m<-Reduce(rbind,apply(tcga.l,1,function(x){
  tcga_fn<-paste("SLK_psi/",x[1],"_psi.txt",
                 sep="",collapse = "")
  print(tcga_fn)
  tcga_psi.m<-read.delim(tcga_fn,sep="\t",header=F)
  tcga_psi<-na.omit(as.numeric(tcga_psi.m[1,]))
  tcga.df<-data.frame(psi=tcga_psi,
                      class=rep(x[1],length(tcga_psi)),
                      type=rep("cancer",length(tcga_psi)))
  normal_fn<-paste("SLK_psi/",x[1],"_psi_normal.txt",
                 sep="",collapse = "")
  print(normal_fn)
  normal_psi.m<-read.delim(normal_fn,sep="\t",header=F)
  normal_psi<-na.omit(as.numeric(normal_psi.m[1,]))
  t<-paste0(x[1],"-normal")
  normal.df<-data.frame(psi=normal_psi,
                        class=rep(x[1],length(normal_psi)),
                        type=rep("normal",length(normal_psi)))
  re.df<-rbind(tcga.df,normal.df)
  return(re.df)
}))
#sava as SLK_barplot.Rdata

total_count.m<-Reduce(rbind,apply(tcga.l,1,function(x){
  tcga_fn<-paste("RBFOX2_exp/",x[1],"_count.txt",
                 sep="",collapse = "")
  print(tcga_fn)
  tcga_count.m<-read.delim(tcga_fn,sep="\t",header=F)
  tcga_count<-na.omit(log2(as.numeric(tcga_count.m[1,])))
  tcga.df<-data.frame(count=tcga_count,
                      class=rep(x[1],length(tcga_count)),
                      type=rep("cancer",length(tcga_count)))
  normal_fn<-paste("RBFOX2_exp/",x[1],"_count_normal.txt",
                   sep="",collapse = "")
  print(normal_fn)
  normal_count.m<-read.delim(normal_fn,sep="\t",header=F)
  normal_count<-na.omit(log2(as.numeric(normal_count.m[1,])))
  t<-paste0(x[1],"-normal")
  normal.df<-data.frame(count=normal_count,
                        class=rep(x[1],length(normal_count)),
                        type=rep("normal",length(normal_count)))
  re.df<-rbind(tcga.df,normal.df)
  return(re.df)
}))
#sava as RBFOX2_barplot.Rdata