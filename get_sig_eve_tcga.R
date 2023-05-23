args = commandArgs(trailingOnly=TRUE)

normal_fn<-paste("/picb/rnasys2/huyue/tcga/filtered_psi/",args[1],"_psi_normal.txt",
               sep="",collapse = "")
print(normal_fn)
normal_psi.m<-read.delim(normal_fn,sep="\t",row.names = 1,header=F)
tcga_fn<-paste("/picb/rnasys2/huyue/tcga/filtered_psi/",args[1],"_psi.txt",
               sep="",collapse = "")
print(tcga_fn)
tcga_psi.m<-read.delim(tcga_fn,sep="\t",row.names = 1,header=F)

sig_even<-Reduce(rbind,lapply(rownames(tcga_psi.m),function(name){
  normal.idx<-which(rownames(normal_psi.m)==name)
  tcga.idx<-which(rownames(tcga_psi.m)==name)
  normal_psi<-na.omit(as.numeric(normal_psi.m[normal.idx,]))
  tcga_psi<-na.omit(as.numeric(tcga_psi.m[tcga.idx,]))
  total_psi<-c(normal_psi,tcga_psi)
  a.idx<-which(normal_psi!=1)
  b.idx<-which(normal_psi!=0)
  c.idx<-which(tcga_psi!=1)
  d.idx<-which(tcga_psi!=0)
  if(length(intersect(a.idx,b.idx))>5&&
     length(intersect(c.idx,d.idx))>10){
    re<-wilcox.test(normal_psi,tcga_psi)
    if(re$p.value<0.05){
      return(c(name,re$p.value,mean(tcga_psi)-mean(normal_psi)))
    }
  }
}))

new_p<-p.adjust(sig_even[,2])
new_sig_even<-sig_even[new_p<0.05,]

re.df<-data.frame(eve=sig_even[,1],pval=sig_even[,2],
                  adj_p=new_p,delta_psi=sig_even[,3],
                  type=rep(as.character(args[1]),nrow(sig_even)))
outn<-paste0(args[1],"_sig_eve_tcga.txt")
write.table(re.df,outn,quote=F,row.names = F,col.names = F)