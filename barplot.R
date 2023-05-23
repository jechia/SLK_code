library(ggplot)
library(ggpubr)

load("D:/2021/SLK/SLK_barplot.RData")
total_psi.m$type<-factor(total_psi.m$type,levels=c("normal","cancer"))
ggplot(total_psi.m,aes(x=class,y=psi,fill=type))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey", "#db524f"))+
  stat_compare_means(aes(group = type),
                     method = "wilcox.test",
                     label="p.signif")+
  theme_classic()


load("D:/2021/SLK/RBFOX2_count.RData")
total_count.m$type<-factor(total_count.m$type,levels=c("normal","cancer"))
ggplot(total_count.m,aes(x=class,y=count,fill=type))+
  geom_boxplot()+
  scale_fill_manual(values = c("grey", "#db524f"))+
  stat_compare_means(aes(group = type),
                     method = "wilcox.test",
                     label="p.signif")+
  theme_classic()