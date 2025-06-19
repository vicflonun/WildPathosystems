setwd("~/Documents/zymoseptoria.microbiome/interactions/paper")
plot.dir=paste(getwd(),"plots",sep = "/")
dir.create(plot.dir)

library(pheatmap)
library(colorspace)
library(ggplot2)

tema=theme(axis.text.x = element_text(color="black",size=8, angle=0,hjust=0.5,vjust=0.5,family = "Arial"),
           axis.text.y = element_text(color="black",size=8,family = "Arial"),
           axis.title = element_text(color="black",size=8,family = "Arial"),
           legend.text = element_text(color = "black",size=8,family = "Arial"),
           legend.key.size = unit(0.3,"cm"),
           legend.title = element_text(color="black",size=8,family = "Arial"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=8, color="black",family = "Arial", margin = margin(0,0,0,0,"cm")),
           strip.text.y = element_text(size=8, color="black",family = "Arial"),
           strip.background = element_rect(fill="white"),
           panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")

#colors
colores=c("#A74476","#EC6C99","#804795","#DF87FF","#524CA1","#9A94EF",
           "#0C2969","#0a468d","#2C79AE","#6EC4FF", "#369DA4","#6EF6FF",
           "#256033","#29A876","#43FFB5","#55A054","#A5F77A","#707032", "#A5A53F",
           "goldenrod","#E3E252","#E16B43","#EC935B","#FFBB7A","#A12F2F","#FF4848")

#Host
h=c("Aegilops.cylindrica","Hordeum.murinum")[2]
#Analysis when Zymoseptoria is in the "t"op or in the "b"ottom
#c=compatible, n=incompatible
z=c("t","b")[2] 

#Load data 
if (h=="Aegilops.cylindrica"){
int=list(
  tc=read.delim(file="all.zt469.top"),tn=read.delim(file="all.zt549.top"),
  bc=read.delim(file="all.zt469.bottom"),bn=read.delim(file="all.zt549.bottom"))
  lim=c("mock","zt469","zt549")
} else if (h=="Hordeum.murinum"){
  int=list(
  tc=read.delim(file="all.zpa796.top"),tn=read.delim(file="all.zpa21.top"),
  bc=read.delim(file="all.zpa796.bottom"),bn=read.delim(file="all.zpa21.bottom"))
  lim=c("mock","zpa796","zpa21")
}

st=read.delim(file="~/Documents/zymoseptoria.microbiome/culture collection/Paper/bacteria_collection_wild_grasses_metadata_v3")
st.h=st[st$host==h,]

colores2=colorRampPalette(colores)
colores2=colores2(44)

names(colores2)=unique(sort(st[,"genus"]))

#Making nomenclature uniform between Malien and my data. 
for (i in names(int)){
  data=int[[i]]
  data[is.na(data$rep.1),"rep.1"]=""
  data=data[data$rep.1!="rm",]
  data[grep("in",data$phenotype),"phenotype"]="I"
  data[grep("prom",data$phenotype),"phenotype"]="P"
  data[data$phenotype %in% "0","phenotype"]="N"
  data[is.na(data$phenotype),"phenotype"]="N"  #Top colony did not grew
  data[data$phenotype %in% "no effect","phenotype"]="N"
  data[data$top.growth %in% "Y","top.growth"]=1
  data[data$top.growth %in% "N","top.growth"]=0
  data[data$cells.under %in% "growth","cells.under"]=1
  data[data$cells.under %in% "no growth","cells.under"]=0
  data[data$cells.under %in% "no data","cells.under"]="ND"
  data[data$intensity %in% "weak","intensity"]="W"
  data[data$intensity %in% "strong","intensity"]="S"
  data[is.na(data$halo.cm),"halo.cm"]=0
  data[is.na(data$colony.cm),"colony.cm"]=0
  data$ratio=data$halo.cm/data$colony.cm
  data[is.na(data$ratio),"ratio"]=0
  data[data$ratio==1,"intensity"]="B"
  data[data$ratio==1,"phenotype"]="N"  #Not considering inhibition only under
  int[[i]]=data
}

#does rows match?
dc=int[[paste(z,"c",sep = "")]] #extract data from list
dn=int[[paste(z,"n",sep = "")]]

#Sort
dc=dc[order(dc$strain.id),]
dn=dn[order(dn$strain.id),]
sum(dn$strain.id==dc$strain.id)==dim(dc)[1] #matches

#create matrix
d=data.frame(strain.id=dc$strain.id, compatible=dc$phenotype, incompatible=dn$phenotype)
d=unique(d) #merges technical replicates 

#For some strains 1 of 3 replicates was different, change it to the common one. 
if (z=="b"){
  d[d$strain.id=="ac37" & d$compatible=="N","compatible"]="I"
  #d[d$strain.id=="hm25" & d$incompatible=="N","incompatible"]="P"
  d=unique(d)
  }
if (z=="t"){
  d[d$strain.id=="hm25" & d$incompatible=="N","incompatible"]="P"
  d=unique(d)
}

rownames(d)=d$strain.id ; d=d[,-1]

#filter strain list, for 1 or 2 strains there is no confrontation data. 
st.h[!st.h$id %in% rownames(d),]
st.h=st.h[st.h$id %in% rownames(d),]
d=d[rownames(d) %in% st.h$id,]
rownames(st.h)=st.h$id

d.st=cbind(d, st.h[rownames(d),c("treatment","genus","phylotype")]) #bind confrontation data and taxonomy

count=data.frame() #Create frequencies
for (i in unique(d.st$treatment)){
  for (j in unique(st.h$genus)){
    
  sub=d.st[d.st$treatment==i & d.st$genus==j,]
  c=summary(factor(sub$compatible,levels=c("I","N","P")))
  plot.dir=paste(getwd(),"plots",sep = "/")
  n=summary(factor(sub$incompatible, levels=c("I","N","P")))
  
  count=rbind(count,data.frame(count=c,freq=c/dim(d.st[d.st$treatment==i,])[1]*100, 
                               zt="compatible", phenotype=c("I","N","P"), treatment=i, genus=j))
  count=rbind(count,data.frame(count=n, freq=c/dim(d.st[d.st$treatment==i,])[1]*100, 
                               zt="incompatible", phenotype=c("I","N","P"), treatment=i, genus=j))
  }
  

}

#plots by treatment
if (z=="b"){nam="Bacteria_vs_Zymo.png"}else if(z=="t"){nam="Zymo_vs_Bacteria.png"}
col= c("#5F1415","grey50","#004B40")

count$treatment=factor(count$treatment, levels = lim)
plot=ggplot(count,aes(x=zt, y=freq, fill=phenotype))+
  geom_bar(stat = "identity")+
  facet_wrap(~treatment, )+
  scale_fill_manual(values=col)+
  scale_x_discrete(label=c("V","A"))+tema+
  ylab("Fequency (%)")+xlab("Tested strain")+ggtitle(nam,h)
plot

png(paste(plot.dir,paste("freq_treatments",h,nam,sep = "."),sep = "/"), 
    width = 1200, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()

#plots by treatment and together all the treatments 
all=count
all$treatment="All"
all$freq=all$count/dim(d.st)[1]*100
all=rbind(all,count)
all$treatment=factor(all$treatment, levels = c("All",lim))
plot=ggplot(all,aes(x=zt, y=freq, fill=phenotype))+
  geom_bar(stat = "identity")+
  facet_wrap(~treatment,nrow = 1)+
  scale_fill_manual(values=col)+
  scale_x_discrete(label=c("V","A"))+tema+
  ylab("Fequency (%)")+xlab("Tested strain")+ggtitle(nam,h)
plot

png(paste(plot.dir,paste("freq_treatments_all",h,nam,sep = "."),sep = "/"), 
    width = 1400, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()

#The same plot but reordered
plot=ggplot(count,aes(x=treatment, y=freq, fill=phenotype))+
  geom_bar(stat = "identity")+
  facet_wrap(~zt)+
  scale_fill_manual(values=col)+tema+
  ylab("Fequency (%)")+xlab("Tested strain")+ggtitle(nam,h)
plot
png(paste(plot.dir,paste(h,"freq_treatments_flip",nam,sep = "."),sep = "/"), 
    width = 1200, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()

#The same plot buth highligting the origin of the strains
plot=ggplot(count,aes(x=zt, y=count, fill=treatment))+
  geom_bar(stat = "identity")+
  facet_wrap(~phenotype)+tema+
  scale_fill_manual(values=c("#54C073","#B80C93","#EF4035"))+
  scale_x_discrete(label=c("V","A"))+
  ylab("No. of strains")+xlab("Tested strain")+ggtitle(nam,h)
plot
png(paste(plot.dir,paste(h,"count_by_treatment",nam,sep = "."),sep = "/"), 
    width = 1200, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()

#Frequency plot by bacterial taxa
#count$genus=factor(count$genus, levels = uniques)
plot=ggplot(count,aes(x=zt, y=count, fill=genus))+
  geom_bar(stat = "identity")+
  facet_wrap(~phenotype)+tema+
  scale_fill_manual(values=colores2[unique(count$genus)])+
  scale_x_discrete(label=c("V","A"))+
  ylab("No. of strains")+xlab("Tested strain")+ggtitle(nam,h)
plot

png(paste(plot.dir,paste(h,"freq4",nam,sep = "."),sep = "/"), 
    width = 1900, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()

#Frequency plot by bacterial taxa in comparison with all the isolates
holi=summary(factor(st.h$genus), maxsum = length(st.h$genus))
single=data.frame(count=holi, freq=holi/dim(st.h)[1]*100,zt="all",phenotype="All", treatment="all",genus=names(holi))
full=rbind(count,single)
full$freq2=full$count/dim(d.st)[1]*100

plot=ggplot(full,aes(x=zt, y=count, fill=genus))+
  geom_bar(stat = "identity")+tema+
  facet_wrap(~phenotype,scales = "free_x", nrow = 1)+
  scale_fill_manual(values=colores2[(unique(full$genus))])+
  scale_x_discrete(label=c("V","A"))+tema+
  ylab("No. of strains")+xlab("Tested strain")+ggtitle(nam,h)
plot

png(paste(plot.dir,paste(h,"freq5",nam,sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()

#The same plot but in proportion of strains
plot=ggplot(full,aes(x=zt, y=freq2, fill=genus))+
  geom_bar(stat = "identity")+
  facet_wrap(~phenotype,scales = "free_x", nrow = 1)+
  scale_fill_manual(values=colores2[(unique(full$genus))])+
  scale_x_discrete(label=c("V","A"))+tema+
  ylab("% of strains")+xlab("Tested strain")+ggtitle(nam,h)
plot

png(paste(plot.dir,paste(h,"freq6",nam,sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
plot ; dev.off()


#### FISHER TEST between Zymoseptoria lineages
detach("package:Rmisc", unload = TRUE)
detach("package:plyr", unload = TRUE)
library(dplyr)
library(tidyr)

count$treatment=factor(count$treatment, levels = lim)
fish= count %>% group_by(zt, phenotype) %>% summarise(cuenta=sum(count))
fish.t=as.data.frame(spread(fish, key=phenotype, value=cuenta))
rownames(fish.t)=fish.t$zt ; fish.t=fish.t[,-1]

fisher.test(fish.t)

#pairwise_fisher_test(as.matrix (data), p.adjust.method = "fdr")

mosaicplot(fish.t, color = TRUE)  

































