##Generates all plots for the confrontation data of Zymoseptoria lineages and fungi microbiome

setwd("~/Documents/zymoseptoria.microbiome/interactions/Fungi")

library(Biostrings)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
der=colorRampPalette(c("#004B40","#009E8E" ,"#89D9CF","white"))
izq=colorRampPalette(c("white","#F7D3D3","#CB6F70","#5F1415"))
colores2=c("#A74476","#EC6C99","#804795","#DF87FF","#524CA1","#9A94EF",
          "#0C2969","#0a468d","#2C79AE","#6EC4FF", "#369DA4","#6EF6FF",
          "#256033","#29A876","#43FFB5","#55A054","#A5F77A","#707032", "#A5A53F",
          "goldenrod","#E3E252","#E16B43","#EC935B","#FFBB7A","#A12F2F","#FF4848")

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

#Load data
data=read.delim("fungi.int.corrected.txt")
control=read.delim("fungi.int.control.txt")
strains=read.delim("fungi_collection_wild_grasses_metadata_v3")

#Summarizes the mean colony sizes and merges the data for filamentous fungi
fc = data %>% group_by(strain) %>% summarise(diam = mean(fungal.control,na.rm = T))
fcom = data %>% group_by(strain) %>% summarise(diam = mean(fungal.test.compatible,na.rm = T))
finc = data %>% group_by(strain) %>% summarise(diam = mean(fungal.test.incompatible,na.rm = T))

print=cbind(fc, compatible_diameter=fcom$diam, incompatible_diameter=finc$diam)
colnames(print)[1:2]=c("Strain","control_diameter")

met=read.delim(file = "fungi_collection_wild_grasses_metadata_v3")

#Filters data for fungi in confrontation but not in culture collection
met[!met$id %in% print$Strain ,]
print[!print$Strain %in% met$id ,]
print2=print[print$Strain %in% met$id ,]
met[!met$id %in% print2$Strain ,]

#####
#Remove strains not identified
print=print[print$Strain %in% met$id ,]
write.table(print,"summary_fungi_interaction.txt",sep = "\t",quote = F, row.names = F)

#Orders the data as in the metadata
fc=fc[fc$strain %in% met$id ,]
fcom=fcom[fcom$strain %in% met$id ,]
finc=finc[finc$strain %in% met$id ,]

#Calculates the ratios
fc$com.ratio=fcom$diam/fc$diam
fc$inc.ratio=finc$diam/fc$diam

fc$com.fc=log2(fcom$diam/fc$diam)
fc$inc.fc=log2(finc$diam/fc$diam)

#Removes NAs
fc[is.na(fc$inc.fc),"inc.fc"]=0

#Defines interaction based on the threshold
fc$f.phe.c=NA
fc[fc$com.fc==0,"f.phe.c"]="N"
fc[fc$com.fc < 0.32 & fc$com.fc > -0.32 ,"f.phe.c"]="N"
fc[fc$com.fc > 0.32  ,"f.phe.c"]="P"
fc[fc$com.fc < -0.32  ,"f.phe.c"]="I"

fc$f.phe.i=NA
fc[fc$inc.fc==0,"f.phe.i"]="N"
fc[fc$inc.fc < 0.32 & fc$inc.fc > -0.32 ,"f.phe.i"]="N"
fc[fc$inc.fc > 0.32  ,"f.phe.i"]="P"
fc[fc$inc.fc < -0.32  ,"f.phe.i"]="I"

#Summarizes the mean colony sizes and merges the data for Zymoseptoria
zcom = data %>% group_by(strain) %>% summarise(diam = mean(zymo.compatible,na.rm = T))
zinc = data %>% group_by(strain) %>% summarise(diam = mean(zymo.incompatible,na.rm = T))
zcom=zcom[zcom$strain %in% met$id ,]
zinc=zinc[zinc$strain %in% met$id ,]

print=cbind(fc[,1], compatible_diameter=zcom$diam, incompatible_diameter=zinc$diam)
#colnames(print)[1:2]=c("Strain")
#write.table(print,"summary_zymoseptoria_interaction.txt",sep = "\t",quote = F, row.names = F)

#Calculates the ratios
zcom[grep("ac",zcom$strain),"ratio.com"]=zcom[grep("ac",zcom$strain),"diam"]/mean(control[control$species=="zt","compatible"])
zcom[grep("hm",zcom$strain),"ratio.com"]=zcom[grep("hm",zcom$strain),"diam"]/mean(control[control$species=="zpa","compatible"])
zcom$fc.com=log2(zcom$ratio.com)

#Calculates the ratios
zinc[grep("ac",zinc$strain),"ratio.inc"]=zinc[grep("ac",zinc$strain),"diam"]/mean(control[control$species=="zt","incompatible"])
zinc[grep("hm",zinc$strain),"ratio.inc"]=zinc[grep("hm",zinc$strain),"diam"]/mean(control[control$species=="zpa","incompatible"])
zinc$fc.inc=log2(zinc$ratio.inc)

#Removes NAs
zinc[is.na(zinc$fc.inc),"fc.inc"]=0
zinc[is.infinite(zinc$fc.inc),"fc.inc"]=-2 #pseudocount, for the zymo strains that did not grew in presence of fungi

zcom[is.na(zcom$fc.com),"fc.com"]=0
zcom[is.infinite(zcom$fc.com),"fc.com"]=-2 #pseudocount, for the zymo strains that did not grew in presence of fungi

#Defines interaction based on the threshold
zcom$z.phe.c=NA
zcom[zcom$fc.com==0,"z.phe.c"]="N"
zcom[zcom$fc.com < 0.32 & zcom$fc.com > -0.32 ,"z.phe.c"]="N"
zcom[zcom$fc.com > 0.32  ,"z.phe.c"]="P"
zcom[zcom$fc.com < -0.32  ,"z.phe.c"]="I"

zinc$z.phe.i=NA
zinc[zinc$fc.inc==0,"z.phe.i"]="N"
zinc[zinc$fc.inc < 0.32 & zinc$fc.inc > -0.32 ,"z.phe.i"]="N"
zinc[zinc$fc.inc > 0.32  ,"z.phe.i"]="P"
zinc[zinc$fc.inc < -0.32  ,"z.phe.i"]="I"

#Generates data for heatmap
library(pheatmap)
fc=rbind(fc,
           data.frame(strain=c(rep("ac118",1),rep("hm39",1),rep("hm36",1)),
                      diam=NA,com.ratio=NA, inc.ratio=NA, com.fc=0,inc.fc=0 ,f.phe.c="N", f.phe.i="N"))

fc.m=as.matrix(fc[c(5,6)])
rownames(fc.m)=fc$strain
b=c(seq(min(fc.m),max(fc.m),0.05))
#pheatmap(fc.m, breaks = b)


#Yeast isolated did not have any effect on Zymoseptoria
zinc=rbind(zinc,data.frame(strain=c(rep("ac118",1),rep("hm39",1),rep("hm36",1)), diam="NA", ratio.inc=1 , fc.inc =0,  z.phe.i="N"))
zcom=rbind(zcom,data.frame(strain=c(rep("ac118",1),rep("hm39",1),rep("hm36",1)), diam="NA", ratio.com=1 , fc.com =0,  z.phe.c="N"))

z.m=as.matrix(cbind(zcom$fc.com,zinc$fc.inc))
colnames(z.m)=c("com.fc","inc.fc")
rownames(z.m)=fc$strain
b=c(seq(min(z.m),max(z.m),0.05),0)
#pheatmap(z.m, breaks = b)


#Generates the three heatmap
library(ggtree, verbose = F)
library(stringr)
library(ape)

tree=read.tree(file="fungi_confrontation_tree.nwk")
#tree=root(tree, node = 84, edgelabel = T)
ggtree(tree)+
  theme_tree2()+
  geom_tiplab(align=F, linesize=.5) 

#remove strains that are not in the tree
rownames(strains)=strains$id
st=strains[tree$tip.label,]
#fc.m=fc.m[tree$tip.label[1:30],]
z.m=z.m[tree$tip.label,]

st$phylo2=paste(st$id,st$genus,str_extract(st$phylotype,"[0-9]"),sep = "_")
order=tree$tip.label
tree$tip.label=as.character(st[order,"phylo2"])

ggtree(tree)+
  theme_tree2()+
  geom_tiplab(align=F, linesize=.5) 

#found min and mx

min=max(c(fc.m[,1],fc.m[,2],z.m[,1],z.m[,2]))
max=min(c(fc.m[,1],fc.m[,2],z.m[,1],z.m[,2]))

#bL=seq(max,min,0.1)
bL=seq(-2,1,0.1)

#rownames(st)=st$phylo2
colnames(fc.m)[1:2]=c("c_f","i_f")
colnames(z.m)[1:2]=c("c_z","i_z")
#heat=cbind(t.2[,1:2],b.2[,1:2])
meta=cbind(id=rownames(z.m),st[rownames(z.m),c("genus","treatment","time","host","phylo2")])
#meta=meta[!is.na(meta$genus),]
meta[meta$treatment == "zt469" | meta$treatment =="zpa796","treatment"]="compatible"
meta[meta$treatment == "zt549" | meta$treatment =="zpa21","treatment"]="incompatible"
lim=c("mock","compatible","incompatible")
meta$treatment=factor(meta$treatment,levels = lim)
#meta$genus=factor(meta$genus,levels = sort(unique(strains$genus)))

fc.m=fc.m[rownames(meta),]
z.m=z.m[rownames(meta),]

rownames(meta)=meta$phylo2
rownames(fc.m)=meta$phylo2
rownames(z.m)=meta$phylo2

#meta=meta[,c(6,1:5)]

library(ggnewscale)
names(colores2)=sort(unique(meta$genus))

meta=meta[,c(6,1:5)]

#tree$tip.label %in% rownames(meta)
p=ggtree(tree,size=0.5, layout = 'circular') %<+% meta
p=p + geom_tippoint(aes(color=genus),size=4)+
  scale_color_manual(values=colores2[unique(meta$genus)])+
  geom_tiplab(align=T, linesize=0.05, size=3,offset = 0.9)
p

ge=data.frame(row.names=meta$phylo2,f= meta$genus) 
p2=gheatmap(p, ge, offset=0.75, width=0.03, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_manual(values=colores2[unique(meta$genus)])
p2

p2.2=p2+new_scale_fill()

ge=data.frame(row.names=meta$phylo2,host= meta$host) 
p2.2=gheatmap(p2.2, ge, offset=0.03, width=0.03, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_manual(values=c("#FFC850","#116875"))
p2.2
p3=p2.2+new_scale_fill()

colores.g=c("#54C073","#B80C93","#EF4035")
names(colores.g)=lim

tre=data.frame(row.names=meta$phylo2,treatment= meta$treatment) 
p4=gheatmap(p3, tre, offset=0.1, width=0.03, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_manual(values = colores.g[lim], breaks = lim)
p4

#max=summary(c(heat$c_t,heat$i_t, heat$c_b,heat$i_b))
#bL=seq(max[1],max[6],0.1)

#fc.m[fc.m < 0.32 & fc.m >= 0]=0
#fc.m[fc.m > -0.32 & fc.m <= 0]=0
#z.m[z.m < 0.32 & z.m >= 0]=0
#z.m[z.m > -0.32 & z.m <= 0]=0

p5=p4+new_scale_fill()
p6=gheatmap(p5, fc.m, offset=0.20, width=0.2, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_gradientn(colours = rev(c(der(sum(bL>0)),izq(sum(bL<0)))), limits=c(-2,1))
p6
p7=p6+new_scale_fill()

plot=gheatmap(p7,z.m, offset=0.45, width=0.2, font.size=2, color = NULL,
              colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_gradientn(colours = rev(c(der(sum(bL>0)),izq(sum(bL<0)))), limits=c(-2,1))
plot

png(paste(getwd(),paste("fungi","confrontation3.png",sep = "_"),sep = "/"), 
    width = 6000, height = 6000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()


#Genrates the bloxpot for the quatitaive analysis

library(tidyr)
ge=data.frame(row.names=meta$phylo2,f= meta$genus,h=meta$host) 
plot=cbind(ge,fc.m)
plot$strain=rownames(plot)

#growth enhancement
tema$legend.position="none"
plot2=plot[plot$c_f >= 0.32 | plot$i_f >= 0.32  ,]
plot3=gather(plot2, key = "strain", value = "logFC", c("c_f","i_f"))
a=ggplot(data=plot3, aes(x=strain, y=logFC))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=plot3, aes(x=strain, y=logFC,fill=f),size=1.5,shape=21, width = 0.2)+
  scale_fill_manual(values=colores2[sort(unique(plot2$f))])+tema+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+
  facet_wrap(~h, ncol = 1)+
  geom_hline(yintercept = 0.32,linetype=2)+
  scale_y_continuous(limits = c(-0.1,2))
a
png(paste(getwd(),"fungal_promotion_2.png",sep = "/"), 
    width = 500, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a); dev.off()

wilcox.test(plot3[plot3$h=="Aegilops.cylindrica" & plot3$strain=="c_f","logFC"],
            plot3[plot3$h=="Aegilops.cylindrica" & plot3$strain=="i_f","logFC"])

wilcox.test(plot3[plot3$h=="Hordeum.murinum" & plot3$strain=="c_f","logFC"],
            plot3[plot3$h=="Hordeum.murinum" & plot3$strain=="i_f","logFC"])

#Growth inhibition of fungi
plot2=plot[plot$c_f <= -0.32 | plot$i_f <= -0.32  ,]
plot3=gather(plot2, key = "strain", value = "logFC", c("c_f","i_f"))
plot3$logFC=plot3$logFC*-1
a=ggplot(data=plot3, aes(x=strain, y=logFC))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=plot3, aes(x=strain, y=logFC,fill=f),size=1.5,shape=21, width = 0.2)+
  scale_fill_manual(values=colores2[sort(unique(plot2$f))])+tema+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+
  facet_wrap(~h, ncol = 1)+
  geom_hline(yintercept = 0.32,linetype=2)+
  scale_y_continuous(limits = c(-0.1,2))
a
png(paste(getwd(),"fungal_inhibition2.png",sep = "/"), 
    width = 500, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a); dev.off()

wilcox.test(plot3[plot3$h=="Aegilops.cylindrica" & plot3$strain=="c_f","logFC"],
            plot3[plot3$h=="Aegilops.cylindrica" & plot3$strain=="i_f","logFC"])

wilcox.test(plot3[plot3$h=="Hordeum.murinum" & plot3$strain=="c_f","logFC"],
            plot3[plot3$h=="Hordeum.murinum" & plot3$strain=="i_f","logFC"])


#Growth inhibition of Zymoseptoria
ge=data.frame(row.names=meta$phylo2,f= meta$genus,h=meta$host) 
plot=cbind(ge,z.m)
plot$strain=rownames(plot)

plot2=plot[plot$c_z <= -0.32 | plot$i_z <= -0.32  ,]
plot3=gather(plot2, key = "strain", value = "logFC", c("c_z","i_z"))
plot3$logFC=plot3$logFC*-1
a=ggplot(data=plot3, aes(x=strain, y=logFC))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=plot3, aes(x=strain, y=logFC,fill=f),size=1.5,shape=21, width = 0.2)+
  scale_fill_manual(values=colores2[sort(unique(plot2$f))])+tema+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+
  facet_wrap(~h, ncol = 1)+
  geom_hline(yintercept = 0.32,linetype=2)+
  scale_y_continuous(limits = c(-0.1,2.2))
a
png(paste(getwd(),"zymo_inhibition2.png",sep = "/"), 
    width = 500, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a); dev.off()


shapiro.test(plot3$logFC)
wilcox.test(plot3[plot3$h=="Aegilops.cylindrica" & plot3$strain=="c_z","logFC"],
           plot3[plot3$h=="Aegilops.cylindrica" & plot3$strain=="i_z","logFC"])

wilcox.test(plot3[plot3$h=="Hordeum.murinum" & plot3$strain=="c_z","logFC"],
            plot3[plot3$h=="Hordeum.murinum" & plot3$strain=="i_z","logFC"])


#Generates collot legend for all plots
holi=strains
holi$point=1
tema$legend.position="right"
a=ggplot(holi, aes(x=host, y=point, fill=genus))+
  geom_point(shape=21)+scale_fill_manual(values = colores2[sort(unique(holi$genus))])+tema
a
png(paste(getwd(),"colores.png",sep = "/"), 
    width = 1300, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a); dev.off()


#Generates the frequency plots
#summarizes data for filamentous fungi
lev=c("I","N","P")
fungal=as.data.frame(cbind(
summary(factor(fc$f.phe.c[grep("ac",fc$strain)],levels = lev))/26*100,
summary(factor(fc$f.phe.i[grep("ac",fc$strain)],levels = lev))/26*100,
summary(factor(fc$f.phe.c[grep("hm",fc$strain)],levels = lev))/16*100,
summary(factor(fc$f.phe.i[grep("hm",fc$strain)],levels = lev))/16*100
))

colnames(fungal)=c("ac.C","ac.I","hm.C","hm.I")
fungal$phenotype=rownames(fungal)

library(tidyr)
fungal=gather(fungal, key="test", value="count", c("ac.C","ac.I","hm.C","hm.I"))
fungal$host=c(rep("ac",6),rep("hm",6))
fungal$strain=c(rep("C",3),rep("I",3),rep("C",3),rep("I",3))

#Summarizes the data for Zymoseptoria
zymo=data.frame(cbind(
summary(factor(zcom$z.phe.c[grep("ac",fc$strain)],levels = lev))/29*100,
summary(factor(zinc$z.phe.i[grep("ac",fc$strain)],levels = lev))/29*100,

summary(factor(zcom$z.phe.c[grep("hm",fc$strain)],levels = lev))/17*100,
summary(factor(zinc$z.phe.i[grep("hm",fc$strain)],levels = lev))/17*100
))
colnames(zymo)=c("ac.C","ac.I","hm.C","hm.I")
zymo$phenotype=rownames(zymo)

#Creates the frequency plots
library(tidyr)
col= c(izq(5)[c(5)],"grey50",der(5)[1])

zymo=gather(zymo, key="test", value="count", c("ac.C","ac.I","hm.C","hm.I"))
zymo$host=c(rep("ac",6),rep("hm",6))
zymo$strain=c(rep("C",3),rep("I",3),rep("C",3),rep("I",3))

plot=ggplot(fungal,aes(x=strain, y=count, fill=phenotype))+
  geom_bar(stat = "identity")+
  facet_wrap(~host,scale="free_y"  )+
  scale_fill_manual(values=col)+
  scale_x_discrete(label=c("V","A"))+tema+
  ylab("Fequency (%)")+xlab("Tested strain")+ggtitle("Efect on Fungi")
plot

png("Effect.fungi_2.png", 
    width = 1000, height = 700, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
print(plot) ; dev.off()


plot=ggplot(zymo,aes(x=strain, y=count, fill=phenotype))+
  geom_bar(stat = "identity")+
  facet_wrap(~host, scale="free_y" )+
  scale_fill_manual(values=col)+
  scale_x_discrete(label=c("V","A"))+tema+
  ylab("Fequency (%)")+xlab("Tested strain")+ggtitle("Efect on Zymoseptoria")
plot

png("Effect.zymo_2.png", 
    width = 1000, height = 700, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
print(plot) ; dev.off()

col= c(izq(5)[c(5)],"grey50",der(5)[1])


#Fisher test. 

lev=c("I","N","P")

#For filamentous fungi
fungal=as.data.frame(cbind(
  summary(factor(fc$f.phe.c[grep("ac",fc$strain)],levels = lev)),
  summary(factor(fc$f.phe.i[grep("ac",fc$strain)],levels = lev)),
  summary(factor(fc$f.phe.c[grep("hm",fc$strain)],levels = lev)),
  summary(factor(fc$f.phe.i[grep("hm",fc$strain)],levels = lev))
))

colnames(fungal)=c("ac.C","ac.I","hm.C","hm.I")
#fungal$phenotype=rownames(fungal)

fu.a=t(fungal)[1:2,]
fu.h=t(fungal)[3:4,]

fisher.test(fu.a)
fisher.test(fu.h)

#For Zymoseptoria
zymo=data.frame(cbind(
  summary(factor(zcom$z.phe.c[grep("ac",fc$strain)],levels = lev)),
  summary(factor(zinc$z.phe.i[grep("ac",fc$strain)],levels = lev)),
  
  summary(factor(zcom$z.phe.c[grep("hm",fc$strain)],levels = lev)),
  summary(factor(zinc$z.phe.i[grep("hm",fc$strain)],levels = lev))
))

colnames(zymo)=c("ac.C","ac.I","hm.C","hm.I")
#zymo$phenotype=rownames(zymo)

zy.a=t(zymo)[1:2,]
zy.h=t(zymo)[3:4,]

fisher.test(zy.a)
fisher.test(zy.h)

mosaicplot(t(zymo),color = T)

#Write summary data: 

write.table(fc,"Effect_of_Zymoseptoria_on_Fungi.txt",sep = "\t",quote = F, row.names = F)

zc=cbind(zcom,zinc[,3:5])


write.table(zc,"Effect_of_Fungi_on_Zymoseptoria.txt",sep = "\t",quote = F, row.names = F)



