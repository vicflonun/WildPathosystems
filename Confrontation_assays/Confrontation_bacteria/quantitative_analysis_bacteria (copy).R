#Generates boxplots and heatmaps for the confrontation data between Zymoseptoria isolates and the bacterial microbiome

#Directories
setwd("~/Documents/zymoseptoria.microbiome/interactions/paper")
plot.dir=paste(getwd(),"plots_test",sep = "/")
dir.create(plot.dir)

#setwd("E:/interactions")
library(pheatmap)
library(colorspace)
library(ggrepel)
library(ggplot2)
library(Rmisc)
library(tidyr)

#plot settings
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

der=colorRampPalette(c("#004B40","#009E8E" ,"#89D9CF","white"))
izq=colorRampPalette(c("white","#F7D3D3","#CB6F70","#5F1415"))

colores1=c("#54C073","#B80C93","#EF4035")

#Select for host 
h=c("Aegilops.cylindrica","Hordeum.murinum")[2]
#Select for analysis when Zymoseptoria is in the "t"op or in the "b"ottom
z=c("t","b")[1] 

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

st.full=read.delim(file="~/Documents/zymoseptoria.microbiome/culture collection/Paper/bacteria_collection_wild_grasses_metadata_v3")
#st.full=st[st$host==h,]

#Assign colors to the genera and creates a collor palette
colores2=colorRampPalette(colores)
colores2=colores2(44)
names(colores2)=unique(sort(st.full[,"genus"]))

holi=st.full
holi$point=1
a=ggplot(holi, aes(x=host, y=point, fill=genus))+
  geom_point(shape=21)+scale_fill_manual(values = colores2[sort(unique(holi$genus))])+tema

png(paste(plot.dir,"colores.png",sep = "/"), 
    width = 1300, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a); dev.off()

#Concatenates data from experiments registered with different codes.  
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

dc=int[[paste(z,"c",sep = "")]] #extract data from list
dn=int[[paste(z,"n",sep = "")]]

#Sorts the data
dc=dc[order(dc$strain.id),]
dn=dn[order(dn$strain.id),]
sum(dn$strain.id==dc$strain.id)==dim(dc)[1]

#create matrix qualitative
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
st=st.full[st.full$id %in% rownames(d),]
d=d[rownames(d) %in% st$id,]
rownames(st)=st$id

#create matrix quantitative 
dd=data.frame(strain.id=unique(dc$strain.id), 
             compatible=summarySE(dc,"ratio","strain.id")[,3], 
             incompatible=summarySE(dn,"ratio","strain.id")[,3])

rownames(dd)=dd$strain.id ; dd=dd[,-1]

#Extract taxa #matching taxa
st=st.full[st.full$id %in% rownames(dd),]
dd=dd[rownames(dd) %in% st$id,]

#make pheatmap
names(colores1)=lim
col=list(genus=colores2[sort(unique(st$genus))], treatment=colores1)

anr=data.frame(row.names = st$id, treatment=st$treatment, genus=st$genus)
anr=anr[order(anr$treatment),]
ord=rownames(anr[order(anr$genus),])

#TRanform to postive for promotion and negative for inhibition 
ddd=dd
ddd[rownames(d[d$compatible=="I",]),"compatible"]=ddd[rownames(d[d$compatible=="I",]),"compatible"]*-1
ddd[rownames(d[d$incompatible=="I",]),"incompatible"]=ddd[rownames(d[d$incompatible=="I",]),"incompatible"]*-1

if (z=="t"){

ddd.t=ddd
#all  ratios of value 1 are NO effect and will be tranformed to O
ddd.t[ddd.t$compatible==1,"compatible"]=0
ddd.t[ddd.t$incompatible==1,"incompatible"]=0
max=summary(c(ddd.t$compatible,ddd.t$incompatible))
bL = seq(max[1], max[6], by = 0.1)

#dev.off()
png(paste(getwd(),paste("plots",paste(h,"Zymoseptoria_vs_bacteria.png",sep="_"),sep = "/"),sep = "/"), 
    width = 1000, height = 2000, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")

pheatmap(ddd[ord,], cluster_rows = F, cluster_cols = F, breaks = bL,
         #color = rev(c(der(90),izq(50))),
         color = rev(c(der(sum(bL>0)),izq(sum(bL<0)))),
         annotation_row = anr, annotation_colors = col, fontsize = 7, border_color = "black")

dev.off()


} else if (z=="b"){
  
  ddd.b=ddd
  #all  ratios of value 1 are NO effect and will be tranformed to O
  ddd.b[ddd.b$compatible==1,"compatible"]=0
  ddd.b[ddd.b$incompatible==1,"incompatible"]=0
  max=summary(c(ddd.b$compatible,ddd.b$incompatible))
  bL = seq(max[1], max[6], by = 0.1)
  
  #dev.off()
  png(paste(getwd(),paste("plots",paste(h,"bacteria_vs_Zymoseptoria.png",sep="_"),sep = "/"),sep = "/"), 
      width = 1000, height = 2000, units = "px", pointsize = 100, bg = "white", res=300, type="cairo")
  pheatmap(ddd[ord,], cluster_rows = F, cluster_cols = F, breaks = bL,
           color = rev(c(der(sum(bL>0)),izq(sum(bL<0)))),
           annotation_row = anr, annotation_colors = col, fontsize = 7, border_color = "black")
  dev.off()
  
  }

#######

# Generates the quantitative boxplots for inhibitory phenotype
if (z=="b"){nam="Bacteria_vs_Zymo_only_inhibition.png"}else if(z=="t"){nam="Zymo_vs_Bacteria_only_inhibition.png"}

#Subset only inhibition
di=ddd[ddd$compatible<0 | ddd$incompatible<0,]

strains=st[st$id %in% rownames(di),]

## boxplots
di$strain.id=rownames(di)
di=di[strains$id,]
#di$strain.id == strains$id
di$treatment=strains$treatment
di$genus=strains$genus
di$phylotype=strains$phylotype

box=gather(di, "zt","ratio", c("compatible","incompatible"))
box$ratio=box$ratio*-1
box=box[box$ratio>1,] #removes points with promotion or no effect


plot=ggplot(data=box, aes(x=zt, y=ratio, fill=zt))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=box, aes(x=zt, y=ratio,fill=genus),size=1.5,shape=21, width = 0.2)+
  scale_fill_manual(values=col$genus[sort(unique(box$genus))])+
  ylab("inhibition halo/colony ratio")+tema+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+scale_y_continuous(limits = c(1,8))
plot

png(paste(plot.dir,paste(h,"subset_2",nam,sep = "_"),sep = "/"), 
    width = 500, height = 500, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()

wilcox.test(ratio~zt, data=box)

box$treatment=factor(box$treatment, levels = lim)
plot=ggplot(data=box, aes(x=zt, y=ratio, fill=zt))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=box, aes(x=zt, y=ratio,fill=genus),size=2,shape=21)+
  scale_fill_manual(values=col$genus[sort(unique(box$genus))])+
  facet_wrap(~treatment)+  ylab("inhibition halo/colony ratio")+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+scale_y_continuous(limits = c(1,8))

plot
png(paste(getwd(),"plots",paste(h,"subset_facet",nam,sep = "_"),sep = "/"), 
    width = 2500, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()

#Test
#shapiro.test(box$ratio)
#wilcox.test(box[box$treatment=="mock" & box$zt=="compatible","ratio"],box[box$treatment=="mock" & box$zt=="incompatible","ratio"])
#wilcox.test(box[box$treatment=="zt469" & box$zt=="compatible","ratio"],box[box$treatment=="zt469" & box$zt=="incompatible","ratio"])
#wilcox.test(box[box$treatment=="zt549" & box$zt=="compatible","ratio"],box[box$treatment=="zt549" & box$zt=="incompatible","ratio"])

library(FSA)
box$join=paste(box$treatment,box$zt,sep = "_")
kruskal.test(ratio ~ join, data = box)
dunnTest(ratio ~ join, data = box, method = "bh")

#Do the same plots for a select group of taxa
genus.dir=paste(plot.dir,"genus",sep = "/") ; dir.create(genus.dir) 

if (z=="b"){
  
nam="Bacteria_vs_Zymo_only_inhibition.png"

if (h=="Hordeum.murinum"){g=c("Pseudomonas","Streptomyces") 
} else { g=c("Microbacterium", "Pseudomonas","Rhodococcus","Paenibacillus","Streptomyces")} 

for (ge in g){

subbox=box[box$genus %in% ge,]
plot=ggplot(data=subbox, aes(x=zt, y=ratio, fill=zt))+
  geom_boxplot(outlier.color=NA, fill="grey50")+
  geom_point(data=subbox, aes(x=zt, y=ratio,fill=genus),size=2,shape=21)+
  geom_line(aes(group=strain.id,color=genus))+
  scale_fill_manual(values=col$genus[sort(unique(subbox$genus))])+
  scale_color_manual(values=col$genus[sort(unique(subbox$genus))])+
  geom_text_repel(label=subbox$phylotype,size=2.5,force_pull = 0.2)+
  facet_wrap(~treatment)+
  ylab("inhibition halo/colony ratio")

png(paste(getwd(),"plots/genus",paste("Box_all_facet",nam,ge,sep = "_"),sep = "/"), 
    width = 2500, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()

#wilcox.test(box[box$treatment=="mock" & box$zt=="compatible","ratio"],box[box$treatment=="mock" & box$zt=="incompatible","ratio"])
#wilcox.test(box[box$treatment=="zt469" & box$zt=="compatible","ratio"],box[box$treatment=="zt469" & box$zt=="incompatible","ratio"])
#wilcox.test(box[box$treatment=="zt549" & box$zt=="compatible","ratio"],box[box$treatment=="zt549" & box$zt=="incompatible","ratio"])

}

}

# Generates the quantitative boxplots for growth enhancement phenotype
if(z=="t"){
  
nam="Zymo_vs_Bacteria_only_promotion.png"
#Subset only inhibition
dp=ddd[ddd$compatible>0 | ddd$incompatible>0,]
strains=st[st$id %in% rownames(dp),]


## boxplots
dp$strain.id=rownames(dp)
dp=dp[strains$id,]
#dp$strain.id == strains$id
dp$treatment=strains$treatment
dp$genus=strains$genus
dp$phylotype=strains$phylotype


box=gather(dp, "zt","ratio", c("compatible","incompatible"))
box=box[box$ratio>1,]

tema$legend.position="none"
plot=ggplot(data=box, aes(x=zt, y=ratio, fill=zt))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=box, aes(x=zt, y=ratio,fill=genus),size=1.5,shape=21, width = 0.2)+
  scale_fill_manual(values=col$genus[sort(unique(box$genus))])+
  ylab("Promotion halo/colony ratio")+tema+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+scale_y_continuous(limits = c(1,8.2))
plot

png(paste(plot.dir,paste(h,"subset_2",nam,sep = "_"),sep = "/"), 
    width = 500, height = 500, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()

shapiro.test(box$ratio)
wilcox.test(box[box$zt=="compatible","ratio"],box[box$zt=="incompatible","ratio"])


box$treatment=factor(box$treatment,levels = lim)
plot=ggplot(data=box, aes(x=zt, y=ratio, fill=zt))+
  geom_boxplot(outlier.colour = "NA",fill="grey50")+
  geom_jitter(data=box, width = 0.2, aes(x=zt, y=ratio,fill=genus),size=2,shape=21)+
  scale_fill_manual(values=col$genus[sort(unique(box$genus))])+
  facet_wrap(~treatment)+tema+
  scale_x_discrete(labels=c("Virulent","Avirulent"))+scale_y_continuous(limits = c(1,8.2))+
  ylab("promotion halo/colony ratio")
plot
png(paste(plot.dir,paste(h,"subset_facet",nam,sep = "_"),sep = "/"), 
    width = 2500, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()


if (h=="Hordeum.murinum"){g=c("Streptomyces","Rhodococcus") 
} else { g=c("Microbacterium","Pseudomonas","Rhodococcus","Curtobacterium","Acidovorax","Streptomyces","Flavobacterium")} 

for (ge in g){
  subbox=box[box$genus %in% ge,]
  plot=ggplot(data=subbox, aes(x=zt, y=ratio, fill=zt))+
    geom_boxplot(outlier.color=NA, fill="grey50")+
    geom_point(data=subbox, aes(x=zt, y=ratio,fill=genus),size=2,shape=21)+
    geom_line(aes(group=strain.id,color=genus))+
    scale_fill_manual(values=col$genus[sort(unique(subbox$genus))])+
    scale_color_manual(values=col$genus[sort(unique(subbox$genus))])+
    geom_text_repel(label=subbox$phylotype,size=2.5,force_pull = 0.2)+
    ylab("promotion halo/colony ratio")
  png(paste(plot.dir,paste("genus/subset",nam,ge,sep = "_"),sep = "/"), 
      width = 1000, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
  print(plot); dev.off()
  


plot=ggplot(data=subbox, aes(x=zt, y=ratio, fill=zt))+
  geom_boxplot(outlier.color=NA, fill="grey50")+
  geom_point(data=subbox, aes(x=zt, y=ratio,fill=genus),size=2,shape=21)+
  geom_line(aes(group=strain.id,color=genus))+
  scale_fill_manual(values=col$genus[sort(unique(subbox$genus))])+
  scale_color_manual(values=col$genus[sort(unique(subbox$genus))])+
  geom_text_repel(label=subbox$phylotype,size=2.5,force_pull = 0.2)+
  facet_wrap(~treatment)+
  ylab("promotion halo/colony ratio")
plot

png(paste(getwd(),"plots/genus",paste(h,"Box_all_facet",nam,ge,sep = "_"),sep = "/"), 
    width = 2500, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()

}

}

#shapiro.test(box$ratio)
#wilcox.test(box[box$treatment=="mock" & box$zt=="compatible","ratio"],box[box$treatment=="mock" & box$zt=="incompatible","ratio"])
#wilcox.test(box[box$treatment=="zt469" & box$zt=="compatible","ratio"],box[box$treatment=="zt469" & box$zt=="incompatible","ratio"])
#wilcox.test(box[box$treatment=="zt549" & box$zt=="compatible","ratio"],box[box$treatment=="zt549" & box$zt=="incompatible","ratio"])

library(FSA)
box$join=paste(box$treatment,box$zt,sep = "_")
kruskal.test(ratio ~ join, data = box)
dunnTest(ratio ~ join, data = box, method = "bh")


##############################################################################

#To create the Three heatmap run the analysis for each host, both Zymoseptoria on top and Zymoseptoria bottom
library(ggtree, verbose = F)
library(stringr)
library(ape)

if (h=="Aegilops.cylindrica"){
  tree=read.tree(file="ac_confrontation_tree.nwk")
  #tree=root(tree, node = 210, edgelabel = T)
}else if (h=="Hordeum.murinum"){
  tree=read.tree(file="hm_confrontation_tree.nwk")
}

#change tip names
rownames(st)=st$id
st$phylo2=paste(st$id,st$genus,str_extract(st$phylotype,"[0-9]"),sep = "_")
order=tree$tip.label
tree$tip.label=as.character(st[order,"phylo2"])

ggtree(tree)+
  theme_tree2()+
  geom_tiplab(align=F, linesize=.5, size=2) 

b.2=ddd.b[rownames(ddd.b) %in% order,]
rownames(b.2)=st[rownames(b.2),"phylo2"]

t.2=ddd.t[rownames(ddd.t) %in% order,]
rownames(t.2)=st[rownames(t.2),"phylo2"]


#all
rownames(st)=st$phylo2
colnames(t.2)[1:2]=c("c_t","i_t")
colnames(b.2)[1:2]=c("c_b","i_b")
heat=cbind(t.2[,1:2],b.2[,1:2])
colnames(heat)=c("b","c","d","e")
meta=cbind(id=rownames(t.2),st[rownames(t.2),c("genus","treatment","time","host","phylo2","phylotype")])
meta$treatment=factor(meta$treatment,levels = lim)
meta$genus=factor(meta$genus,levels = sort(unique(st$genus)))

library(ggnewscale)

p=ggtree(tree,size=0.5, layout = 'circular') %<+% meta
p=p + geom_tippoint(aes(color=genus), size=2.5)+
  scale_color_manual(values=colores2[unique(st$genus)])+
  geom_tiplab(align=T, linesize=0.05, size=2,offset = 0.6)
p

ge=data.frame(row.names=meta$id,f= meta$genus) 
p2=gheatmap(p, ge, offset=0.55, width=0.03, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_manual(values=colores2[unique(st$genus)])
p2
#p2=p
p3=p2+new_scale_fill()


names(colores1)=lim

tre=data.frame(row.names=meta$id,a= meta$treatment) 
p4=gheatmap(p3, tre[], offset=0.05, width=0.03, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_manual(values = colores1[lim], breaks = lim)
p4

max=summary(c(heat$a,heat$b, heat$c,heat$d))
bL=seq(max[1],max[6],0.1)

p5=p4+new_scale_fill()
p6=gheatmap(p5, heat[,1:2], offset=0.15, width=0.15, font.size=2, color = NULL,
            colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_gradientn(colours = rev(c(der(sum(bL>0)),izq(sum(bL<0)))), limits=c(max[1],max[6]))

p7=p6+new_scale_fill()

plot=gheatmap(p7, heat[,3:4], offset=0.35, width=0.15, font.size=2, color = NULL,
              colnames_angle=0, hjust=0.5, colnames_position="bottom")+
  scale_fill_gradientn(colours = rev(c(der(sum(bL>0)),izq(sum(bL<0)))), limits=c(max[1],max[6]))
plot

png(paste(getwd(),"plots",paste(h,"confrontation_tree.png",sep = "_"),sep = "/"), 
    width = 3500, height = 3500, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot); dev.off()


######




