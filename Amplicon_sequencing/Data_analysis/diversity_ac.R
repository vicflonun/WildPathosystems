#Run diversity analysis for Aegilops cylindrica infected plants. 
#Create panels for Figure 2 Supplementary Figure 4 a 7. 

#Packages
library(vegan) 
library(ggplot2)
library(Rmisc)
library(FSA)

# Set PARAMETERS
amplicon=c("16S","ITS2")[1]
method=c("clust","denoise")[1]
data=c("zymo")[1]  

#Directories
#Set result DIRECTORIES
setwd(paste("~/Documents/grass.microbiome.2",amplicon,"result_final",sep = "/"))
result.dir=paste("~/Documents/grass.microbiome.2","Results_final",amplicon,sep = "/") ; dir.create(result.dir)
plot.dir=paste(result.dir,"plot_paper_treatment_ac",sep = "/");dir.create(plot.dir)

#Graphics
colores=c("grey50","#29854E","#B8C466","#A02355","#C55DA1","#234EA0","#008CCF")
colores.t=c("#54C073","#B80C93","#EF4035")
colores.h2=c("#29854E","#74A557","#B8C466")
colores.h=c("#FFC850","#116875")
colores.i=c("#76e5a5","#41c6a8","#18a6a2","#2a6678","#2f4858")
scaleFUN <- function(x) sprintf("%.2f", x) #Put 2 extra decimal zeros

tema=theme(axis.text.x = element_text(color="black",size=8, angle=0,hjust=0.5,vjust=0.5,family = "Arial"),
           axis.text.y = element_text(color="black",size=8,family = "Arial"),
           axis.title = element_text(color="black",size=8,family = "Arial"),
           legend.text = element_text(color = "black",size=8,family = "Arial"),
           legend.key.size = unit(0.3,"cm"),
           legend.title = element_text(color="black",size=8,family = "Arial"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=6, color="black",family = "Arial", margin = margin(0,0,0,0,"cm")),
           strip.text.y = element_text(size=6, color="black",family = "Arial"),
           panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")

# Load metadata
met=read.delim(file=paste("~/Documents/grass.microbiome.2/metadata.table",data,"all","txt",sep = "."),sep = ",")
rownames(met)=met$sample.id
met[met$time=="t0","treatment"]="Mock"
met$treatment=factor(met$treatment,levels = c("Mock","Zt469","Zt549","Zpa796","Zpa21","IPO323","IPO323koAVRStb6"))
met$host.id=factor(met$host.id,levels = c("Ac","Hm","OB","CS","CSkoStb6"))
met$time=factor(met$time, levels = c("t0","t2","t3","t4","t7","t8","t11"))

# Load data
if (method == "clust" | method == "denoise"){
  if (amplicon=="16S"){
    #16s
    otu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(result.dir,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
    core_t=c(50)
  } else if (amplicon=="ITS2"){
    #ITS2
    otu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(result.dir,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
    core_t=c(20)
  }
}

# Data FILTERING 
met.2=met[met$host.id %in% c("Ac"),]
otu.2=otu[,colnames(otu) %in% met.2$sample.id]
met.2=met.2[rownames(met.2) %in% colnames(otu.2),]

otu.2=otu.2[,colSums(otu.2)>100]
otu.2=otu.2[rowSums(otu.2)>0,]
met.2=met.2[colnames(otu.2),]
tax=tax[rownames(otu.2),]
  
#Normalize by SUBSAMPLING
set.seed(23171341)
print(min(colSums(otu.2)))
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) 
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#Calculate alpha-diversity
ric=colSums(otu.s>0)
shan=diversity(t(otu.s),index = "shannon")

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan, total=colSums(otu.2))
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c("plate.number","plate.id","well","sample.id","host",
                            "sample.name", "host.id","treatment","time","individual")])


a=ggplot(alpha, aes(x=treatment, y=richness, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=c(colores.t))+ylab("Observed richness")+
  tema+ggtitle(amplicon)
a
shapiro.test(alpha$richness)
kruskal.test(richness~treatment, data = alpha)
dunnTest(richness~treatment, data = alpha)


png(paste(plot.dir,paste("p1",method,amplicon,"richness.treatment.png",sep = "."),sep = "/"), 
    width = 1000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

a=ggplot(alpha, aes(x=treatment, y=shannon, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=c(colores.t))+ylab("Shannon index")+
  tema+ggtitle(amplicon)+scale_y_continuous(limits = c(1,6))
a
shapiro.test(alpha$shannon)
kruskal.test(shannon~treatment, data = alpha)
dunnTest(shannon~treatment, data = alpha)


png(paste(plot.dir,paste("p1",method,amplicon,"shannon.treatment.png",sep = "."),sep = "/"), 
    width = 800, height = 600, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

a=ggplot(alpha, aes(x=time, y=richness, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=c(colores.t[c(1:3,2:3)]))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
  tema+ggtitle(amplicon)+ylab("Observed richness")+
  facet_wrap(~host.id*treatment,scales="free_x",nrow = 1)
a

png(paste(plot.dir,paste("p1",method,amplicon,"richness.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

a=ggplot(alpha, aes(x=time, y=shannon, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=c(colores.t[c(1:3,2:3)]))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
  tema+ggtitle(amplicon)+ylab("Shannon")+
  facet_wrap(~host.id*treatment,scales="free_x",nrow = 1)
a

png(paste(plot.dir,paste("p1",method,amplicon,"shannon.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


a=ggplot(alpha, aes(x=treatment, y=richness, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=c(colores.t[c(1:3,2:3)]))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
  tema+ggtitle(amplicon)+ylab("Shannon")+
  facet_wrap(~host.id*time,scales="free_x",nrow = 1)
a

png(paste(plot.dir,paste("p1",method,amplicon,"richness2.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


a=ggplot(alpha, aes(x=treatment, y=shannon, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=c(colores.t[c(1:3,2:3)]))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
  tema+ggtitle(amplicon)+ylab("Shannon")+
  facet_wrap(~host.id*time,scales="free_x",nrow = 1)
a

png(paste(plot.dir,paste("p1",method,amplicon,"shannon2.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

#Analysis between times by treatment
for ( i in unique(alpha$host.id)){
  for (j in unique(alpha[alpha$host.id==i, "treatment"])){
    sub=alpha[alpha$host.id==i & alpha$treatment==j,]
    if (j!="Mock"){
      sub=rbind(sub,alpha[alpha$host.id==i & alpha$treatment=="Mock" & alpha$time=="t0",])
    }
    test=kruskal.test(shannon~time, data = sub)
    print(paste(i,j))
    print(test$p.value)
    test=dunnTest(shannon~time, data = sub, method = "bh")$res
    print(test[test[,4]<=0.05,])
  }
}

#Analysis between treatments by time
for ( i in unique(alpha$host.id)){
  for (j in unique(alpha[alpha$host.id==i, "time"])){
    if(j != "t0"){
      sub=alpha[alpha$host.id==i & alpha$time==j,]
      test=kruskal.test(shannon~treatment, data = sub)
      print(paste(i,j))
      print(test$p.value)
      test=dunnTest(shannon~treatment, data = sub, method = "bh")$res
      print(test[test[,4]<=0.05,])
      
    } else if (j =="t0"){print("no treatments at t0")}
  }
  
}

#Analysis with no Zymoseptoria reads. 
if (amplicon=="ITS"){
  
  z=tax[tax$genus %in% "g:Zymoseptoria", "otu.id"]
  otu.z=otu.s[!(rownames(otu.s) %in% z),]
  tax.z=tax[!(tax$genus %in% "g:Zymoseptoria"),]
  
  #Calculate alpha-diversity
  ric=colSums(otu.z>0)
  shan=diversity(t(otu.z),index = "shannon")
  
  #Generate plot data
  alpha.z=data.frame(richness=ric, shannon=shan, total=colSums(otu.2))
  sum(rownames(alpha.z)==rownames(met.2))==dim(alpha.z)[1]
  alpha.z=cbind(alpha.z, met.2[,c("plate.number","plate.id","well","sample.id","host",
                              "sample.name", "host.id","treatment","time","individual")])
  
  
  a=ggplot(alpha.z, aes(x=treatment, y=richness, fill=treatment))+
    geom_boxplot(outlier.colour = NA)+
    geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
    scale_fill_manual(values=c(colores.t[c(1:3,2:3)]))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
    tema+ggtitle(amplicon)+ylab("richness")+
    facet_wrap(~host.id*time,scales="free_x",nrow = 1)
  a
  
  png(paste(plot.dir,paste("p1",method,amplicon,"richness3.all.png",sep = "."),sep = "/"), 
      width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
  print(a) ; dev.off()
  
  
  a=ggplot(alpha.z, aes(x=treatment, y=shannon, fill=treatment))+
    geom_boxplot(outlier.colour = NA)+
    geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
    scale_fill_manual(values=c(colores.t[c(1:3,2:3)]))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
    tema+ggtitle(amplicon)+ylab("shannon")+
    facet_wrap(~host.id*time,scales="free_x",nrow = 1)
  a
  
  
  png(paste(plot.dir,paste("p1",method,amplicon,"shannon3.all.png",sep = "."),sep = "/"), 
      width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
  print(a) ; dev.off()
  
  
}


#Beta diversity
library(metagenomeSeq)

#normalization of count reads using CSS 
mp=newMRexperiment(as(otu.2, "matrix"))
mr=MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),norm=TRUE, log=TRUE )

#calcultae distances
set.seed(23171341)
scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
scaling2=metaMDS(scaling) ; scaling2$stress
scaling3=data.frame(scaling2$points) #select coordinates
scaling3=cbind(scaling3,met.2[rownames(scaling3),])

a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=treatment,fill=treatment))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=3, alpha=0.85, color="black")+tema+
  scale_fill_manual(values=c(colores.t))+ylab("NMDS2")+xlab("NMDS1")+
  scale_color_manual(values=c(colores.t))+
  scale_shape_manual(values=c(21:25,21))+
  ggtitle(amplicon, "NMDS-Bray")

a

png(paste(plot.dir,paste(method,amplicon,"wild..CSS.NMDS.treatment.png",sep = "."),sep = "/"), 
    width = 900, height = 700, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=treatment,fill=time))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=3, alpha=1, color="black")+tema+
  scale_fill_manual(values=c(colores.i))+ylab("NMDS2")+xlab("NMDS1")+
  scale_color_manual(values=c(colores.i))+
  scale_shape_manual(values=c(21:25,21))+
  ggtitle(amplicon, "NMDS-Bray")

a

png(paste(plot.dir,paste(method,amplicon,"wild..CSS.NMDS.time.png",sep = "."),sep = "/"), 
    width = 1100, height = 900, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


#Permanova
scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
per=adonis2(scaling~time*treatment, data = met.2, permutations = 1000)
per
write.table(per, paste(plot.dir,paste(amplicon,method,"permanova.CSS.host.treatment.time.txt",sep="."),sep="/"),sep = "\t")


  for (j in levels(met.2$time)[c(2,4,5,7)]){
    print(j)
    met.s=met.2[met.2$time==j,]
    mp.s=newMRexperiment(as(otu.2[,rownames(met.s)], "matrix"))
    mr=MRcounts( cumNorm( mp.s, p=cumNormStat( mp.s ) ),norm=TRUE, log=TRUE )
    set.seed(23171341)
    scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
    per=adonis2(scaling~treatment, data = met.s, permutations = 1000)
    print(per)
    write.table(per, paste(plot.dir,paste(amplicon,method,j,"wild.permanova.treatment.txt",sep="."),sep="/"),sep = "\t")
  }


#t7

#normalization of count reads using CSS 
met.s=met.2[met.2$time=="t7",]
mp.s=newMRexperiment(as(otu.2[,rownames(met.s)], "matrix"))
mr=MRcounts( cumNorm( mp.s, p=cumNormStat( mp.s ) ),norm=TRUE, log=TRUE )

#calcultae distances
set.seed(23171341)
scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
scaling2=metaMDS(scaling) ; scaling2$stress
scaling3=data.frame(scaling2$points) #select coordinates
scaling3=cbind(scaling3,met.2[rownames(scaling3),])

a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=treatment,fill=treatment))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=3, alpha=0.85, color="black")+tema+
  scale_fill_manual(values=c(colores.t))+ylab("NMDS2")+xlab("NMDS1")+
  scale_color_manual(values=c(colores.t))+
  scale_shape_manual(values=c(21:25,21))+
  ggtitle(amplicon, "NMDS-Bray")

a

png(paste(plot.dir,paste(method,amplicon,"wild.t7.CSS.NMDS.treatment.png",sep = "."),sep = "/"), 
    width = 1100, height = 900, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

if (amplicon == "ITS"){
  
  z=tax[tax$genus %in% "g:Zymoseptoria", "otu.id"]
  otu.z=otu.2[!(rownames(otu.2) %in% z),]
  tax.z=tax[!(tax$genus %in% "g:Zymoseptoria"),]
  
  #normalization of count reads using CSS 
  mp=newMRexperiment(as(otu.z, "matrix"))
  mr=MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),norm=TRUE, log=TRUE )
  
  #calcultae distances
  set.seed(23171341)
  scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
  scaling2=metaMDS(scaling) ; scaling2$stress
  scaling3=data.frame(scaling2$points) #select coordinates
  scaling3=cbind(scaling3,met.2[rownames(scaling3),])
  
  a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=treatment,fill=treatment))+
    geom_hline(yintercept = 0, linetype=2)+
    geom_vline(xintercept = 0,linetype=2)+
    geom_point(size=3, alpha=0.85, color="black")+tema+
    scale_fill_manual(values=c(colores.t))+ylab("NMDS2")+xlab("NMDS1")+
    scale_color_manual(values=c(colores.t))+
    scale_shape_manual(values=c(21:25,21))+
    ggtitle(amplicon, "NMDS-Bray")
  
  a
  
  png(paste(plot.dir,paste(method,amplicon,"wild..CSS.NMDS.treatment_no zymo.png",sep = "."),sep = "/"), 
      width = 1100, height = 900, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
  print(a) ; dev.off()
  
  
  
  #Permanova
  scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
  per=adonis2(scaling~time*treatment, data = met.2, permutations = 1000)
  per
  write.table(per, paste(plot.dir,paste(amplicon,method,"permanova.CSS.host.treatment.time_nozymo.txt",sep="."),sep="/"),sep = "\t")
  
  
  for (j in levels(met.2$time)[c(2,4,5,7)]){
    print(j)
    met.s=met.2[met.2$time==j,]
    mp.s=newMRexperiment(as(otu.z[,rownames(met.s)], "matrix"))
    mr=MRcounts( cumNorm( mp.s, p=cumNormStat( mp.s ) ),norm=TRUE, log=TRUE )
    set.seed(23171341)
    scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
    per=adonis2(scaling~treatment, data = met.s, permutations = 1000)
    print(per)
    write.table(per, paste(plot.dir,paste(amplicon,method,j,"wild.permanova.treatment_nozymo.txt",sep="."),sep="/"),sep = "\t")
  }
  
}


#Barplot
detach("package:Rmisc", unload = TRUE)
detach("package:plyr", unload = TRUE)
#barplot
colores2=c("#A74476","#F56AB0",
           "#804795","#DF87FF",
           "#524CA1","#9A94EF",
           "#2C79AE","#6EC4FF",
           "#369DA4","#6EF6FF",
           "#29A876","#43FFB5",
           "#55A054","#A5F77A",
           "#A5A53F","#F1F18A",
           "#AA723D","#FFBB7A",
           "#A12F2F","#FF4848")
library(dplyr)
library(tidyr)
t="class"
otu.s2=as.data.frame(cbind(otu.s,tax))
otu.t=gather(data=otu.s2, key = "sample", value = "abs", colnames(otu.2))
otu.t[is.na(otu.t[,t]),t]="Unclassified"
for (i in met.2$sample.id){
  otu.t[otu.t$sample == i, c( "host.id","time","treatment","plate.id")] = met.2[i,c("host.id","time","treatment","plate.id" )]
}

otu.t0 = otu.t %>% group_by(host.id,treatment,time) %>% summarise(absab = sum(abs))  #Sum all the reads in each combination of conditions
otu.t1 = otu.t %>% group_by(otu.t[,t],host.id,treatment,time) %>% summarise(absab = sum(abs)) #Summ all the reads in each combination of conditions and taxa
colnames(otu.t1)[1]=c("rank")
sum(otu.t1$host.id==otu.t0$host.id, na.rm = T)==dim(otu.t1)[1] #All combinations of conditions should match (TRUE)
sum(paste(otu.t1$host.id,otu.t1$rank)==paste(otu.t0$host.id,otu.t1$rank))==dim(otu.t1)[1] #All combinations of conditions should match (TRUE)
otu.t1$relab=otu.t1$absab/otu.t0$absab*100 #Generate relative abundance data

#Reorder taxa if needed
orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
low = orden[orden$relab < 0.25, "rank"]
otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "Low abundant"
n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 
#orden=orden[order(otu.t$relab,decreasing = T),]

tema$legend.position="none"
a=ggplot(data=otu.t1, aes(y=relab, x=time, fill=rank))+
  geom_bar(stat="identity",linewidth = 0.95,size=0.5)+tema+
  scale_fill_manual(values = c(colores2,"grey80","grey20","black"))+
  ylab("Relative abundance (%)")+
  facet_wrap(~host.id*treatment, scales="free", ncol = 3)
a

png(paste(plot.dir,paste(method,amplicon,t,"barplot.all_no_legend.png",sep = "."),sep = "/"), 
    width = 2000, height = 1500,units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


tema$legend.position="none"
a=ggplot(data=otu.t1, aes(y=relab, x=treatment, fill=rank))+
  geom_bar(stat="identity",linewidth = 0.95,size=0.5)+tema+
  scale_fill_manual(values = c(colores2,"grey80","grey20","black"))+
  ylab("Relative abundance (%)")+
  facet_wrap(~host.id*time, scales="free", ncol = 5)
a

png(paste(plot.dir,paste(method,amplicon,t,"barplot.all_no_legend2.png",sep = "."),sep = "/"), 
    width = 2500, height = 1500,units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


tema$legend.position="right"
a=ggplot(data=otu.t1, aes(y=relab, x=time, fill=rank))+
  geom_bar(stat="identity",linewidth = 0.95,size=0.5)+tema+
  scale_fill_manual(values = c(colores2,"grey80","grey20","black"))+
  ylab("Relative abundance (%)")+
  facet_wrap(~host.id*treatment, scales="free", ncol = 3)
a

png(paste(plot.dir,paste(method,amplicon,t,"barplot.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 2000,units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()



#Differentially abundant OTUs

#### Enrichments
#Packages
library(ANCOMBC)
library(phyloseq)
library(vegan)
library(ggplot2)


#Define levels of the variables
time=c("t0","t2","t4","t7","t11")
host=c("Ac")
treat=c("Mock","Compa","Incom")


met.3=met.2
met.3$treatment=as.character(met.3$treatment)
met.3[met.3$treatment=="Zt469","treatment"]="Compa"
met.3[met.3$treatment=="Zt549","treatment"]="Incom"
met.3[met.3$treatment=="Zpa21","treatment"]="Compa"
met.3[met.3$treatment=="Zpa796","treatment"]="Incom"
met.3$treatment=factor(met.3$treatment,levels = c("Mock","Compa","Incom"))


#Create phyloseq object and set levels
g=phyloseq(otu_table(as.matrix(otu.2), taxa_are_rows = T), sample_data(met.3), tax_table(as.matrix(tax)))
sample_data(g)$treatment = factor(sample_data(g)$treatment, levels = treat)
sample_data(g)$time = factor(sample_data(g)$time,levels = time)
sample_data(g)$host.id = factor(sample_data(g)$host.id,levels = host)
ranks=c("domain","phylum","class","order","family","genus")

g=subset_samples(g, time %in% c("t4","t7"))

#Set results directory
dir.1=paste(plot.dir,"enrichments",sep = "/") #name for the main folder
dir.create(dir.1)

#Select taxonomy to test
t=NULL ; tn="otu.id"
#t=ranks[6] ; tn=t

#Frequency cut off for ANCOMBC2
cut=core_t/100

#Select the hosts for comparison and the treatments to analyze
#for more than 3 levels some adjustments need to be done
hosts=host
treats=treat
times=time

#### Comparison between treatment despite time ####

dir.2=paste(dir.1,"between_treatments_all_times",sep = "/") #Dont change this name
dir.create(dir.2)
  
  #Subset phyloseq objects
  g1=subset_samples(g, host.id %in% hosts)
  
    tg=ancombc2(data=g1,
                assay_name = "wild_grass",
                p_adj_method = "BH",
                prv_cut = core_t/100,
                lib_cut = 0, 
                group = "treatment",
                struc_zero = F,
                fix_formula = "treatment",
                alpha = 0.05,
                verbose = T,
                s0_perc=0.05,
                tax_level = t
                #max_iter = 100,conserve =T ,global = T,pairwise = T,dunnet = T, formula = form,
    )
    
    final=data.frame(tg$res, enriched1="no_diff",enriched2="no_diff")
    if (tn=="otu.id"){final=cbind(final,tax[tg$res[,1],])}
    
   
    final[final$diff_treatmentCompa==TRUE & final$lfc_treatmentCompa < 0, "enriched1"]=treat[1]
    final[final$diff_treatmentCompa==TRUE & final$lfc_treatmentCompa > 0, "enriched1"]=treat[2]
    
    final[final$diff_treatmentIncom==TRUE & final$lfc_treatmentIncom < 0, "enriched2"]=treat[1]
    final[final$diff_treatmentIncom==TRUE & final$lfc_treatmentIncom > 0, "enriched2"]=treat[3]
    
    write.table(final, paste(dir.2,paste("4-7",tn,"DA_ANCOMBC2.txt",sep = "_"),sep = "/"),
                quote = F, sep = "\t", row.names = F, col.names = T)

    
    







