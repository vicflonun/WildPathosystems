#Run diversity analysis and plots for figure 1 and Supplementary Figure 2

#PACKAGES 
library(vegan) 
library(ggplot2)
library(Rmisc)
library(FSA)

# Set PARAMETERS
amplicon=c("16S","ITS2")[2]
method=c("clust","denoise")[1]
data=c("zymo")[1]  

#Set DIRECTORIES
setwd(paste("~/Documents/grass.microbiome.2",amplicon,"result_final",sep = "/"))
result.dir=paste("~/Documents/grass.microbiome.2","Results_final",amplicon,sep = "/") ; dir.create(result.dir)
plot.dir=paste(result.dir,"plot_paper_mock",sep = "/");dir.create(plot.dir)

#Set GRAPHIC parameters 
colores=c("grey50","#29854E","#B8C466","#A02355","#C55DA1","#234EA0","#008CCF")
colores.t=c("#54C073","#B80C93","#EF4035")
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

# Load METADATA
met=read.delim(file=paste("~/Documents/grass.microbiome.2/metadata.table",data,"all","txt",sep = "."),sep = ",")
rownames(met)=met$sample.id
met[met$time=="t0","treatment"]="Mock"
met$treatment=factor(met$treatment,levels = c("Mock","Zt469","Zt549","Zpa796","Zpa21","IPO323","IPO323koAVRStb6"))
met$host.id=factor(met$host.id,levels = c("Ac","Hm","OB","CS","CSkoStb6"))
met$time=factor(met$time, levels = c("t0","t2","t3","t4","t7","t8","t11"))

# Load DATA
if (method == "clust" | method == "denoise"){
  if (amplicon=="16S"){
    #16S
    otu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(result.dir,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
    core_t=c(50)
  } else if (amplicon=="ITS2"){
    #ITS2
    otu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(result.dir,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
    core_t=c(25)
    
    #Remove the few Zymoseptoria reads as they might represent contamination in Mock samples
    z=tax[tax$genus %in% "g:Zymoseptoria", "otu.id"]
    otu=otu[!(rownames(otu) %in% z),]
    tax=tax[!(tax$genus %in% "g:Zymoseptoria"),]
    
  }
}

# Data FILTERING 
met.2=met[met$host.id %in% c("Ac","Hm") & met$treatment %in% c("Mock"),] #filter samples to analyze
otu.2=otu[,colnames(otu) %in% met.2$sample.id]
met.2=met.2[rownames(met.2) %in% colnames(otu.2),]

otu.2=otu.2[,colSums(otu.2)>100] #remove low count samples
otu.2=otu.2[rowSums(otu.2)>0,] #remove zero count OTU
met.2=met.2[colnames(otu.2),] 
tax=tax[rownames(otu.2),]
  
#Normalize by SUBSAMPLING
set.seed(23171341)
print(min(colSums(otu.2))) #min read count
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) 
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#Calculate alpha-diversity
ric=colSums(otu.s>0) #richness
shan=diversity(t(otu.s),index = "shannon") #Shannon

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan, total=colSums(otu.2))
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c("plate.number","plate.id","well","sample.id","host","sample.name", "host.id","treatment","time","individual")])

#Plot alpha diversity
a=ggplot(alpha, aes(x=host.id, y=richness, fill=host.id))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=host.id))+
  scale_fill_manual(values=c(colores.h))+ylab("Observed richness")+
  tema+ggtitle(amplicon)
a
shapiro.test(alpha$richness)
wilcox.test(alpha[alpha$host.id=="Ac","richness"], alpha[alpha$host.id=="Hm","richness"])

png(paste(plot.dir,paste("p1",method,amplicon,"richness.host.id.png",sep = "."),sep = "/"), 
    width = 1000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

a=ggplot(alpha, aes(x=host.id, y=shannon, fill=host.id))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=host.id))+
  scale_fill_manual(values=c(colores.h))+ylab("Shannon index")+
  tema+ggtitle(amplicon)+scale_y_continuous(limits = c(1,5))
a
shapiro.test(alpha$shannon)
wilcox.test(alpha[alpha$host.id=="Ac","shannon"], alpha[alpha$host.id=="Hm","shannon"])

png(paste(plot.dir,paste("p1",method,amplicon,"shannon.host.id.png",sep = "."),sep = "/"), 
    width = 800, height = 600, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

a=ggplot(alpha, aes(x=host.id, y=richness, fill=host.id))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=host.id))+
  scale_fill_manual(values=c(colores.h))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
  tema+ggtitle(amplicon)+ylab("Observed richness")+
  facet_wrap(~time,scales="free_x",nrow = 1)
a


png(paste(plot.dir,paste("p1",method,amplicon,"richness.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


a=ggplot(alpha, aes(x=host.id, y=shannon, fill=host.id))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=host.id))+
  scale_fill_manual(values=c(colores.h))+#geom_hline(yintercept = 100, linetype=2,color="grey50")+
  tema+ggtitle(amplicon)+ylab("Observed shannon")+
  facet_wrap(~time,scales="free_x",nrow = 1)
a

png(paste(plot.dir,paste("p1",method,amplicon,"shannon.all.png",sep = "."),sep = "/"), 
    width = 2000, height = 800, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

#Analysis by time
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

#Analysis by host
  for (j in unique(alpha[alpha$host.id==i, "time"])){
      sub=alpha[alpha$time==j,]
      test=wilcox.test(shannon~host.id, data = sub)
      print(paste(j))
      print(test$p.value)
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

a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=treatment,fill=host.id))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=3, alpha=0.85, color="black")+tema+
  scale_fill_manual(values=c(colores.h))+ylab("NMDS2")+xlab("NMDS1")+
  scale_color_manual(values=c(colores.h))+
  scale_shape_manual(values=c(21:25,21))+
  ggtitle(amplicon, "NMDS-Bray")

a

png(paste(plot.dir,paste(method,amplicon,"wild..CSS.NMDS.treatment.png",sep = "."),sep = "/"), 
    width = 900, height = 700, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()


a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=treatment,fill=time))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=3, alpha=0.85, color="black")+tema+
  scale_fill_manual(values=hcl.colors(5,palette = "Hawaii",rev = T)[1:5])+
  ylab("NMDS2")+xlab("NMDS1")+
  #scale_color_manual(values=hcl.colors(5,palette = "Hawaii",rev = T)[1:5])+
  scale_shape_manual(values=c(21:25,21))+
  ggtitle(amplicon, "NMDS-Bray")

a

png(paste(plot.dir,paste(method,amplicon,"wild..CSS.NMDS.time.png",sep = "."),sep = "/"), 
    width = 1100, height = 900, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

#Permanova
scaling=vegdist(t(mr), method = "bray", binary = F) #calculate distance
per=adonis2(scaling~host.id*time, data = met.2, permutations = 1000)
per
write.table(per, paste(plot.dir,paste(amplicon,method,"permanova.CSS.host.treatment.time.txt",sep="."),sep="/"),sep = "\t")


##VEEN diagrams
#Packages
library(phyloseq)
library(microbiome)

dir.1=paste(plot.dir,"core",sep = "/") #name for the main folder
dir.create(dir.1)

#Define levels of the variables
time=c("t0","t2","t4","t7","t11")
host=c("Ac","Hm")
treat=c("Mock")

#Define prevalence threshold (%)
core_t=c(20,50)[1]

#Create phyloseq object and set levels
g=phyloseq(otu_table(as.matrix(otu.s), taxa_are_rows = T), sample_data(met.2), tax_table(as.matrix(tax)))
sample_data(g)$treatment = factor(sample_data(g)$treatment, levels = treat)
sample_data(g)$time = factor(sample_data(g)$time,levels = time)
sample_data(g)$host.id = factor(sample_data(g)$host.id,levels = host)
ranks=c("domain","phylum","class","order","family","genus")

#Ac core
g1=subset_samples(g, host.id %in% "Ac")
g2 <- prune_taxa(taxa_sums(g1) > 0, g1)
g3 <- microbiome::transform(g2, "compositional")
core.taxa.standard <- core_members(g3, detection = 0.001, prevalence = core_t/100)
ac_core=core.taxa.standard

counts=otu_table(g3)
summary(colSums(counts[ac_core,])*100)
sort(rowSums(counts[ac_core,])*100)


#hm core
g1=subset_samples(g, host.id %in% "Hm")
g2 <- prune_taxa(taxa_sums(g1) > 0, g1)
g3 <- microbiome::transform(g2, "compositional")

core.taxa.standard <- core_members(g3, detection = 0.001, prevalence = core_t/100)
hm_core=core.taxa.standard

counts=otu_table(g3)
summary(colSums(counts[hm_core,])*100)
sort(rowSums(counts[hm_core,])*100)


#Venn diagram
#library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)

s=list(ac=ac_core,hm=hm_core)
venn=Venn(s)
plotdata=process_data(venn)

a=ggvenn(s, fill_color = colores.h,
  #stroke_color = "black"
  stroke_size = 0.0, set_name_size = 4, text_size = 7)
a
png(paste(dir.1,paste(amplicon,core_t,"core.host.png",sep = "."),sep = "/"), 
    width = 2000, height = 1500, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

tax[overlap(venn, 1:2),] # core all
write.table(tax[overlap(venn, 1:2),], paste(dir.1,paste(amplicon,core_t,"core_all",sep="."),sep="/"),sep = "\t")

tax[discern(venn, 1) ,] #unique ac
write.table(tax[discern(venn, 1) ,], paste(dir.1,paste(amplicon,core_t,"unique_AC_Mock",sep="."),sep="/"),sep = "\t")

tax[discern(venn, 2) ,] #unique hm
write.table(tax[discern(venn, 2) ,], paste(dir.1,paste(amplicon,core_t,"unique_HM_mock",sep="."),sep="/"),sep = "\t")


library(UpSetR)
library(colorspace)

input=c(
  Ac=dim(tax[discern(venn, 1) ,])[1],
  Hm=dim(tax[discern(venn, 2) ,])[1],
  "Ac&Hm"=dim(tax[overlap(venn, 1:2),])[1]
)

a=upset(fromExpression(input), 
      nintersects = NA, 
      nsets = 3, 
      #keep.order = T,
      #show.numbers = T,
      order.by = "degree", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 0.9, 
      point.size = 2, 
      line.size = 1,
      main.bar.color = c("black"),
      sets.bar.color =colores.h
)
a

png(paste(dir.1,paste(amplicon,core_t,"core.host_upset.png",sep = "."),sep = "/"), 
    width = 700, height = 500, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

#Barplot
detach("package:Rmisc", unload = TRUE)
detach("package:plyr", unload = TRUE)
#barplot
colores2=c(#"#A74476","#F56AB0",
           #"#804795","#DF87FF",
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
t="phylum"
otu.s2=as.data.frame(cbind(otu.s,tax))
otu.t=gather(data=otu.s2, key = "sample", value = "abs", colnames(otu.2))
otu.t[is.na(otu.t[,t]),t]="Unclassified"

for (i in met.2$sample.id){
  otu.t[otu.t$sample == i, c( "host.id","time","treatment","plate.id")] = met.2[i,c("host.id","time","treatment","plate.id" )]
}

otu.t0 = otu.t %>% group_by(host.id,treatment) %>% summarise(absab = sum(abs))
otu.t1 = otu.t %>% group_by(otu.t[,t],host.id,treatment) %>% summarise(absab = sum(abs)) 

colnames(otu.t1)[1]=c("rank")
sum(otu.t1$host.id==otu.t0$host.id, na.rm = T)==dim(otu.t1)[1]
sum(paste(otu.t1$host.id,otu.t1$rank)==paste(otu.t0$host.id,otu.t1$rank))==dim(otu.t1)[1]
otu.t1$relab=otu.t1$absab/otu.t0$absab*100

#Reorder taxa if needed
orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
low = orden[orden$relab < 0.1, "rank"]
otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "Low abundant"
n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 
#orden=orden[order(otu.t$relab,decreasing = T),]

tema$legend.position="right"
otu.t1$rank=gsub("p:","",otu.t1$rank)
a=ggplot(data=otu.t1, aes(y=relab, x=host.id, fill=rank))+
  geom_bar(stat="identity",linewidth = 0.95,size=0.5)+tema+
  scale_fill_manual(values = c(colores2,"grey20","grey80"))+
  ylab("Relative abundance (%)")#+
  #facet_wrap(~host.id*treatment, scales="free", ncol = 3)
a

png(paste(plot.dir,paste(method,amplicon,t,"barplot.onlyhost_legend.png",sep = "."),sep = "/"), 
    width = 1800, height = 1000,units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

tema$legend.position="none"
a=ggplot(data=otu.t1, aes(y=relab, x=host.id, fill=rank))+
  geom_bar(stat="identity",linewidth = 0.95,size=0.5)+tema+
  scale_fill_manual(values = c(colores2,"grey20","grey80"))+
  ylab("Relative abundance (%)")#+
#facet_wrap(~host.id*treatment, scales="free", ncol = 3)
a

png(paste(plot.dir,paste(method,amplicon,t,"barplot.onlyhost_nolegend.png",sep = "."),sep = "/"), 
    width = 400, height = 500,units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
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
host=c("Ac","Hm")
treat=c("Mock","Compa","Incom")


#Create phyloseq object and set levels
g=phyloseq(otu_table(as.matrix(otu.2), taxa_are_rows = T), sample_data(met.2), tax_table(as.matrix(tax)))
sample_data(g)$treatment = factor(sample_data(g)$treatment, levels = treat)
sample_data(g)$time = factor(sample_data(g)$time,levels = time)
sample_data(g)$host.id = factor(sample_data(g)$host.id,levels = host)
ranks=c("domain","phylum","class","order","family","genus")

#Set results directory
dir.1=paste(plot.dir,"enrichments",sep = "/") #name for the main folder
dir.create(dir.1)

#Select taxonomy to test
t=NULL ; tn="otu.id"
t=ranks[6] ; tn=t

#Frequency cut off for ANCOMBC2
cut=core_t/100

#Select the hosts for comparison and the treatments to analyze
#for more than 3 levels some adjustments need to be done
hosts=host
treats=treat[1]

#### Comparison between cultivars despite time ####
#Create directory for this comparison 
dir.2=paste(dir.1,"between_hosts_all_times",sep = "/") #Dont change this name
dir.create(dir.2)

for (i in 1:length(treats)){
  
  #Create directory for this comparison 
  dir.3=paste(dir.2,paste(paste(hosts,collapse = "_"),treats[i],sep = "_in_"),sep = "/") #Dont change this name
  dir.create(dir.3, showWarnings = F)
  
  #Subset phyloseq objects
  g1=subset_samples(g, host.id %in% hosts & treatment %in% treats[i])
  
    tg=ancombc2(data=g1,
                assay_name = "wild_grass",
                p_adj_method = "BH",
                prv_cut = cut,
                lib_cut = 0, 
                group = "host.id",
                struc_zero = F,
                fix_formula = "host.id",
                alpha = 0.05,
                verbose = T,
                s0_perc=0.05,
                tax_level = t
                #max_iter = 100,conserve =T ,global = T,pairwise = T,dunnet = T, formula = form,
    )
    
    final=data.frame(tg$res[,c(1,3,9,11,13)], enriched="no_diff")
    if (tn=="otu.id"){final=cbind(final,tax[tg$res[,1],])}
    
    colnames(final)[1:5]=c("taxon","lfc","p_una","p_adj","diff")
    final[final$diff==TRUE & final$lfc < 0, "enriched"]=hosts[1]
    final[final$diff==TRUE & final$lfc > 0, "enriched"]=hosts[2]
    
    write.table(final, paste(dir.3,paste(tn,"DA_ANCOMBC2.txt",sep = "_"),sep = "/"),
                quote = F, sep = "\t", row.names = F, col.names = T)
    

    a=ggplot(data=final, aes(x=lfc, y=-log10(p_adj), color=enriched))+
      geom_point(size=3)+scale_color_manual(values=c(colores.h,"grey50"))+
      geom_hline(yintercept = -log10(0.05), linetype=2)+
      geom_vline(xintercept = c(0.5,-0.05), linetype=2)+
      tema+
      scale_x_continuous(limits = c(-3.5,3.5))
    
    a
    png(paste(dir.2,paste(core_t,tn,"DA_ANCOMBC2_volcano.png",sep = "_"),sep = "/"), 
        width = 1100, height = 700, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
    print(a) ; dev.off()
    
    
    
}




