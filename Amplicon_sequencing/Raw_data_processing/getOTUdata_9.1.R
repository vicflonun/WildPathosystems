##### Generate files for data analysis #####

library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Biostrings)
library(splitstackshape)
library(Rmisc)
library(tidyr)
library(dplyr)
library(ggpubfigs)

# Set PARAMETERS
amplicon=c("16S","ITS2")[1]
method=c("clust","denoise")[1]
bootstrap=0.80 #not implemented yet
data=c("zymo")[1]  #what sub sample would you need
remove(c,e,u)

if ( amplicon == "16S"){ #minimum read count for a valid sample
  filter=350
} else if (amplicon == "ITS2"){
  filter=100 #minimum read count for a valid sample
}

#Set result DIRECTORIES
setwd(paste("~/Documents/grass.microbiome.2",amplicon,"result_final",sep = "/"))
result.dir=paste("~/Documents/grass.microbiome.2","Results_final",amplicon,sep = "/")
dir.create(result.dir)



# Select the FILES
if (method == "clust"){
  table=c("otutab.txt")
  fasta=c("otus.fasta")
  taxa=c(rdp="otus.sintax2")
} else if (method == "denoise"){
  table=c("zotutab.txt")
  fasta=c("zotus.fasta")
  taxa=c(rdp="zotus.sintax")
}

# Load OTU table
otu=read.delim(paste(getwd(),table,sep = "/"), header = T) #load
otu$X.OTU.ID=as.numeric(gsub("OTU","",otu$X.OTU.ID)) #change OTU names and sort
otu=otu[order(otu$X.OTU.ID),]
otu$X.OTU.ID=paste("OTU",otu$X.OTU.ID,sep="_")
rownames(otu)=otu$X.OTU.ID ; otu=otu[,-1]

# Load TAXA classification
ranks=c("otu.id","domain","phylum","class","order","family","genus")
sintax=read.delim(paste(getwd(),taxa,sep = "/"), header = F)
sintax=sintax[,-c(2,3,5)] #Format taxa
sintax=splitstackshape::cSplit(indt = sintax, splitCols = "V4", sep = ",", type.convert = F)
sintax=data.frame(sintax[,c(1:7)])
colnames(sintax)=ranks 

sintax$otu.id=as.numeric(gsub("OTU","",sintax$otu.id)) #change OTU names and sort
sintax=sintax[order(sintax$otu.id),]
sintax$otu.id=paste("OTU",sintax$otu.id,sep="_")
rownames(sintax)=sintax$otu.id ; sintax=sintax[,c(2:7,1)]
sintax=sintax[rownames(otu),]

# Load SEQUENCES
seq=readDNAStringSet(paste(getwd(),fasta,sep = "/"))
names(seq)=gsub("OTU","OTU_",names(seq)) #change names

#Load METADATA
meta=read.delim(file=paste("~/Documents/grass.microbiome.2/metadata.table",data,"all","txt",sep = "."),sep = ",")
rownames(meta)=meta$sample.id
meta=meta[colnames(otu),] #remove samples that did not got sequences
meta=meta[order(meta$sample.id),] #order samples
otu=otu[,order(colnames(otu))]
rownames(meta) == colnames(otu)

# FILTERING SAMPLES
total=colSums(otu) #calculate total read per sample 
total.t=data.frame(total,seq="umgc")
ggplot(data=total.t, aes(x=log10(total)))+
  geom_histogram()+ scale_x_continuous()+ylab("Number of samples")+xlab("log10(Total reads)")+
  #scale_x_continuous(limits = c(0,6), breaks = seq(0,6,0.5))+
  #scale_y_continuous(limits = c(0,90))+
  ggtitle(paste("Total samples ",dim(otu)[2]))

## Identify  non microbial sequences #adjustments needed if you are using another database
if ( amplicon == "16S"){
  e=unique(c(sintax[grep("Eukaryota",sintax$domain),"otu.id"],sintax[grep("Mitochondria",sintax$family), "otu.id"]))
  c=sintax[grep("Chloroplast",sintax$class),"otu.id"]
  print(paste("Chloroplast", length(c), "OTUs"))
  print(sum(rowSums(otu[c,])) / sum(total) *100)
  print(paste("Eukaryota", length(e), "OTUs"))
  print(sum(rowSums(otu[e,]),na.rm = T) / sum(total) * 100)
  
  u=sintax[is.na(sintax$domain),"otu.id"]
  print(paste("Unclassified", length(u), "OTUs"))
  print(sum(rowSums(otu[u,]),na.rm = T) / sum(total) * 100)
  
} else if (amplicon == "ITS2"){
  
  e=unique(sintax[grep("Viridiplantae",sintax$domain),"otu.id"])
  ec=unique(sintax[grep("Chlorophyta",sintax$phylum),"otu.id"])
  e=e[!e %in% ec]
  
  print(paste("Viridiplantae", length(e), "OTUs"))
  print(sum(rowSums(otu[e,])) / sum(total) * 100)
  c=c()
  
  u=sintax[is.na(sintax$domain),"otu.id"]
  print(paste("Unclassified", length(u), "OTUs"))
  print(sum(rowSums(otu[u,]),na.rm = T) / sum(total) * 100)
  
}


writeXStringSet(seq[u[1:100]], paste(result.dir,"unknown16s.fasta", format = "fasta"))


plant=t(t(otu[c(e,c),]) / (total) * 100)
print("Plant read distribution")
summary(colSums(plant))

meta.p=meta
plant=data.frame(total, host=meta.p$host)
summarySE(data=plant, "total", "host", na.rm = T)
write.table(summarySE(data=plant, "total", "host"),
            file = paste(result.dir,"plant.reads.txt",sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

# Remove not desired sequences and unclasified sequences
otu.a=otu[!rownames(otu) %in% unique(c(c,e,u)),] 
meta.a=meta[meta$sample.id %in% colnames(otu.a),]

## Recalculate total 
total.a=colSums(otu.a)

## Check read samples
read=data.frame(total=total.a,
                sample.name=meta.a$sample.name,
                host=meta.a$host,
                host.id=meta.a$host.id,
                plate=meta.a$plate.id,
                name=meta.a$treatment, 
                row.names =names(total.a))

plot=ggplot(data=read[,], aes(x=host, y=total, fill=host))+
  geom_boxplot()+
  scale_x_discrete(limits=c("Aegilops.cylindrica","Hordeum.murinum","Triticum.aestivum","Soil","Peat"),
                   labels=c("AC","HM","TA","soil","peat"))
  
plot

png(paste(result.dir,paste(method,amplicon,"reads.png",sep = "."),sep = "/"), 
    width = 1500, height = 750, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()

summarySE(data=read, "total", "host")
write.table(summarySE(data=read, "total", "host"),
            file = paste(result.dir,"noplant.reads.txt",sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

plot=ggplot(data=read, aes(x=log10(total), fill=host))+
  geom_histogram()+ scale_x_continuous()+ scale_fill_manual(values = friendly_pal("muted_nine"))+
  geom_vline(xintercept = log10(filter))
plot

#plot=ggplot(data=read, aes(x=log10(total), fill=host.id))+
 # geom_histogram()+ scale_x_continuous()+ scale_fill_manual(values = friendly_pal("muted_nine"))+
#  geom_vline(xintercept = log10(filter))
#plot

png(paste(result.dir,paste(method,amplicon,"reads2.png",sep = "."),sep = "/"), 
    width = 1000, height = 500, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()


plot=ggplot(data=read, aes(x=log10(total), fill=host.id))+
  geom_histogram()+ scale_x_continuous()+scale_fill_manual(limits=c("CA","CAKoStb6"), values = friendly_pal("muted_nine"))+
  geom_vline(xintercept = log10(filter))
plot



print(paste("Losing samples due to low read count", "less than", filter))
summary(as.factor(read[read$total<filter,"host"]))

## Remove low read samples
otus=otu.a[,colSums(otu.a)>=filter]
reads=read[colnames(otus),]
metas=meta[meta.a$sample.id %in% colnames(otus),]

##Removing low abundant a frequent OTUs
# Removing OTUs with cero reads in all samples, there might be some because we removed low read samples
otus=otus[rowSums(otus==0)<dim(otus)[2],]

#OTUs with just 1 read #singletons removed previously but after the filtering some migh be singletons now
otus=otus[rowSums(otus)>1,]
otus=otus[rowSums(otus>0)>1,]


## Optional: Compare the read number against the richness
rarecurve(t(otus), step = 2000, xlab = "Sample Size", ylab = "Species", label = F)
reads$richness=colSums(otus>0)

png(paste(result.dir,paste(method,amplicon,"rarecurve.png",sep = "."),sep = "/"), 
    width = 2000, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
rarecurve(t(otus), step = 2000, xlab = "Sample Size", ylab = "Species", label = F)
dev.off()


#reads$treatment=meta$treatment
plot=ggplot(data=reads,aes(x=total,y=richness))+
  geom_point(aes(color=host),size=2,alpha=0.5)+
  geom_smooth(method = "lm",color="black")+
  scale_x_continuous(breaks = seq(0,max(read$total),5000))
plot
png(paste(result.dir,paste(method,amplicon,"rerecurve2.png",sep = "."),sep = "/"), 
    width = 2000, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()


## Extract sequences
seqs=seq[rownames(otus)]

## Extract taxa
taxas=sintax[rownames(otus),]

## Save the data
write.table(otus,
            file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")

write.table(metas,
            file=paste(result.dir,paste("metadata.table",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")

write.table(taxas,
            file=paste(result.dir,paste("taxa.table",data,names(taxa),"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")

writeXStringSet(seqs, 
                file=paste(result.dir,paste("otus",data,"fasta",sep = "."),sep = "/"),
                format = "fasta")


############################




