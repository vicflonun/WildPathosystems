#Create randomized networks for the biotrophic infection of plants using data from 
#both virulent and avirulent isolates

#Packages
library(vegan) 
library(ggplot2)
library(Rmisc)
library(FSA)
library(cowplot)
library(SpiecEasi)
library(igraph)
library(Matrix)
library(GGally) 
library(ggrepel)


# Set PARAMETERS
amplicon=c("ITS2")[1]
method=c("clust","denoise","dada2")[1]
bootstrap=0.80 #not implemented yet
data=c("zymo")[1]  #what sub sample would you need

#Directories
#Set result DIRECTORIES
setwd(paste("~/Documents/grass.microbiome.2",amplicon,"result_final",sep = "/"))
result.dir=paste("~/Documents/grass.microbiome.2","Results_final",amplicon,sep = "/") ; dir.create(result.dir)
plot.dir=paste(result.dir,"plot_paper_networks_sum_final_test_new",sep = "/");dir.create(plot.dir)

#Graphics
colores=c("#A74476","#F56AB0",
          "#804795","#DF87FF",
          "#524CA1","#9A94EF",
          "#2C79AE","#6EC4FF",
          "#369DA4","#6EF6FF",
          "#29A876","#43FFB5",
          "#55A054","#A5F77A",
          "#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A",
          "#A12F2F","#FF4848")[c(2,6,10,14,18,4,8,12,16,20)]
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
met$host.id=factor(met$host.id,levels = c("Ac","Hm","OB","CS","CSkoStb6","CA","CAkoStb6"))
met$time=factor(met$time, levels = c("t0","t2","t3","t4","t7","t8","t11"))


#Network
# Set PARAMETERS
host=c("Aegilops.cylindrica","Hordeum.murinum")
ti=c("t0","t2","t4","t7","t11")[c(3:4)]


#data

#Directories
amplicon=c("16S","ITS2")[1]
result.dir=paste("~/Documents/grass.microbiome.2","Results_final",amplicon,sep = "/")
potu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
colnames(potu)=gsub("X","",colnames(potu))
rownames(potu)=paste("P",rownames(potu),sep = "")
ptax=read.delim(paste(result.dir,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
rownames(ptax)=paste("P",rownames(ptax),sep = "")
ptax[is.na(ptax$domain),"domain"]="d:unidentified_prokaryota"

amplicon=c("16S","ITS2")[2]
result.dir=paste("~/Documents/grass.microbiome.2","Results_final",amplicon,sep = "/") 
fotu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
colnames(fotu)=gsub("X","",colnames(fotu))
rownames(fotu)=paste("F",rownames(fotu),sep = "")
ftax=read.delim(paste(result.dir,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
rownames(ftax)=paste("F",rownames(ftax),sep = "")
ftax[is.na(ftax$domain),"domain"]="d:unidentified_eukaryota"
ftax[ftax$domain=="d:unidentified","domain"]="d:unidentified_eukaryota"

potu=potu[,colnames(potu) %in% colnames(fotu)]
fotu=fotu[,colnames(fotu) %in% colnames(potu)]

zall=paste("F",ftax[ftax$genus %in% "g:Zymoseptoria","otu.id"],sep = "")


for (ho in host){
  
  result.metrics=data.frame()
  result.zymo=data.frame()
  
  if (ho=="Aegilops.cylindrica"){
    treat=c("Mock","Zt469","Zt549")
    fotu=fotu[!rownames(fotu) %in% zall[-1],]
    ftax=ftax[!rownames(ftax) %in% zall[-1],]
    other=zall[-1] 
    zotu="FOTU_1"
    
  } else if (ho=="Hordeum.murinum") {
    treat=c("Mock","Zpa796","Zpa21") 
    
    FOTU_1=colSums(fotu[zall[c(2:5)],])
    fotu.2=fotu[!rownames(fotu) %in% zall,]
    fotu=rbind(FOTU_1=FOTU_1,fotu.2)
    ftax=ftax[!rownames(ftax) %in% zall[-1],]
    other=zall[-1] 
    zotu="FOTU_1"
    
  }
  #empty objects for metrics   
  dens=c();diam=c();modu=c();tran=c();mdeg=c();aspa=c();node=c();edge=c();
  bf=c();bfp=c();n.hub=c();bb=c();bbp=c();ff=c();ffp=c()
  
  #############for(tr in treat){
    
    
    samples=met[met$host %in% ho & met$time %in% ti,"sample.id"]
  
    print(samples)
    potu.2=potu[,colnames(potu) %in% samples]
    fotu.2=fotu[,colnames(fotu) %in% samples]
    
    potu.2=potu.2[rowSums(potu.2 >=1) >= length(samples)*0.50,]
    fotu.2=fotu.2[rowSums(fotu.2 >=1) >= length(samples)*0.20,]
    
    otu=rbind(potu.2,fotu.2)
    tax=rbind(ptax,ftax)
    
    #improve taxonomy
    for (i in 6:1){
      tax[is.na(tax[,i]),i:6]=tax[is.na(tax[,i]),i-1]
    }
    
    sum(rowSums(otu)==0)
    sum(colSums(otu)==0)
    
    met.2=met[met$host.id %in% c("Ac","Hm"),]
    otu.2=otu[,colnames(otu) %in% samples]
    met.3=met.2[rownames(met.2) %in% colnames(otu.2),]
    
    otu.2=otu.2[,colSums(otu.2)>1]
    otu.2=otu.2[rowSums(otu.2)>0,]
    met.3=met.3[colnames(otu.2),]
    
    tax.2=tax[rownames(otu.2),]
    
    
    sparcc.1 <- sparcc(t(otu.2))
    sparcc.1$Cor[sparcc.1$Cor > -0.6 & sparcc.1$Cor < 0.6] <- 0
    sparcc.graph.1 <- sparcc.1$Cor
    diag(sparcc.graph.1) <- 0
    sparcc.graph.1 <- Matrix(sparcc.graph.1, sparse=TRUE)
    g <- graph_from_adjacency_matrix(sparcc.graph.1, mode = "undirected", weighted = T)
    vertex_attr(g)$name=rownames(otu.2)
    g=delete_vertices(g, V(g)[degree(g) == 0])
    vertex_attr(g, name="kingdom", index = V(g)) <- tax[V(g)$name,1]
    vertex_attr(g, name="zymo", index = V(g)) <- ifelse(tax[V(g)$name,6] == "g:Zymoseptoria", "Zymoseptoria","Other")
    edge_attr(g, name="peso",index=E(g)) <- ifelse(E(g)$weight>0, "grey50", "red")
    
    if (ho != "Hordeum.murinum"){
    a=ggnet2(g, node.color="kingdom", edge.color = "peso",
             mode="fruchtermanreingold", edge.size = 0.7, node.size = "zymo")+
      geom_point(aes(fill=color,size=size),shape=21)+
      scale_size_manual(values = c(6,4), limits=c("Zymoseptoria","Other"))+
      scale_fill_manual(values =c(colores), 
                        limits=c("d:Bacteria","d:Archaea","d:Fungi","d:unidentified_eukaryota","d:Viridiplantae"))
    #geom_text_repel(label=tax[V(g)$name,6],size=2)
    png(paste(plot.dir,paste(method,ho,"network","png",sep = "."),sep = "/"), 
        width = 1000, height = 600, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(a) ; dev.off()
    
    z=zotu
    if (sum(z %in% names(V(g)))>0){
  
    #Extract zymo vertices
    z=z[z %in% names(V(g))]
    edges=igraph::as_data_frame(g, what = c("edges"))
    from=edges[edges$from %in% z,"to"]
    to=edges[edges$to %in% z,"from"]
    
    
    g2=induced_subgraph(g,c(from,to,z))
  
      a=ggnet2(g2, node.color="kingdom", edge.color = "peso",
             mode="fruchtermanreingold", edge.size = 0.7, node.size = "zymo")+
      geom_point(aes(fill=color,size=size),shape=21)+
      scale_size_manual(values = c(6,4), limits=c("Zymoseptoria","Other"))+
      scale_fill_manual(values =c(colores), 
                        limits=c("d:Bacteria","d:Archaea","d:Fungi","d:unidentified_eukaryota","d:Viridiplantae"))+
      geom_text_repel(label=paste(tax[V(g2)$name,6],tax[V(g2)$name,7]),size=2)
    
    png(paste(plot.dir,paste(method,ho,"subnetwork","png",sep = "."),sep = "/"), 
        width = 1000, height = 600, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(a) ; dev.off()
    
    
    
    #Remove non zymo edges
    edges2=igraph::as_data_frame(g2, what = c("edges"))
    edges2[edges2$from %in% z,"Ez"]=1
    edges2[edges2$to %in% z,"Ez"]=1
    seq=1:length(edges2$from)
    g3=delete_edges(g2,seq[is.na(edges2$Ez)])
    plot=ggnet2(g3, node.color="kingdom", edge.color = "peso",
                mode="fruchtermanreingold", edge.size = 0.7, node.size = "zymo")+
      geom_point(aes(fill=color,size=size),shape=21)+
      scale_size_manual(values = c(6,4), limits=c("Zymoseptoria","Other"))+
      scale_fill_manual(values =colores, 
                        limits=c("d:Bacteria","d:Archaea","d:Fungi","d:unidentified_eukaryota","d:Viridiplantae"))+
      geom_text_repel(label=paste(tax[V(g2)$name,6],tax[V(g2)$name,7]),size=2)
    
    #png(paste(plot.dir,paste(method,ho,"subnetwork_onlyzymo","png",sep = "."),sep = "/"), 
     #   width = 1000, height = 600, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    #print(plot) ; dev.off()
    }
    }
    
    #Select the most connected OTUS
    grado=degree(g, v = V(g), mode = c("all"),loops = F, normalized = FALSE)
    grado=grado[grado>summary(grado)[5]]
    
    #Select the most central OTUS
    central=betweenness(g, v = V(g), directed = F, weights = abs(E(g)$weight), normalized = FALSE)
    central=central[central>summary(central)[5]]
    
    #Select the hub OTUS
    hub=names(grado[names(grado) %in% names(central)])
    
    #
    
    #Generate the result data.frames
    
    grafica=data.frame(degree=degree(g, v = V(g), mode = c("all"),loops = F, normalized = FALSE),betweenesscentrality=betweenness(g, v = V(g), directed = F, weights = abs(E(g)$weight), normalized = FALSE))
    grafica=cbind(grafica,tax.2[rownames(grafica),]) #add taxonomy
    grafica$hub=0 ; grafica$central=0 ; grafica$grado=0 #add extra columns
    grafica[hub,"hub"]=1 ; grafica[names(central),"central"]=1 ; grafica[names(grado),"grado"]=1 #indicate the hub, highly conncted and cetral OTUs
    
    write.table(grafica, paste(plot.dir,paste(method,ho,"hubs","txt",sep = "."),sep = "/"), quote=F, sep = "\t")
    #write.table(igraph::as_data_frame(g3, what = c("edges")), paste(plot.dir,paste(method,ho,"zymo_edges","txt",sep = "."),sep = "/"), quote=F, sep = "\t")
    
    
    ## RANDOM NETWORKS ##
    
    #Generate the result data.frames
    result=data.frame(row.names = rownames(grafica), 
                      PRS=rep(0, times=dim(grafica)[1]), #vertex present in random sample
                      HRS=rep(0, times=dim(grafica)[1]), #hub present in random sample
                      GRS=rep(0, times=dim(grafica)[1]), #highly connected vertex present in random sample
                      CRS=rep(0, times=dim(grafica)[1]), #highly central vertex present in random sample
                      NRS=rep(0, times=dim(grafica)[1])) #vertex non-hub but present in random sample
    
    z.edges=data.frame()
    
    #Set seeds (every loop will subsample with a different seed)
    initial.seed=1587571314
    the.seed=71314
    
    for (seed in seq(the.seed, initial.seed, 15880000 )){ 
      
      print((seed))
      random=data.frame(row.names = rownames(grafica), 
                        PRS=rep(0, times=dim(grafica)[1]),
                        HRS=rep(0, times=dim(grafica)[1]),
                        GRS=rep(0, times=dim(grafica)[1]),
                        CRS=rep(0, times=dim(grafica)[1]),
                        NRS=rep(0, times=dim(grafica)[1]))
    
      set.seed(seed)
      r.samples=sample(x=samples, size=10)
      
      potu.2=potu[,colnames(potu) %in% r.samples]
      fotu.2=fotu[,colnames(fotu) %in% r.samples]
      
      set.seed(seed)
      potu.2=potu.2[sample(x=rownames(potu.2),size= dim(potu.2)[1]/2),]
      fotu.2=fotu.2[sample(x=rownames(fotu.2),size= dim(fotu.2)[1]/2),]
      
      potu.2=potu.2[rowSums(potu.2 >=1) >= length(r.samples)*0.50,]
      fotu.2=fotu.2[rowSums(fotu.2 >=1) >= length(r.samples)*0.20,]
      
      rotu=rbind(potu.2,fotu.2)
      tax=rbind(ptax,ftax)
      
      #improve taxonomy
      for (i in 6:1){
        tax[is.na(tax[,i]),i:6]=tax[is.na(tax[,i]),i-1]
      }
      
      
      rotu.2=rotu[,colnames(rotu) %in% r.samples]
      met.3=met.2[rownames(met.2) %in% colnames(rotu.2),]
      
      rotu.2=rotu.2[,colSums(rotu.2)>1]
      rotu.2=rotu.2[rowSums(rotu.2)>0,]
      met.3=met.3[colnames(rotu.2),]
      
      tax.2=tax[rownames(rotu.2),]
      
      
      sparcc.1 <- sparcc(t(rotu.2))
      sparcc.1$Cor[sparcc.1$Cor > -0.6 & sparcc.1$Cor < 0.6] <- 0
      sparcc.graph.1 <- sparcc.1$Cor
      diag(sparcc.graph.1) <- 0
      sparcc.graph.1 <- Matrix(sparcc.graph.1, sparse=TRUE)
      r.g <- graph_from_adjacency_matrix(sparcc.graph.1, mode = "undirected", weighted = T)
      vertex_attr(r.g)$name=rownames(rotu.2)
      r.g=delete_vertices(r.g, V(r.g)[degree(r.g) == 0])
      
      #Select the most connected OTUS
      r.grado=degree(r.g, v = V(r.g), mode = c("all"),loops = F, normalized = FALSE)
      r.grado=r.grado[r.grado>summary(r.grado)[5]]
      
      #Select the most central OTUS
      r.central=betweenness(r.g, v = V(r.g), directed = F, weights = abs(E(r.g)$weight), normalized = FALSE)
      r.central=r.central[r.central>summary(r.central)[5]]
      
      #Select the hub OTUS
      r.hub=names(r.grado[names(r.grado) %in% names(r.central)])
      
      #Determine the presence of the vertices
      #PRS
      random[rownames(random) %in% rownames(rotu.2),"PRS"]=1
      
      #HRS
      random[rownames(random) %in% r.hub,"HRS"]=1
      #GRS
      random[rownames(random) %in% names(r.grado),"GRS"]=1
      #CRS
      random[rownames(random) %in% names(r.central),"CRS"]=1
      #NRS
      random[random$PRS==1 & random$HRS==0, "NRS"]=1
      
      #Generate the counts 
      result=result+random
      
      #Calculate metrics
      dens=c(dens,edge_density(r.g, loops = FALSE)) 
      diam=c(diam,diameter(r.g, directed=F, weights=NA))
      modu=c(modu,modularity(r.g, membership(cluster_fast_greedy(r.g,weights = abs(E(r.g)$weight))), weights = abs(E(r.g)$weight)))
      tran=c(tran,transitivity(r.g, type = "undirected", vids = NULL, weights = abs(E(r.g)$weight), isolates = c("NaN")))
      mdeg=c(mdeg,mean(degree(r.g, v = V(r.g), mode = c("all"),loops = F, normalized = FALSE)))
      #aspa=c(aspa,mean_distance(r.g, directed = TRUE, unconnected = TRUE))
      node=c(node,length(V(r.g)))
      edge=c(edge,length(E(r.g)))
      n.hub=c(n.hub,length(r.hub))
      
      #Determine bacterial-fungal edges
      ver=igraph::as_data_frame(r.g, what = c("edges"))
      one=ver[grep("POTU_",ver$from),]
      two=dim(one[grep("FOTU_",one$to),])[1]
      three=ver[grep("FOTU_",ver$from),]
      four=dim(three[grep("POTU_",three$to),])[1]
      
      bf=c(bf,c(two+four))
      bfp=c(bfp,c((two+four)/dim(ver)[1]*100))
      
      #Determine bacterial-bacteria edges
      one=ver[grep("POTU_",ver$from),]
      two=dim(one[grep("POTU_",one$to),])[1]
      bb=c(bb,c(two))
      bbp=c(bbp,c((two)/dim(ver)[1]*100))
      
      #Determine fungal-fungal edges
      one=ver[grep("FOTU_",ver$from),]
      two=dim(one[grep("FOTU_",one$to),])[1]
      ff=c(ff,c(two))
      ffp=c(ffp,c((two)/dim(ver)[1]*100))
      
      #Determine zymo negative and positive interactions
      #print(sum(z %in% names(V(r.g)))>0)
      z=zotu
      if (sum(z %in% row.names(rotu.2))>0 & sum(z %in% names(V(r.g)))>0){
        
        #Extract zymo vertices
        z=z[z %in% names(V(r.g))]
        edges=igraph::as_data_frame(r.g, what = c("edges"))
       
        #rEMOVE OTHER zYMO CORRELATIONS
        edges=edges[!edges$from %in% other,]
        edges=edges[!edges$to %in% other,]
        #
        
        
        from=edges[edges$from %in% z,"to"]
        to=edges[edges$to %in% z,"from"]
      
        
        r.g2=induced_subgraph(r.g,c(from,to,z))
        
        #Remove non zymo edges
        edges2=igraph::as_data_frame(r.g2, what = c("edges"))
        edges2= rbind(edges2[edges2$from %in% z,], edges2[edges2$to %in% z,])
        #all=100/dim(edges2)[1]
        all=100/length(E(r.g))
        #all=1/mean(degree(r.g, v = V(r.g), mode = c("all"),loops = F, normalized = FALSE))
        z.all=data.frame(total=dim(edges2)[1],negative=sum(edges2$weight<0),positive=sum(edges2$weight>0),
                         total.p=dim(edges2)[1]*all, negative.p=sum(edges2$weight<0)*all, positive.p=sum(edges2$weight>0)*all,
                         host=ho)
        
        print("present")
        
      } else if ((sum(z %in% row.names(rotu.2))>0 & sum(z %in% names(V(r.g)))==0)) {
        
        z.all=data.frame(total=0,negative=0,positive=0,
                         total.p=0, negative.p=0, positive.p=0,
                         host=ho)
        
      } else {
        
        
        z.all=c()
        
      }
      
      z.edges=rbind(z.edges,z.all)
      
    }
    
    #Calculate frequencies
    result$HF=round(result$HRS/result$PRS*100,digits=1)
    result$GF=round(result$GRS/result$PRS*100,digits=1)
    result$CF=round(result$CRS/result$PRS*100,digits=1)
    
    #Merge dataframes
    grafica$random.hub=grafica$hub
    grafica$random.grado=grafica$grado
    grafica$random.central=grafica$central
    sum(rownames(grafica)==rownames(result))==dim(grafica)[1]
    
    #Indicate if a vertex is a hub, central or connected based on their frequencies
    grafica[result$HF<50,"random.hub"]=0
    grafica[result$GF<50,"random.grado"]=0
    grafica[result$CF<50,"random.central"]=0
    
    #Save results
    
    #Each vertex with: centrality, degree, taxonomy, hub, central, connected (both in a single network and across random networks)
    write.table(grafica, file = paste(plot.dir,paste(ho,"R_nodes","txt",sep = "."),sep = "/"),
                row.names = F, quote = F, sep = "\t" ,col.names = T)
    
    #hub, connected and central frequencies in random networks
    write.table(result, file = paste(plot.dir,paste(ho,"random_networks","txt",sep = "."),sep = "/"),
                row.names = T, quote = F, sep = "\t", col.names = T)
    
    #plot
    
    
    
    #Generate plots
    freq=result
    freq$OBH="NO HUB"
    freq[hub,"OBH"]="HUB"
    freq[,c("HRS","NRS")]=freq[,c("HRS","NRS")]/freq[,c("PRS")]*100
    
    a=ggplot(data=freq, aes(x=HRS,  fill=OBH))+
      geom_histogram(alpha=.50,position="identity", size=0.9)+
      scale_fill_manual(values=colores.h,limits=c("NO HUB","HUB"))+
      ggtitle("Hub frequency in random networks (all otus)")+
      xlab("frequency of an OTU being a hub")+ylab("OTU Count")+
      geom_vline(xintercept = 50)
    
    png(paste(plot.dir,paste(method,ho,"hubfreq","png",sep = "."),sep = "/"), 
        width = 1000, height = 600, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(a) ; dev.off()
    
    
    metrics=data.frame(host=rep(ho,times=length(seq(the.seed, initial.seed, 15880000 ))),
                       density=dens,diameter=diam,modularity=modu,clustering=tran,
                       av.degree=mdeg,vertex=node,edges=edge,
                       p.e.edges=bfp, no.hubs=n.hub,p.p.edges=bbp,e.e.edges=ffp)
    
    result.metrics=rbind(result.metrics,metrics)
    result.zymo=rbind(result.zymo,z.edges)
    
    
  ##########}
  
  write.table(result.metrics, file = paste(plot.dir,paste(ho,"metrics","txt",sep = "."),sep = "/"),
              row.names = F, quote = F, sep = "\t" ,col.names = T)
  
  
  #result.metrics$treatment=factor(result.metrics$treatment, levels=treat)
  
  a1=ggplot(result.metrics,aes(x=host,y=diameter))+geom_boxplot()
  a2=ggplot(result.metrics,aes(x=host,y=modularity))+geom_boxplot()
  a3=ggplot(result.metrics,aes(x=host,y=clustering))+geom_boxplot()
  a4=ggplot(result.metrics,aes(x=host,y=av.degree))+geom_boxplot()
  a5=ggplot(result.metrics,aes(x=host,y=vertex))+geom_boxplot()
  a6=ggplot(result.metrics,aes(x=host,y=edges))+geom_boxplot()
  a7=ggplot(result.metrics,aes(x=host,y=no.hubs))+geom_boxplot()
  a8=ggplot(result.metrics,aes(x=host,y=p.e.edges))+geom_boxplot()
  a9=ggplot(result.metrics,aes(x=host,y=p.p.edges))+geom_boxplot()
  a10=ggplot(result.metrics,aes(x=host,y=e.e.edges))+geom_boxplot()
  
  a=plot_grid(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
  png(paste(plot.dir,paste(ho,"metrics","png",sep = "."),sep = "/"), 
      width = 2000, height = 1500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(a) ; dev.off()
  
  
  
  #result.zymo$treatment=factor(result.zymo$treatment, levels=treat)
  a1=ggplot(result.zymo,aes(x=host,y=total.p))+geom_boxplot()
  a2=ggplot(result.zymo,aes(x=host,y=negative.p))+geom_boxplot()
  a3=ggplot(result.zymo,aes(x=host,y=positive.p))+geom_boxplot()
  a=plot_grid(a1,a2,a3)
  
  write.table(z.edges, file = paste(plot.dir,paste(ho,"random_zymo_edges","txt",sep = "."),sep = "/"),
              row.names = T, quote = F, sep = "\t", col.names = T)
  
  png(paste(plot.dir,paste(ho,"zymo_total","png",sep = "."),sep = "/"), 
      width = 2000, height = 1500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(a) ; dev.off()
  
  
  
}


