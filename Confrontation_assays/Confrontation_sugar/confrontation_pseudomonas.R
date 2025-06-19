#Generates the plot for figure 6. The effect of sugar in the confrontation test and invertase activity

setwd("C:/Users/victorfn/Downloads/")

library(ggplot2)
library(FSA)



colores.h=c("#FFC850","#116875","#FF4848","#41c6a8","black")
tema=theme(axis.text.x = element_text(color="black",size=9, angle=0,hjust=0.5,vjust=0.5,family = "Arial"),
           axis.text.y = element_text(color="black",size=9,family = "Arial"),
           axis.title = element_text(color="black",size=9,family = "Arial"),
           legend.text = element_text(color = "black",size=9,family = "Arial"),
           legend.key.size = unit(0.3,"cm"),
           legend.title = element_text(color="black",size=9,family = "Arial"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=9, color="black",family = "Arial", margin = margin(0.1,0.1,0.1,0.1,"cm")),
           strip.text.y = element_text(size=9, color="black",family = "Arial"),
           strip.background = element_rect(fill="gray70"),
           panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "solid"),
           panel.grid.major.y = element_line(color = "gray70"),panel.grid.minor.x = element_blank(),
           legend.position = "right")


#Load data
data=read.delim("Results_carbon_tidy.txt")[1:6]
data$pathogen=factor(data$pathogen, levels = c("Zt469", "Zt549", "Zpa796","Zpa21" ))
data$od=factor(data$od, levels=c("2.7","1.35"))

#sanity review
ggplot(data=data, aes(x=pathogen, y=halo, fill=od))+
  geom_boxplot()+
  facet_wrap(~bacteria, nrow = 1)

ggplot(data=data, aes(x=pathogen, y=colony, fill=od))+
  geom_boxplot()+
  facet_wrap(~bacteria, nrow = 1)

ggplot(data=data, aes(x=pathogen, y=ttc, fill=od))+
  geom_boxplot()+
  facet_wrap(~bacteria, nrow = 1)

#Calculation of ratios
data[data$halo==0,"halo"]=data[data$halo==0,"colony"]
data[data$ttc==0,"ttc"]=data[data$ttc==0,"colony"]
data$ratio_pro=data$halo/data$colony
data$ratio_ttc=data$ttc/data$colony

#Growth enhancement
a=ggplot(data=data, aes(x=pathogen, y=ratio_pro, fill=od))+
  geom_boxplot(size=0.7)+scale_fill_manual(values = colores.h[c(1,2)])+
  facet_wrap(~bacteria, nrow = 1)+ylab("Ratio_enhancement")+tema
a

png("enhancement.png", 
    width = 3000, height = 600, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

#Invertase activity by TTC
a=ggplot(data=data, aes(x=pathogen, y=ratio_ttc, fill=od))+
  geom_boxplot(size=0.7)+scale_fill_manual(values = colores.h[c(1,2)])+
  facet_wrap(~bacteria, nrow = 1)+ylab("Ratio_invertase")+tema
a
png("invertase.png", 
    width = 3000, height = 600, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

#Statistical test

#Growth enhancement #between zymoseptoria lineages
for (i in c("ac26", "ac30","ac66", "ac9")){
  sub=data[data$bacteria==i,]
  #kruskal.test(ratio_ttc~pathogen, data = sub)
  print(i)
  print(dunnTest(ratio_pro~pathogen, data = sub)$res)
  
}

#TTC#between zymoseptoria lineages
for (i in c("ac26", "ac30","ac66", "No-bacteria")){
  sub=data[data$bacteria==i,]
  #kruskal.test(ratio_ttc~pathogen, data = sub)
  print(i)
  print(dunnTest(ratio_ttc~pathogen, data = sub)$res)
  
}


#in Zt469 betwween OD concentrations

for (i in c("ac26", "ac30","ac66","ac9", "No-bacteria")){
  
  print(i)
  sub=data[data$bacteria==i & data$pathogen=="Zt469",]
  print(wilcox.test(sub[sub$od==2.7,"ratio_pro"], sub[sub$od==1.35,"ratio_pro"]))
  print(wilcox.test(sub[sub$od==2.7,"ratio_ttc"], sub[sub$od==1.35,"ratio_ttc"]))
  
  
}


#in Zt549 betwween OD concentrations
for (i in c("ac66","ac9", "No-bacteria")){
  
  print(i)
  sub=data[data$bacteria==i & data$pathogen=="Zt549",]
  print(wilcox.test(sub[sub$od==2.7,"ratio_pro"], sub[sub$od==1.35,"ratio_pro"]))
  print(wilcox.test(sub[sub$od==2.7,"ratio_ttc"], sub[sub$od==1.35,"ratio_ttc"]))
  
  
}

#in Zpa796 betwween OD concentrations
for (i in c("ac9", "No-bacteria")){
  
  print(i)
  sub=data[data$bacteria==i & data$pathogen=="Zpa796",]
  print(wilcox.test(sub[sub$od==2.7,"ratio_pro"], sub[sub$od==1.35,"ratio_pro"]))
  print(wilcox.test(sub[sub$od==2.7,"ratio_ttc"], sub[sub$od==1.35,"ratio_ttc"]))
  
  
}















