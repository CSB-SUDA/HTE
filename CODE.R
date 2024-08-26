###part1 8组火山图####
#G:\P1\1.8组火山图

library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(Seurat)
library(CellChat)


load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz) <- "RNA"
Rz


##R110
diff_R110_Ctrl <- FindMarkers(Rz,ident.1=c("R1-10"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R110_Ctrl)
write.csv(diff_R110_Ctrl,file="R110_Ctrl_DEG_UP.csv")

diff_R110_Ctrl <- FindMarkers(Rz,ident.1=c("R1-10"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R110_Ctrl,file="R110_Ctrl_DEG_UP.csv")


###R1

diff_R1_Ctrl <- FindMarkers(Rz,ident.1=c("R1-10","R1-30","R1-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R1_Ctrl)
write.csv(diff_R1_Ctrl,file="R1_Ctrl_DEG_UP.csv")


diff_R1_Ctrl <- FindMarkers(Rz,ident.1=c("R1-10","R1-30","R1-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R1_Ctrl,file="R1_Ctrl_DEG_UP.csv")
Rz

##R130
diff_R130_Ctrl <- FindMarkers(Rz,ident.1=c("R1-30"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R130_Ctrl)
write.csv(diff_R130_Ctrl,file="R130_Ctrl_DEG_UP.csv") 

diff_R130_Ctrl <- FindMarkers(Rz,ident.1=c("R1-30"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R130_Ctrl,file="/public/home/SZUYanwy/LWS0726/R130_Ctrl_DEG_UP.csv")
Rz
##R150
diff_R150_Ctrl <- FindMarkers(Rz,ident.1=c("R1-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R150_Ctrl)
write.csv(diff_R150_Ctrl,file="R150_Ctrl_DEG_UP.csv") 

diff_R150_Ctrl <- FindMarkers(Rz,ident.1=c("R1-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R150_Ctrl,file="/public/home/SZUYanwy/LWS0726/R150_Ctrl_DEG_UP.csv")
Rz



##R2510
diff_R2510_Ctrl <- FindMarkers(Rz,ident.1=c("R25-10"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R2510_Ctrl)
write.csv(diff_R2510_Ctrl,file="R2510_Ctrl_DEG_UP.csv")


diff_R2510_Ctrl <- FindMarkers(Rz,ident.1=c("R25-10"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R2510_Ctrl,file="/public/home/SZUYanwy/LWS0726/R2510_Ctrl_DEG_UP.csv")
Rz


##R2530
diff_R2530_Ctrl <- FindMarkers(Rz,ident.1=c("R25-30"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R2530_Ctrl)
write.csv(diff_R2530_Ctrl,file="R2530_Ctrl_DEG_UP.csv")


diff_R2530_Ctrl <- FindMarkers(Rz,ident.1=c("R25-30"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R2530_Ctrl,file="/public/home/SZUYanwy/LWS0726/R2530_Ctrl_DEG_UP.csv")
Rz

##R2550
diff_R2550_Ctrl <- FindMarkers(Rz,ident.1=c("R25-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R2550_Ctrl)
write.csv(diff_R2550_Ctrl,file="R2550_Ctrl_DEG_UP.csv")

diff_R2550_Ctrl <- FindMarkers(Rz,ident.1=c("R25-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R2550_Ctrl,file="/public/home/SZUYanwy/LWS0726/R2550_Ctrl_DEG_UP.csv")
Rz

##R25
diff_R25_Ctrl <- FindMarkers(Rz,ident.1=c("R25-10","R25-30","R25-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0.1,assay = "RNA",slot = "counts",logfc.threshold = 0.25,only.pos = TRUE,test.use = "bimod")
dim(diff_R25_Ctrl)
write.csv(diff_R25_Ctrl,file="R25_Ctrl_DEG_UP.csv")


diff_R25_Ctrl <- FindMarkers(Rz,ident.1=c("R25-10","R25-30","R25-50"),group.by = "orig.ident",ident.2="Ctrl",min.pct = 0,assay = "RNA",slot = "counts",logfc.threshold = 0,only.pos = FALSE,test.use = "bimod")
write.csv(diff_R25_Ctrl,file="R25_Ctrl_DEG_UP.csv")
Rz
#huahuoshantu

setwd("G:\\yan\\23.8.7\\差异新的")


object.markers <- read.csv("R25_HUOSHAN.CSV",row.names = 1)
object.markers$names <- rownames(object.markers)
library(dplyr)
object.markers <- object.markers %>%
  mutate(Difference = pct.1 - pct.2)
library(ggplot2)
library(ggrepel)

ggplot(object.markers, aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=0.5, color="#999999") + 
  theme_classic()




object.markers$group=0

for (i in 1:nrow(object.markers)){
  if (object.markers$avg_log2FC[i] >= 1 & object.markers$Difference[i] >= 0.2){
    object.markers$group[i]='up'
  }
  else if(object.markers$avg_log2FC[i] <= -1 & object.markers$Difference[i] <= -0.2 ){
    object.markers$group[i]='down'
  }
  else {
    object.markers$group[i]='no'
  }
}
write.csv(object.markers,file="R25_hutuqian.csv")

setwd("G:\\P1\\1.8组火山图")
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
object.markers <- read.csv("R10_hutuqian.csv",row.names = 1) 

p <-  ggplot(object.markers, aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=3,aes(color=group)) + 
  scale_color_manual(values=c('blue','grey','red'))+
  #  geom_label_repel(data=subset(object.markers, group !='no'), aes(label=names), segment.size = 0.25, size=2.5)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()
p
library(ggplot2)

object.markers$gene <- rownames(object.markers)
object.markers$label=ifelse(object.markers$marker == 1, as.character(object.markers$gene), '')
p <- p+
  theme(panel.grid=element_blank()) +
  geom_text_repel(aes(x = object.markers$Difference,                   # geom_text_repel 标记函数
                      y = object.markers$avg_log2FC,          
                      label=object.markers$label),                       
                  max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
                  size=6,family = "B", colour="black",                            # 字体大小
                  box.padding=unit(0.5,'lines'),           # 标记的边距
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                   # 标记线条的颜色
                  show.legend=FALSE)+
  theme(
    axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(size = 20,family = 'B'),
        legend.title = element_text(size = 20,family = 'B'))

p1 <- p+
  theme(panel.grid=element_blank()) +
  theme(
    axis.text.y=element_text(colour="black",size=12,family = 'B'),)+
  theme(axis.text.x = element_text(colour="black",family = "B",size=12),)+
  labs(x=NULL,y=NULL) + 
  theme(axis.title.x = element_text(colour="black",size=12,family = "B"),)+
  theme(axis.title.y = element_text(colour="black",size=12,family = "B"),)+
  theme(legend.text = element_text(size = 20,family = 'B'),
        legend.title = element_text(size = 20,family = 'B'))


p1



library(eoffice)
topptx(p,"R150HUOSHAN_2024.pptx",width = 7.81,height = 6.56)  

###p1 venn####
#G:\P1\1.8组火山图\venn
setwd("G:\\yan\\23.8.7\\差异新的\\venn")
getwd()


###p1 8zudekegg####

setwd("G:\\yan\\23.8.7\\差异_照搬10.30的\\大圆")
#G:\P1\1.8组气泡图KEGG\大圆
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)
##kegg
diff<-rownames(diff_R110_Ctrl) 
gene<-bitr(t(diff),fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
gene <- gene$ENTREZID
#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 1, pvalueCutoff=1)

hh <- as.data.frame(kk)#自己记得保存结果哈！
write.csv(hh,file="R2550_keggres.csv",row.names = F)

##qipaotu
library(ggplot2)
hh <- hh[1:36,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
p <- ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="",title="KEGG Pathway Enrichment")+
  theme_bw()+
  theme(text=element_text(face="bold",family = 'A'))
p+ theme(
  axis.text.y=element_text(colour="black",size=18,face="bold",family = 'A'),)+
  theme(axis.text.x = element_text(size=15,face="bold",family = "A"),)+
  theme(legend.text = element_text(size = 20,face = 'bold',family = 'A'),
        legend.title = element_text(size = 20,face = 'bold',family = 'A'))
#现在算和之前不完成一样，以文件夹的为准，气泡图出了后用keggres去公众号画分类的那个图
###改kegg的dotplot####
setwd("G:\\P1/1.8组气泡图KEGG/xindediantu/")
a <- read.csv("工作簿1.csv")
c <- unique(a$ID)

a$clu <-  factor(a$clu,levels = c("R1-10","R1-30","R1-50","R25-10","R25-30","R25-50"))
a$order=factor(rev(as.integer(rownames(a))),labels = rev(a$Description))

p12 <- ggplot(a,aes(y=order,x=clu))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12

library(eoffice)
topptx(p12,"KEGG的通路图-yangben.pptx",width = 7,height = 6)





###2.单细胞尺度异质性####



load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz) <- "RNA"
Rz

newname <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12")
names(newname) <- levels(Rz)
Rz1 <- RenameIdents(Rz,newname)
DimPlot(Rz1,reduction = "umap",ncol=2,pt.size = 0.1,label = T)
save(Rz1,file="G:\\yan\\22.11.14\\Rz_0.5.RData")
table(Rz1$sample)

Rz1$xin <- Rz1@active.ident
write.csv(Rz1$sample,file="G:\\P2\\P2.1UMAP图\\yangben_gai.csv")
b <- Rz1$sample
c <- read.csv("G:\\P2\\P2.1UMAP图\\yangben_gai.csv",row.names = 1)
unique(Rz1$sample)
Rz1$xinde <- c
Idents(Rz1) <- "sample"
table(Rz1$xinde)
Rz1$xinde <- factor(Rz1$xinde, levels = c("Ctrl","Ctrl1","Ctrl2","R1-10","R1-30","R1-50","R25-10","R25-30","R25-50"))
cols <- c("#D3D3D3","white","white","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
names(cols) <-c("Ctrl","Ctrl1","Ctrl2","R1-10","R1-30","R1-50","R25-10","R25-30","R25-50")
p1 <- DimPlot(Rz1,reduction = "umap",ncol=3,pt.size = 0.1,label = F,cols = cols,split.by = "xinde")
p1
library(ggplot2)
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
setwd("G:\\P2\\P2.1UMAP图")
topptx(p2,"SAMPLE.pptx",width =8.89,height = 6.82)


##平均蓝：#1F78B4  平均红：#ED5F5F
##anzhaocluster
unique(Rz1$xin)
Idents(Rz1) <- "xin"

##1新特性
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
##2商务
cols <- c("194f97","#555555","#bd6b08","#00686b","#c82d31","#625ba1","#898989","#9c9800","#007f54","#a195c5","#103667","#f19272")

##3雅致
cols <- c("#3682be","#45a776","#f05326","#eed777","#334f65","#b3974e","#38cb7d","#ddae33","#844bb3","#93c555","#5f6694","#df3881")
##4复古
cols <- c("#0780cf","#765005","#fa6d1d","#0e2c82","#b6b51f","#da1f18","#701866","#f47a75","#009db2","#024b51","#0780cf","#765005")
###xuanze1
names(cols) <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12")
p1 <- DimPlot(Rz1,reduction = "umap",ncol=2,pt.size = 0.1,label = F,cols = cols)
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2

library(eoffice)
getwd()
topptx(p2,"group.pptx",width =4.85,height = 3.5)

###P2.2cluster比例变化图####
setwd("G:\\P2\\P2.2BILITU")
library(ggplot2)
##huitu

a <- read.csv("huitu_bili1-4.csv")
allcolour=c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
a$Cluster <- factor(a$Cluster, levels = c("cluster 1","cluster 2","cluster 3","cluster 4"))
p1 <- ggplot(a,aes(Sample,Rate))+geom_bar(aes(fill=Sample),stat='identity', position='dodge')+facet_wrap(~Cluster,nrow = 1)+theme_classic()+scale_fill_manual(values=allcolour)+theme(axis.text.x = element_text(angle=90,
                                                                                                                                                                                                  vjust=1,hjust=1,colour="black",face="bold",family = "A"))                                                                       
P2 <- p1 +  theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))


library(eoffice)
getwd()
topptx(P2,"huitu_bili1-4.pptx",width =11.18,height = 3.84)
#5-8
setwd("G:\\P2\\P2.2BILITU")
library(ggplot2)
##huitu

a <- read.csv("huitu_bili5-8.csv")
allcolour=c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
a$Cluster <- factor(a$Cluster, levels = c("cluster 5","cluster 6","cluster 7","cluster 8"))
p1 <- ggplot(a,aes(Sample,Rate))+geom_bar(aes(fill=Sample),stat='identity', position='dodge')+facet_wrap(~Cluster,nrow = 1)+theme_classic()+scale_fill_manual(values=allcolour)+theme(axis.text.x = element_text(angle=90,
                                                                                                                                                                                                                 vjust=1,hjust=1,colour="black",face="bold",family = "A"))                                                                       
P2 <- p1 +  theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))

P2
library(eoffice)
getwd()
topptx(P2,"huitu_bili5-8.pptx",width =11.18,height = 3.84)


###9-12


setwd("G:\\P2\\P2.2BILITU")
library(ggplot2)
##huitu

a <- read.csv("huitu_bili9-12.csv")
allcolour=c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
a$Cluster <- factor(a$Cluster, levels = c("cluster 9","cluster 10","cluster 11","cluster 12"))
p1 <- ggplot(a,aes(Sample,Rate))+geom_bar(aes(fill=Sample),stat='identity', position='dodge')+facet_wrap(~Cluster,nrow = 1)+theme_classic()+scale_fill_manual(values=allcolour)+theme(axis.text.x = element_text(angle=90,
                                                                                                                                                                                                                 vjust=1,hjust=1,colour="black",face="bold",family = "A"))                                                                       
P2 <- p1 +  theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))

P2
library(eoffice)
getwd()
topptx(P2,"huitu_bili9-12.pptx",width =11.18,height = 3.84)

###P2.3整体比例图####
##huitu
#y样本xcluster
setwd("G:\\P2\\P3.3比例整合图")
a <- read.csv("huitu_bili.csv")


a <- read.csv("huitu_bili.csv")
allcolour=c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
a$Sample <-  factor(a$Sample,levels = c("Ctrl","R1-10","R1-30","R1-50","R25-10","R25-30","R25-50"))



p1 <- ggplot(a) + 
  geom_bar(aes(x =Cluster, y= Rate, fill =Sample),stat = "identity",position = "fill",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black",size=0.5, linetype="solid"))+
  scale_x_discrete(limits = c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12"))+
  #scale_y_discrete(limits = c("Ctrl","R1-10","R1-30","R1-50","R25-10","R25-30","R25-50"))+
  RotatedAxis()
 

P2 <- p1 + theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(axis.text.x = element_text(angle=90,
                                                              vjust=0.3,hjust=1,colour="black",face="bold",family = "Times New Roman"))+
  theme(legend.text = element_text(colour="black",face="bold",family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
P2

library(eoffice)
topptx(P2,"P2.3整体比例图.pptx",width =6.6,height = 6.12)
###p2.4多重火山图####


####meihuatu

library(scRNAtoolVis)
setwd("G:\\P2\\P2.4多组火山图")

a <- read.csv("差异的.csv")

a$cluster <-  factor(a$cluster,levels = c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12"))
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
library(Seurat)
p <- jjVolcano(diffData = a,
          tile.col = cols,
          size  = 3.5, 
          fontface = 'italic', family = "Times New Roman",
          flip = T,show.legend = F)+NoLegend()
p
                                                                   
p <- p + 
  theme(axis.title.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "Times New Roman"),)

p
library(eoffice)
topptx(p,"P2.4多组火山.pptx",width = 19.19,height =10.16)  
#获取一下颜色和大小

#2.6
setwd("G:\\P2\\P2.6dotplot")
load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz1) <- "RNA"
Rz1




top5 <- read.csv("G:\\P2\\P2.4多组火山图\\差异的.csv")

top51 <-unique(top5$gene)
p1 <- DotPlot(Rz1, features =top51, cols = c("white","red"), dot.scale = 8, 
) + RotatedAxis()
p1


p <- DotPlot(Rz1, features = top51) +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "Times New Roman"))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
p
p3 <- p+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))
p3
library(eoffice)
topptx(p,"G:\\P2\\P2.4多组火山图\\为了获取点P进另一个图.pptx",width = 18.00,height = 6.35)



### P2.5所有top5的marker小提琴图####


load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz1) <- "RNA"
Rz1



##clu1
library(MySeuratWrappers)
markers <- c("H4C3","H2AC20","H1-5","CXCL8","H1-4")
my36colors <-c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91") 
library(ggplot2)
p1 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p1 <- p1+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p1




##clu2

markers <- c("NR2F2","DKK1","GINS2","NEAT1","UPP1")
p2 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")

p2 <- p2+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p2

##clu3

markers <- c("IFI27","IFITM3","ISG15","IFI6","IFITM2")
p3 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")

p3 <- p3+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p3





##clu4

markers <- c("EIF4EBP1","RPLP0","GAS5","RSRC2","RWDD1")
p4 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")



p4 <- p4+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p4

##clu5
markers <- c("CCNB1","CDC20","PTTG1","HMMR","UBE2S")
p5 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")




p5 <- p5+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p5

##clu6
markers <- c("UBE2C","ARL6IP1","AURKA","CKS2","CCNB1")
p6 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")

p6 <- p6+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p6



##clu7
markers <- c("IFIT2","IFIT3","OASL","CCL5","ISG15")
p7 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p7 <- p7+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p7



##clu8
markers <- c("CXCL8","CALR","H1-4","H1-5","H1-2")
p8 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p8 <- p8+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p8


##clu9
markers <- c("PRDX1","C1orf56","CKS2","CTNNB1","TIMP1")
p9 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p9 <- p9+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p9


##clu10
markers <- c("CDKN1A","MDM2","TP53I3","FDXR","GDF15")
p10 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p10 <- p10+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p10

##clu11
markers <- c("MALAT1","MT-CO1","MT-CO2","MT-ND4","NEAT1")
p11 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p11 <- p11+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p11

##clu12
markers <- c("DDIT3","PHLDA2","HSPA1A","SNHG12","HMOX1")
p12 <- VlnPlot(Rz1, features = markers, stacked=T, pt.size= 0, cols = my36colors, direction = "vertical")
p12 <- p12+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))

p12

setwd("G:\\P2\\P2.5每组的top5marker_小提琴图")

library(patchwork)
p44 <- p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12

library(eoffice)
topptx(p44,"P2.5每组的top5marker_小提琴图.pptx",width = 17.19,height = 12) 
#调整字体为新罗马就好了再PPT里


###P2.6dotp[lot美化####

setwd("G:\\P2\\P2.6dotplot")
load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz1) <- "RNA"
Rz1



top5 <- c("H4C3","H2AC20","H1-5","CXCL8","H1-4",
          "NR2F2","DKK1","GINS2","NEAT1","UPP1",
          "IFI27","IFITM3","ISG15","IFI6","IFITM2",
          "EIF4EBP1","RPLP0","GAS5","RSRC2","RWDD1",
          "CCNB1","CDC20","PTTG1","HMMR","UBE2S",
          "UBE2C","ARL6IP1","AURKA","CKS2","CCNB1",
          "IFIT2","IFIT3","OASL","CCL5","ISG15",
          "CXCL8","CALR","H1-4","H1-5","H1-2",
          "PRDX1","C1orf56","CKS2","CTNNB1","TIMP1",
          "CDKN1A","MDM2","TP53I3","FDXR","GDF15",
          "MALAT1","MT-CO1","MT-CO2","MT-ND4","NEAT1",
          "DDIT3","PHLDA2","HSPA1A","SNHG12","HMOX1")

top51 <-unique(top5)
p1 <- DotPlot(Rz1, features =top51, cols = c("white","red"), dot.scale = 8, 
) + RotatedAxis()
p1


p <- DotPlot(Rz1, features = top51) +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "Times New Roman"))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
p
p3 <- p+theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'A'),)+
  theme(axis.title.x = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.title.y = element_text(colour="black",face="bold",family = "A"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black",face="bold",family = "A"),)+
  theme(legend.text = element_text(colour="black",face="bold",family = 'A'),
        legend.title = element_text(colour="black",face = 'bold',family = 'A'))
p3
library(eoffice)
topptx(p,"G:\\P2\\P2.6dotplot\\P2.6dotplot.pptx",width = 12.8,height = 6.35)
#后面根据基迪奥 微信收藏，加了泼墨的图

###P2.7每组marker的kegg####
setwd("G:\\P2\\P2.7每组marker的kegg")


##qipaotu
#cluster1
library(ggplot2)
hh <- read.csv(file="0_jieguo.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p1 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 1")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p1 <- p1+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p1

library(eoffice)
topptx(p1,"CLU0-CLU1.pptx",width = 6,height = 5)
#cluster2
library(ggplot2)
hh <- read.csv(file="1_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p2 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 2")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p2 <- p2+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p2

library(eoffice)
topptx(p2,"CLU2.pptx",width = 6,height = 5)
#cluster3
library(ggplot2)
hh <- read.csv(file="2_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p3 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 3")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p3 <- p3+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p3
library(eoffice)
topptx(p3,"CLU3.pptx",width = 6,height = 5)

#cluster4
library(ggplot2)
hh <- read.csv(file="3_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p4 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 4")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p4 <- p4+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p4
library(eoffice)
topptx(p4,"CLU4.pptx",width = 6,height = 5)
#cluster5
library(ggplot2)
hh <- read.csv(file="4_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p5 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 5")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p5 <- p5+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p5
library(eoffice)
topptx(p5,"CLU5.pptx",width = 6,height = 5)

#cluster6
library(ggplot2)
hh <- read.csv(file="5_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p6 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 6")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p6 <- p6+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p6
library(eoffice)
topptx(p6,"CLU6.pptx",width = 6,height = 5)

#cluster7
library(ggplot2)
hh <- read.csv(file="6_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p7 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 7")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p7 <- p7+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p7
library(eoffice)
topptx(p7,"CLU7.pptx",width = 6,height = 5)

#cluster8
library(ggplot2)
hh <- read.csv(file="7_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p8 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 8")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p8 <- p8+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p8
library(eoffice)
topptx(p8,"CLU8.pptx",width = 6,height = 5)

#cluster9
library(ggplot2)
hh <- read.csv(file="8_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p9 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 9")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p9 <- p9+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p9
library(eoffice)
topptx(p9,"CLU9.pptx",width = 6,height = 5)
#cluster10
library(ggplot2)
hh <- read.csv(file="9_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p10 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 10")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p10 <- p10+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p10
library(eoffice)
topptx(p10,"CLU10.pptx",width = 6,height = 5)

#cluster11
library(ggplot2)
hh <- read.csv(file="10_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p11 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 11")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p11 <- p11+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p11
library(eoffice)
topptx(p11,"CLU11.pptx",width = 6,height = 5)


#cluster12
library(ggplot2)
hh <- read.csv(file="11_JIEGUO.csv")
library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:10,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="cluster 12")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",size=15,face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12

library(eoffice)
topptx(p12,"CLU12.pptx",width = 6,height = 5)

###G:\P2\P2.7-更新点图形式的KEGG####
setwd("G:\\P2\\P2.7-更新点图形式的KEGG")
a <- read.csv("每个的前五个_调整顺序.csv")
c <- unique(a$ID)

a$clu <-  factor(a$clu,levels = c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12"))
a$order=factor(rev(as.integer(rownames(a))),labels = rev(a$Description))

p12 <- ggplot(a,aes(y=order,x=clu))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12

library(eoffice)
topptx(p12,"KEGG的通路图-.pptx",width = 7,height = 6)


###p2.8zhuanluyinzi####
#8.20
library(Seurat)
library(dorothea)
library(tidyverse)
library(stringi)
load("G:\\yan\\23.8.21\\FUWUQI\\sce.RData")

DefaultAssay(object = sce) <- "dorothea"
table(Idents(sce))


sce <- ScaleData(sce)

sce@assays$dorothea@scale.data[1:4,1:4]

gc()
gene_cell_exp_scaledata <- AverageExpression(sce,
                                             group.by = 'seurat_clusters',
                                             slot = 'scale.data') 


a <- gene_cell_exp_scaledata$dorothea

top5 <- c("E2F4","E2F1","TFDP1","ZFX","NR2C2",
          "MYCN","E2F2","CUX1","MYC","AR",
          "HNF1A","NR5A1","THAP11","MXI1","CDX2",
          "THAP11","HNF1A","NR5A1","CDX2","MXI1",
          "LEF1","TCF7L2","FOXM1","SOX2","HSF1",
          "FOXM1","NFYB","MYBL2","ZNF384","NR2C2",
          "IRF9","STAT2","IRF1","STAT1","IRF2",
          "TCF7L2","CEBPA","JUN","PPARA","CEBPD",
          "FOSL1","ZEB2","FOXP1","LEF1","FOS",
          "TP53","TP63","TWIST1","KLF4","SOX9",
          "PRDM14","PDX1","PAX6","THAP11","SPI1",
          "ZNF263","ZEB1","ZBTB7A","SNAI2","MNT")
TOP5 <- unique(top5)


dim(a)

library(pheatmap)
az <- top5
az <- unique(az)
az <- data.frame(az)
rownames(az) <- az[,1]
az1 <- a[rownames(az),]
setwd("G:\\P2\\p2.8zhuanluyinzi")
write.csv(az1,file="huitu.csv")
az1 <- read.csv("huitu.csv",row.names = 1)
getwd()
pheatmap(az1,file="xianp1.tiff",scale = "row",show_rownames=1,color = mycol2, show_colnames=1,cluster_rows= F,cluster_cols = FALSE,  border=FALSE,cellheight=8,cellwidth = 12,fontsize=5,treeheight_row=0)#??????ͼ
#复刻了，下面是更新的


library(ComplexHeatmap)
library(Seurat)
library(dplyr)
library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)
#自定义配色：
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

#归一化：
colnames(az1) <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12")

aver_dtt <- t(scale(t(az1)))
#添加行列注释：

cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
names(cols) <- colnames(az1)

#列注释：
cell <- data.frame(colnames(aver_dtt))
colnames(cell) <- 'cell'

col_anno <- HeatmapAnnotation(df = cell,
                              show_annotation_name = F,
                              gp = gpar(col = 'white', lwd = 2),
                              col = list(cell = cols))

#行注释：
row_cols <- setNames(rep(cols, each = 5), rownames(aver_dtt))
row_cols

row_anno <- rowAnnotation(foo = anno_text(rownames(aver_dtt),
                                          location = 0,
                                          just = "left",
                                          gp = gpar(fill = row_cols,
                                                    col = "black",
                                                    fontface = 'italic'),
                                          width = max_text_width(rownames(aver_dtt))*1.25))

#热图绘制：
p1 <- Heatmap(aver_dtt,
        name = 'expression',
        col = mycol2,
        cluster_columns = F,
        cluster_rows = F,
        column_names_side = c('top'),
        row_names_gp = gpar(fontsize = 12, fontface = 'italic'),
        rect_gp = gpar(col = "white", lwd = 1.5),
        top_annotation = col_anno) + row_anno 
p1
library(eoffice)
topptx(p1,"P1.pptx",width = 6.8,height = 8.63)
#要加个bar值，用常规的pheatmap重新画，
mycol <- colorRampPalette(c("#0da9ce", "white", "#e74a32"))(50)
p2 <- pheatmap(az1,scale = "row",show_rownames=1,color = mycol, show_colnames=1,cluster_rows= F,cluster_cols = FALSE,  border=FALSE,cellheight=8,cellwidth = 12,fontsize=5,treeheight_row=0)#??????ͼ
p2
topptx(p2,"P2.pptx",width = 6.8,height = 8.63)
#把P2的bar口过去了，直接用p1PPT就好

### p3.1inferrcnv heatmap####



setwd("G:\\yan\\22.12.5\\23.11.22huatu_cnv\\zuizhong")
#G:\P3\p3.1infercnvretu已经弄好了，换了样本泼墨的颜色和字体

###G:\P3\p3.2样本之间的小提琴图和中位数####
setwd("G:\\P3\\p3.2样本之间的小提琴图和中位数")
#G:\yan\23.9.4
a <- read.csv("cell_shuxing+cnv_yangben.csv")




library(ggplot2)
library(ggpubr)
a$Sample <-  factor(a$Sample,levels = c("R1-50","R1-30","R1-10","Ctrl","R25-10","R25-30","R25-50"))

cols <- c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
names(cols) <-c("Ctrl","R1-10","R1-30","R1-50","R25-10","R25-30","R25-50")

p4=ggplot(a,aes(x=Sample,y=Log2,fill=Sample)) +
  geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓???)
  #"trim"如果为TRUE(默认???),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部???
  #geom_boxplot(width=0.3,position=position_dodge(0.9))+ #绘制箱线???
  #  scale_fill_manual(values= c("#FF0000","#0000CD"))+ #设置填充的颜???
  labs(title="", x="", y = "log2(CNV_level)") +
  
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",), # 将图表标题居???
        axis.text.x=element_text(size=15), #设置x轴刻度标签的字体显示倾斜角度???45度，并向下调???1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(),#设置x轴标题的字体属???
        axis.title.y=element_text(), #设置y轴标题的字体属???   
        
        panel.grid.major = element_blank(), #不显示网格线
       # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("Ctrl",c("R25-30","R25-50")))
        
     # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("R1-10","R1-30"),c("R1-30","R1-50"),c("R25-10","R25-30"),c("R25-30","R25-50"),c("R1-50","R25-50"))
                                                               
        ) #不显示网格线
#绘图
library(Seurat)
p5 <- p4+stat_summary(fun.data = "mean_se", geom = "point", show.legend = F)+NoLegend()
p6 <- p5+ geom_signif(comparisons = list(c("R1-10","Ctrl"),c("R1-30","R1-50"),c("R25-10","Ctrl"),c("R25-30","R25-50"),c("R1-50","R25-50")),
                test = "t.test",
                step_increase = 0.1,
                tip_length = 0,
                textsize = 8,
                map_signif_level = TRUE)
p6
library(eoffice)
topptx(p6,"需要下移星号换色和加平均值.pptx",width = 6,height = 6.8)
#上了泼墨。换了轮廓

### G:\P3\p3.3cluster之间的比较####
setwd("G:\\P3\\p3.3cluster之间的比较和中位数柱状图")

c <- read.csv("所有IR细胞加上一个IR细胞的整体.csv")

c$IR_cluster <-  factor(c$IR_cluster,levels =c("IR_cluster 9","IR_cluster 3","IR_cluster 5","IR_cluster 6","IR_cluster 4","all_IR cell","IR_cluster 2","IR_cluster 7","IR_cluster 1","IR_cluster 10","IR_cluster 8","IR_cluster 11","IR_cluster 12"))


p4=ggplot(c,aes(x=IR_cluster,y=Log2,fill=IR_cluster)) +
  #geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓???)
  #"trim"如果为TRUE(默认???),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部???
  geom_boxplot(width=0.3,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线???
   scale_fill_manual(values= c("#eddd86","#f8cb7f","#7cd6cf","#9192ab","#f89588","#FF6347","#76da91","#7898E1","#63b2ee","#9987ce","#efa666","#87CEEB","#76da91"))+ #设置填充的颜???
  labs(title="", x="", y = "log2(CNV_level)") +
  
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",), # 将图表标题居???
        axis.text.x=element_text(colour="black",angle=90,vjust = 0.2), #设置x轴刻度标签的字体显示倾斜角度???45度，并向下调???1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(),#设置x轴标题的字体属???
        axis.title.y=element_text(), #设置y轴标题的字体属???   
        
        panel.grid.major = element_blank(), #不显示网格线
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("Ctrl",c("R25-30","R25-50")))
        
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("R1-10","R1-30"),c("R1-30","R1-50"),c("R25-10","R25-30"),c("R25-30","R25-50"),c("R1-50","R25-50"))
        
  )
p4
#cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
 #             1          2        3          4        5         6         7          8        9         10         11       12

p5 <- p4+NoLegend()
p5

df <- c
library(dplyr)

df %>%
  group_by(IR_cluster) %>%
  summarise(median_V=median(Log2)) -> df



df

p5 <-  p5+geom_line(data=df,aes(x=IR_cluster,y=median_V,group=1),color="orange",size=0.2)+geom_point(data=df,aes(x=IR_cluster,y=median_V))
p5



p6 <- p5+ geom_signif(comparisons = list(c("all_IR cell","IR_cluster 4"),c("all_IR cell",c("IR_cluster 2"))),
                      test = "t.test",
                      step_increase = 0.1,
                      tip_length = 0,
                      textsize = 8,
                      map_signif_level = TRUE)
p6
#p5和P6都保存一次，把P6的显著性偷出来



library(eoffice)
topptx(p6,"tup6.pptx",width =5,height = 6.63)

### G:\P3\p3.4umap高低加上圈####
setwd("G:\\P3\\p3.4umap高低加上圈")


load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz) <- "RNA"
Rz


library(ggplot2)
##1新特性
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
##2商务
cols <- c("194f97","#555555","#bd6b08","#00686b","#c82d31","#625ba1","#898989","#9c9800","#007f54","#a195c5","#103667","#f19272")

##3雅致
cols <- c("#3682be","#45a776","#f05326","#eed777","#334f65","#b3974e","#38cb7d","#ddae33","#844bb3","#93c555","#5f6694","#df3881")
##4复古
cols <- c("#0780cf","#765005","#fa6d1d","#0e2c82","#b6b51f","#da1f18","#701866","#f47a75","#009db2","#024b51","#0780cf","#765005")
###xuanze1
names(cols) <- c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10","cluster 11","cluster 12")
DimPlot(Rz2,reduction = "umap",ncol=2,pt.size = 0.1,label = T,cols = cols)
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
P2 <- p2+NoLegend()
#框的方法不合理，自己画
library(eoffice)
topptx(P2,"画圈.pptx",width =5.5,height = 4.5)

#自己加了个圈可以了


### G:\P3\p3.5定义了高低后的UMAP可以和34一起####
setwd("G:\\P3\\p3.5定义了高低后的UMAP可以和34一起")
Rz2

cols=c('#F5A889','#ACD6EC')
p1 <- DimPlot(Rz2,reduction = "umap",ncol=2,pt.size = 0.1,label = T,cols = cols)
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
P2 <- p2+NoLegend()
P2
#框的方法不合理，自己画
library(eoffice)
topptx(P2,"gaodi.pptx",width =5.5,height = 4.5)
### G:\P3\p3.6DOTPLOT####


Rz2

marker1 <- c("CD44",
             "CDCP1",
             "KLF4",
             "MET",
             "TP53",
             "EGFR",
             "KRAS",
             "ERBB2",
             "BRAF",
             "DDR2",
             "STK11",
             "ERO1A",
             "MUC1",
             "TUBA1A",
             "WSB1",
             "CCND1",
             "VEGFA",
             "BCAM",
             "NDRG1",
             "GBP1",
             "MARCKS",
             "SPHK1",
             "CDK5RAP2",
             "PGAM1",
             "ARF5",
             "NBEAL1",
             "PFN2",
             "AP2M1",
             "CHCHD10",
             "TGFB1",
             "LAMB3",
             "SLC7A5",
             "PMEPA1",
             "F3",
             "FLNA",
             "B4GALT1",
             "MKI67",
             "RNASET2")


p <- DotPlot(Rz2, features = marker1) +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle=90,vjust=0.3,hjust=1,colour="black"))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#ACD6EC','#F5A889')) #颜色
p
p3 <- p+theme(
  axis.text.y=element_text(colour="black"),)+
  theme(axis.title.x = element_text(colour="black"),)+
  theme(axis.title.y = element_text(colour="black"),)+
  theme(axis.text.x = element_text(angle=90,vjust=0.3,hjust=1,colour="black"),)+
  theme(legend.text = element_text(colour="black"),
        legend.title = element_text(colour="black"))
p3

library(eoffice)
getwd()
topptx(p3,"dotplot.pptx",width =9.5,height =3.32)

### G:\P3\P3.7高低之间的火山图####
#G:\yan\23.8.7\HUITUQIAN1 但是我改了
setwd("G:\\P3\\P3.7高低之间的火山图")

#火山图不合适，因为几个下调的基因同时▲PCT也是正的，所以重新画柱状那个


a <- read.csv("hutuqian1.csv",row.names = 1)

library(Seurat)
p <- jjVolcano(diffData = a,
               tile.col = cols,
               size  = 4, 
               fontface = 'italic', family = "Times New Roman",
               flip = T,show.legend = F)+NoLegend()
p


library(eoffice)
topptx(p,"HUOSHAN.pptx",width =8.85,height = 3.32)
#在PPT里改了改

### G:\P3\P3.8GSEA####
#G:\yan\22.12.19
setwd("G:\\P3\\P3.8GSEA")
a <- read.csv("GSEA输入.csv")



library(stringr)
library(clusterProfiler)
df = read.csv(file="GSEA输入.csv")
head(df)#查看前面几行
dim(df)#数据总共几行几列
df <- df[,c(1:2)]
geneList<-df$log2FC
names(geneList)=df$gene
geneList=sort(geneList,decreasing = T)  
kegmt<-read.gmt("h.all.v2022.1.Hs.symbols.gmt") 
a <- GSEA(geneList, TERM2GENE=kegmt,
          pvalueCutoff = 1)
a@result 
gsea_results_df <- a@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'GESA_resul_HALLMARK.csv',row.names = FALSE)     

library(ggplot2)
dotplot(a) #出点图 
dotplot(a,color="pvalue")
library(enrichplot)
p1 <- gseaplot2(a, subplots = 1, c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","KEGG_OOCYTE_MEIOSIS","KEGG_PROTEASOME","KEGG_RIBOSOME"),pvalue_table =F,base_size = 11)
library(Seurat)
p1+NoLegend()
library(eoffice)
getwd()
topptx(p1,"不要下面的_P.pptx",width =12,height =6.8)
#另外弄一个不带P值的
getwd()
topptx(p1,"不要下面的.pptx",width =5.5,height =5.32)
###3.8-让我改为两个数据库的柱状图####
setwd("G:\\P3\\P3.8-改为两个柱状图/")
a <- read.csv("HALLMARK画图.csv")
library(ggplot2)
library(tidyverse)
a$order=factor(rev(as.integer(rownames(a))),labels = rev(a$Description))

library(Seurat)
p <- ggplot(data = a, aes(x = NES, y = order)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "", title = "Enriched pathways of three changing gene clusters") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() +
  NoLegend()

p
topptx(p,"hallmark.pptx",width =9,height =6)
#kegg
a <- read.csv("kegg画图.csv")
library(ggplot2)
library(tidyverse)
a$order=factor(rev(as.integer(rownames(a))),labels = rev(a$Description))

library(Seurat)
p <- ggplot(data = a, aes(x = NES, y = order)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "", title = "Enriched pathways of three changing gene clusters") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() +
  NoLegend()

p
topptx(p,"kegg.pptx",width =9,height =6)

### 3.8-让我改为1个数据库的柱状图 ####
setwd("G:\\P3\\P3.8-改为两个柱状图/合为一个/")
a <- read.csv("GSEA输入.csv")



library(stringr)
library(clusterProfiler)
df = read.csv(file="GSEA输入.csv")
head(df)#查看前面几行
dim(df)#数据总共几行几列
df <- df[,c(1:2)]
geneList<-df$log2FC
names(geneList)=df$gene
geneList=sort(geneList,decreasing = T)  
kegmt<-read.gmt("c2.cp.keggheHALLMERK.gmt") 
a <- GSEA(geneList, TERM2GENE=kegmt,
          pvalueCutoff = 1)
a@result 
gsea_results_df <- a@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'GESA_resul_HALLMARKHEKEGG.csv',row.names = FALSE)     
a <- read.csv("合为一个的画图.csv")
library(ggplot2)
library(tidyverse)
a$order=factor(rev(as.integer(rownames(a))),labels = rev(a$Description))

library(Seurat)
p <- ggplot(data = a, aes(x = NES, y = order)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "", title = "Enriched pathways of three changing gene clusters") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() +
  NoLegend()

p
getwd()
topptx(p,"heweiyige.pptx",width =8,height =4)



###P3.7.5上调的348的KEGG####
# G:\yan\23.7.31\higcnv vs lowcnv

setwd("G:\\P3\\P3.7.5上调的348的KEGG")
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("348updeg_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:20,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*p.adjust))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"TOP20_UP348GENE.pptx",width =6.45,height =5.69)

### G:\P3\P3.9cancermarker和UPDEGs交集####
#ZUOWANLE PPT

### G:\P3\P3.10流式分选图####
找聂老师

### G:\P3\p3.13highcnv降维图####
setwd("G:\\P3\\p3.13highcnv降维图")
#0.2 resolution
#不能要加对照的，不然说不清楚
library(Seurat)
load("G:\\yan\\23.2.27\\highCNV_Rz2.RData")
high_Rz
table(high_Rz@active.ident)

newname <- c("subcluster 1","subcluster 2","subcluster 3","subcluster 4","subcluster 5","subcluster 6","subcluster 7")
names(newname) <- levels(high_Rz)
high_Rz1 <- RenameIdents(high_Rz,newname)
DimPlot(high_Rz1,reduction = "umap",ncol=2,pt.size = 0.1,label = T)
high_Rz <- high_Rz1
p1 <- DimPlot(high_Rz,reduction="tsne",label = T)

library(ggplot2)
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
P2 <- p2+NoLegend()
P2

library(eoffice)
getwd()
topptx(P2,"p3.13highcnv降维图.pptx",width =5.6,height =5.5)


### G:\P3\p3.14HIGHCNV比例柱状图####
setwd("G:\\P3\\p3.14HIGHCNV比例柱状图")
#2.27wenjian li de huitubili.csv
a <- read.csv("huitu_bili.csv")
library(ggplot2)
##huitu
a <- read.csv("huitu_bili.csv")


library(ggplot2)
library(ggplot2)
allcolour=c("#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")

windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))

a <- read.csv("huitu_bili.csv")

a$Sample <-  factor(a$Sample,levels = c("R1-10","R1-30","R1-50","R25-10","R25-30","R25-50"))



p1 <- ggplot(a) + 
  geom_bar(aes(x =Cluster, y= num, fill =Sample),stat = "identity",position = "fill",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black",size=0.5, linetype="solid"))+
  #scale_y_discrete(limits = c("Ctrl","R1-10","R1-30","R1-50","R25-10","R25-30","R25-50"))+
  RotatedAxis()

p1

library(eoffice)
getwd()
topptx(p1,"p3.14HIGHCNV比例柱状图.pptx",width =5.6,height =6.5)
#以PPT为准

### P3.15highcnv注释图####



library(Seurat)
load("G:\\yan\\23.2.27\\highCNV_Rz2.RData")
high_Rz
table(high_Rz@active.ident,high_Rz$sample)
allcolour=c("#FF6767","#4DB3E6","#336699","#993333","#E39FF6","#E39FF6","#E39FF6")
newname <- c("H_R25-10&H_R25-30","H_R1-10&H_R1-30","H_R1-50","H_R25-50","Mixture","Mixture","Mixture")
names(newname) <- levels(high_Rz)
high_Rz1 <- RenameIdents(high_Rz,newname)
DimPlot(high_Rz1,reduction = "umap",ncol=2,pt.size = 0.1,label = T,cols = allcolour)
high_Rz <- high_Rz1
p1 <- DimPlot(high_Rz,reduction="tsne",label = T,cols = allcolour)


windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
P2 <- p2+NoLegend()
P2

library(eoffice)
getwd()
topptx(P2,"G:\\P3\\P3.15highcnv注释图\\tu.pptx",width =5.6,height =5.5)

###p3.16两个50带HIGHvsCTRL ####
#G:\yan\23.2.27\50 VS 2ctrl 的两个
#这两个文件是上下调都有，我们最后只保留了上调的用于画图

###G:\P3\p3.17KEGG（可以只做上调） ####
setwd("G:\\P3\\p3.17KEGG（可以只做上调）")

library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)
##kegg
#r150
diff <- read.csv("上调R150.csv")

gene<-bitr(t(diff),fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
gene <- gene$ENTREZID
#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 1, pvalueCutoff=1)

hh <- as.data.frame(kk)#自己记得保存结果哈！
write.csv(hh,file="上调R150_keggres.csv",row.names = F)
#r2550
diff <- read.csv("上调R2550.csv")

gene<-bitr(t(diff),fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
gene <- gene$ENTREZID
#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 1, pvalueCutoff=1)

hh <- as.data.frame(kk)#自己记得保存结果哈！
write.csv(hh,file="上调R2550_keggres.csv",row.names = F)
### G:\P4\P4.1单次照射的图####
setwd("G:\\P4\\P4.1单次照射的图")
#G:\yan\23.4.24\fenkaide
library(monocle)
load("G:\\P4\\P4.1单次照射的图\\分开的\\R1\\r1order5liangci.RData")

library(monocle)
p1 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")+scale_color_viridis_c()


p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

p2


library(eoffice)
getwd()
topptx(p1,"G:\\P4\\P4.1单次照射的图\\R1nishixutu.pptx",width =5,height =5.5)

#加每个cluster细胞的饼图
az <- table(HSMM$State,HSMM$sample)
write.csv(az,file="G:\\P4\\P4.1单次照射的图\\R1-state和sample的比例.csv")

plot_cell_trajectory(HSMM, color_by = "State")

#补充图
#1分样本
setwd("G:\\P4\\P4.1单次照射的图\\R1nishixu\\补充图")
tiqu <- HSMM$sample
allcolour=c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
p1 <- plot_cell_trajectory(HSMM,color_by = "sample")+facet_wrap(~sample,nrow=1)+ scale_color_manual(breaks = c("R1-10", "R1-30", "R1-50"), values=c("#99CCFF","#0099CC","#336699")) 
p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))

p2
topptx(p2,"G:\\P4\\P4.1单次照射的图\\R1nishixu\\补充图\\分样本.pptx",width =12,height =5.5)
#2波浪


library(ggpubr)
df <- pData(HSMM) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容

p1 <- ggplot(df, aes(Pseudotime, colour = sample, fill=sample)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
cols <- c("#99CCFF","#0099CC","#336699")
p2 <- p1+   scale_color_discrete(type = cols) +
   scale_fill_discrete(type = cols)
p2

p2 <- p2 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))

p2
topptx(p2,"G:\\P4\\P4.1单次照射的图\\R1nishixu\\补充图\\密度图.pptx",width =5.34,height =4.58)
#热图

setwd("G:\\P4\\P4.1单次照射的图\\R1nishixu\\基因热图")
my_pseudotime_de <- read.csv("G:\\P4\\P4.1单次照射的图\\分开的\\R1\\Time_diff_all.csv")
my_pseudotime_de1 <- subset(my_pseudotime_de,qval <= 0.01)
load("G:\\P4\\P4.1单次照射的图\\分开的\\R1\\r1order5liangci.RData")
rownames(my_pseudotime_de1) <- my_pseudotime_de1$gene_short_name
table(my_pseudotime_de1$use_for_ordering)

library(dplyr)
p=plot_pseudotime_heatmap(HSMM[row.names(my_pseudotime_de1),], num_clusters=6, show_rownames=F,cores = 4, return_heatmap=T)
topptx(p2,"G:\\P4\\P4.1单次照射的图\\R1nishixu\\基因热图.pptx",width =5.09,height =6.48)
#用tiff图插入到PPT中自己画了个

#KEGG 三个cluster的
#用微信收藏把三个画到一起
z <- read.csv("cluster123整合到一起画图.csv")
library(ggplot2)
# 先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)
z$Description <- factor(z$Description, levels = z$Description)
library(Seurat)
p <- ggplot(data = z, aes(x = -log10(pvalue), y = rev(Description), fill = gene.cluster)) +
  scale_fill_manual(values =c( '#FF83FF','#FF9289','#82B7FF')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "", title = "Enriched pathways of three changing gene clusters") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()

p


library(dplyr)
library(eoffice)
topptx(p,"G:\\P4\\P4.1单次照射的图\\R1nishixu\\基因热图\\基因类群里的富集.pptx",width =5.09,height =6.48)

#单次三个genecluster真正的富集分析，分别的

setwd("G:\\P4\\P4.1单次照射的图\\R1nishixu\\基因热图\\补充文件123的KEGG")
#要每个的前30
#clu2
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("CLU=2_1500_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:30,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"TOP30_clu2GENE.pptx",width =6.45,height =5.69)
#clu1
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("CLU=1_4721_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:30,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"TOP30_clu1GENE.pptx",width =6.45,height =5.69)
#clu3
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("CLU=3_1959_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:30,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"TOP30_clu3GENE.pptx",width =6.45,height =5.69)

###多次的一整套####



setwd("G:\\P4\\P4.1单次照射的图")
#G:\yan\23.4.24\fenkaide
library(monocle)
load("G:\\P4\\P4.2多次照射的图\\R25/R25liangciorder_root3.RData")

library(monocle)
p1 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")+scale_color_viridis_c()


p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

p2


library(eoffice)
getwd()
topptx(p1,"G:\\P4\\P4.2多次照射的图\\R25nishixu\\R25nishixutu.pptx",width =5,height =5.5)

#加每个cluster细胞的饼图
az <- table(HSMM$State,HSMM$sample)
write.csv(az,file="G:\\P4\\P4.2多次照射的图\\R25nishixu\\R25-state和sample的比例.csv")

plot_cell_trajectory(HSMM, color_by = "State")

#补充图
#1分样本
setwd("G:\\P4\\P4.2多次照射的图\\R25nishixu\\补充图")
tiqu <- HSMM$sample
allcolour=c("#D3D3D3","#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
p1 <- plot_cell_trajectory(HSMM,color_by = "sample")+facet_wrap(~sample,nrow=1)+ scale_color_manual(breaks = c("R25-10", "R25-30", "R25-50"), values=c("#FFCCCC","#FF6666","#993333")) 
p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))

p2
topptx(p2,"G:\\P4\\P4.2多次照射的图\\R25nishixu\\补充图\\分样本.pptx",width =12,height =5.5)
#2波浪


library(ggpubr)
df <- pData(HSMM) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容

p1 <- ggplot(df, aes(Pseudotime, colour = sample, fill=sample)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
cols <- c("#FFCCCC","#FF6666","#993333")
p2 <- p1+   scale_color_discrete(type = cols) +
  scale_fill_discrete(type = cols)
p2

p2 <- p2 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))

p2
topptx(p2,"G:\\P4\\P4.2多次照射的图\\R25nishixu\\补充图\\密度图.pptx",width =5.34,height =4.58)
#热图

setwd("G:\\P4\\P4.2多次照射的图\\R25nishixu\\基因热图")
my_pseudotime_de <- read.csv("G:\\P4\\P4.2多次照射的图\\R25\\Time_diff_all.csv")
my_pseudotime_de1 <- subset(my_pseudotime_de,qval <= 0.01)
load("G:\\P4\\P4.2多次照射的图\\R25/R25liangciorder_root3.RData")
rownames(my_pseudotime_de1) <- my_pseudotime_de1$gene_short_name
table(my_pseudotime_de1$use_for_ordering)

library(dplyr)
p=plot_pseudotime_heatmap(HSMM[row.names(my_pseudotime_de1),], num_clusters=10, show_rownames=F,cores = 4, return_heatmap=T)
#topptx(p2,"G:\\P4\\P4.2多次照射的图\\R25nishixu\\基因热图.pptx",width =5.09,height =6.48)
#用tiff图插入到PPT中自己画了个

#KEGG 三个cluster的
#用微信收藏把三个画到一起
z <- read.csv("R25-cluster143整合到一起画图.csv")

# 先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)
z$Description <- factor(z$Description, levels = z$Description)

p <- ggplot(data = z, aes(x = -log10(pvalue), y = rev(Description), fill = gene.cluster)) +
  scale_fill_manual(values =c('#FF9289', '#FF81D7','#BEC100')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "", title = "Enriched pathways of three changing gene clusters") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()

p


library(dplyr)
library(eoffice)
topptx(p,"G:\\P4\\P4.2多次照射的图\\R25nishixu\\基因热图\\基因类群里的富集.pptx",width =5.09,height =6.48)

#单次三个genecluster真正的富集分析，分别的

setwd("G:\\P4\\P4.2多次照射的图\\R25nishixu\\基因热图\\R25补充文件143的KEGG")
#要每个的前30,但1只有18个 所以1是17个 43是30个
#clu1
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("CLU=1_981_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:17,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"R25-TOP17_clu1GENE.pptx",width =6.45,height =5.69)
#clu4
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("CLU=4_3421_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:30,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"R25-TOP30_clu4GENE.pptx",width =6.45,height =5.69)
#clu3
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)

library(ggplot2)
hh <- read.csv("CLU=3_2981_keggres.csv")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))



library(DOSE)
hh$GeneRatio <- parse_ratio(hh$GeneRatio)
hh <- hh[1:30,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p12 <- ggplot(hh,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=1*pvalue))+# 淇敼鐐圭殑澶у皬
  scale_color_gradient(low="red",high = "blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="Gene Ratio",y="",title="")+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p12 <- p12+ theme(
  axis.text.y=element_text(colour="black",face="bold",family = 'Times New Roman'),)+
  theme(axis.text.x = element_text(colour="black",face="bold",family = "Times New Roman"),)+
  theme(legend.text = element_text(colour="black",face = 'bold',family = 'Times New Roman'),
        legend.title = element_text(colour="black",face = 'bold',family = 'Times New Roman'))
p12
library(eoffice)
getwd()
topptx(p12,"R25-TOP30_clu3GENE.pptx",width =6.45,height =5.69)


###G:\P4\P4.3.5--RAD51和其他基因？####
setwd("G:\\P4\\P4.3.5--RAD51和其他基因？")
load("G:\\yan\\24.1.22\\high_Rz.RData")
p1 <- DotPlot(high_Rz, features = c("RAD51","BRCA2","RPA2"), cols = c("white","red"),group.by = "sample", dot.scale = 8, 
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P4\\P4.3.5--RAD51和其他基因？")
topptx(p2,"G:\\P4\\P4.3.5--RAD51和其他基因？\\HR.pptx",width =7.4,height =5.4)
#NHEJ
p1 <- DotPlot(high_Rz, features = c("XRCC6","XRCC5"), cols = c("white","red"),group.by = "sample", dot.scale = 8, 
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P4\\P4.3.5--RAD51和其他基因？")
topptx(p2,"G:\\P4\\P4.3.5--RAD51和其他基因？\\nhej.pptx",width =7.4,height =5.4)


###G:\P4\P4.4.修复通路 #####
az <- read.csv("G:\\P4\\P4.4.修复通路\\HUATU.csv")
library(ggplot2)
library(ggpubr)
colnames(az)
p4=ggplot(az,aes(x=SAMPLE,y=KEGG_HOMOLOGOUS_RECOMBINATION,fill=SAMPLE)) +
  geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓???)
  #"trim"如果为TRUE(默认???),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部???
  geom_boxplot(width=0.15,position=position_dodge(0.9))+ #绘制箱线???
  scale_fill_manual(values= c("#0099CC","#336699","#FF6666","#993333"))+ #设置填充的颜???
  labs(title="KEGG_HOMOLOGOUS_RECOMBINATION", x="", y = " ") +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",), # 将图表标题居???
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black"), #设置x轴刻度标签的字体显示倾斜角度???45度，并向下调???1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(),#设置x轴标题的字体属???
        axis.title.y=element_text(), #设置y轴标题的字体属???
        
        panel.grid.major = element_blank(), #不显示网格线
       # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("R1-30","R1-50"),c("R25-30","R25-50"))
                                                               
        ) #不显示网格线
#绘图
p6 <- p4+ geom_signif(comparisons = list(c("H_R1-30","H_R1-50"),c("H_R25-30","H_R25-50")),
                         test = "t.test",
                         step_increase = 0.1,
                         tip_length = 0,
                         textsize = 8,
                         map_signif_level = TRUE)+NoLegend()
p6

library(eoffice)
getwd()
setwd("G:\\P4\\P4.4.修复通路")
topptx(p6,"HR.pptx",width =4.5,height =5.5)


###G:\P5\P5.1相互作用数量 ####
setwd("G:\\P5\\P5.1相互作用数量")




library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(Seurat)
library(CellChat)



gc()
load("G:/yan/23.5.22/ctrl/CTRL_cellchat.RData")
load("G:/yan/23.5.22/R110/R110_cellchat.RData")
load("G:/yan/23.5.22/R130/R130_cellchat.RData")
load("G:/yan/23.5.22/R150/R150_cellchat.RData")
load("G:/yan/23.5.22/R2510/R2510_cellchat.RData")
load("G:/yan/23.5.22/R2530/R2530_cellchat.RData")
load("G:/yan/23.5.22/R2550/R2550_cellchat.RData")
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")

az <- cellchat_ctrl@idents
#ctrl-shuliang
groupSize <- as.numeric(table(cellchat_ctrl@idents))

P1 <- netVisual_circle(cellchat_ctrl@net$count, vertex.weight = groupSize,color.use=cols,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"ctrl-shuliang.pptx",width =5,height =5)
#R1-10-shuliang
groupSize <- as.numeric(table(cellchat_R110@idents))

P1 <- netVisual_circle(cellchat_R110@net$count, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"R110-shuliang.pptx",width =5,height =5)
#R1-30-shuliang
groupSize <- as.numeric(table(cellchat_R130@idents))

P1 <- netVisual_circle(cellchat_R130@net$count, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"R130-shuliang.pptx",width =5,height =5)

#R1-50-shuliang
groupSize <- as.numeric(table(cellchat_R150@idents))

P1 <- netVisual_circle(cellchat_R150@net$count, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"R150-shuliang.pptx",width =5,height =5)




#R25-10-shuliang
groupSize <- as.numeric(table(cellchat_R2510@idents))

P1 <- netVisual_circle(cellchat_R2510@net$count, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"R2510-shuliang.pptx",width =5,height =5)
#R25-30-shuliang
groupSize <- as.numeric(table(cellchat_R2530@idents))

P1 <- netVisual_circle(cellchat_R2530@net$count, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"R2530-shuliang.pptx",width =5,height =5)

#R25-50-shuliang
groupSize <- as.numeric(table(cellchat_R2550@idents))

P1 <- netVisual_circle(cellchat_R2550@net$count, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Number of interactions")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.1相互作用数量")
topptx(P1,"R2550-shuliang.pptx",width =5,height =5)


###G:\P5\P5.2相互作用强度 ####
setwd("G:\\P5\\P5.2相互作用强度")




library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(Seurat)
library(CellChat)



gc()
load("G:/yan/23.5.22/ctrl/CTRL_cellchat.RData")
load("G:/yan/23.5.22/R110/R110_cellchat.RData")
load("G:/yan/23.5.22/R130/R130_cellchat.RData")
load("G:/yan/23.5.22/R150/R150_cellchat.RData")
load("G:/yan/23.5.22/R2510/R2510_cellchat.RData")
load("G:/yan/23.5.22/R2530/R2530_cellchat.RData")
load("G:/yan/23.5.22/R2550/R2550_cellchat.RData")
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")

az <- cellchat_ctrl@idents
#ctrl-qiangdu
groupSize <- as.numeric(table(cellchat_ctrl@idents))

P1 <- netVisual_circle(cellchat_ctrl@net$weight, vertex.weight = groupSize,color.use=cols,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")


library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"ctrl-qiangdu.pptx",width =5,height =5)
#R1-10-qiangdu
groupSize <- as.numeric(table(cellchat_R110@idents))

P1 <- netVisual_circle(cellchat_R110@net$weight, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Interaction weights/strength")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"R110-qiangdu.pptx",width =5,height =5)
#R1-30-qiangdu
groupSize <- as.numeric(table(cellchat_R130@idents))

P1 <- netVisual_circle(cellchat_R130@net$weight, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Interaction weights/strength")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"R130-qiangdu.pptx",width =5,height =5)

#R1-50-qiangdu
groupSize <- as.numeric(table(cellchat_R150@idents))

P1 <- netVisual_circle(cellchat_R150@net$weight, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Interaction weights/strength")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"R150-qiangdu.pptx",width =5,height =5)




#R25-10-qiangdu
groupSize <- as.numeric(table(cellchat_R2510@idents))

P1 <- netVisual_circle(cellchat_R2510@net$weight, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Interaction weights/strength")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"R2510-qiangdu.pptx",width =5,height =5)
#R25-30-qiangdu
groupSize <- as.numeric(table(cellchat_R2530@idents))

P1 <- netVisual_circle(cellchat_R2530@net$weight, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Interaction weights/strength")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"R2530-qiangdu.pptx",width =5,height =5)

#R25-50-qiangdu
groupSize <- as.numeric(table(cellchat_R2550@idents))

P1 <- netVisual_circle(cellchat_R2550@net$weight, vertex.weight = groupSize,color.use=cols,
                       weight.scale = T, label.edge= F,
                       title.name = "Interaction weights/strength")
library(eoffice)
getwd()
setwd("G:\\P5\\P5.2相互作用强度")
topptx(P1,"R2550-qiangdu.pptx",width =5,height =5)


setwd("G:\\P5\\P5.3数量强度比较")


object.list <- list(R110 = cellchat_R110, R130 = cellchat_R130,R150 = cellchat_R150,R2510 = cellchat_R2510,R2530 = cellchat_R2530,R2550 = cellchat_R2550)

#run netAnalysis_computeCentrality
object.list = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})

# merge data
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cols <- c("#99CCFF","#0099CC","#336699","#FFCCCC","#FF6666","#993333")
gg1 <- compareInteractions(cellchat, color.use=cols,show.legend = F, measure = "count",group = c(1,2,3,4,5,6))

gg2 <- compareInteractions(cellchat, color.use=cols,show.legend = F, group = c(1,2,3,4,5,6), measure = "weight")
gg3 <- gg1 + gg2

library(eoffice)
getwd()
setwd("G:\\P5\\P5.3数量强度比较")
topptx(gg3,"P5.3数量强度比较.pptx",width =9.68,height =6.68)
save(cellchat,file="G:\\P5\\P5.3数量强度比较\\整合了六个的cellchat.RData")
###G:\P5\p5.4r1强度变化热图 ####
setwd("G:\\P5\\p5.4r1强度变化热图")
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use = cols,comparison = c(2,3))
gg2
library(eoffice)
getwd()
setwd("G:\\P5\\p5.4r1强度变化热图")
topptx(gg2,"P5.3数量强度比较1.pptx",width =5.89,height =5.58)
#R130VSR1-10
setwd("G:\\P5\\p5.4r1强度变化热图\\R130VSR110")
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use = cols,comparison = c(1,2))
gg2
library(eoffice)
getwd()
setwd("G:\\P5\\p5.4r1强度变化热图\\R130VSR110")
topptx(gg2,"P5.3数量强度比较1_r130vsr110.pptx",width =5.89,height =5.58)

### G:\P5\P5.5R25数量变化热图####

setwd("G:\\P5\\P5.5R25数量变化热图")
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
gg2 <- netVisual_heatmap(cellchat, measure = "count",color.use = cols,comparison = c(5,6))
gg2
library(eoffice)
getwd()
setwd("G:\\P5\\P5.5R25数量变化热图")
topptx(gg2,"P5.5R25数量变化热图.pptx",width =5.89,height =5.58)
#R2530VSR2510
setwd("G:\\P5\\P5.5R25数量变化热图\\R25-30VSR2510/")
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
gg2 <- netVisual_heatmap(cellchat, measure = "count",color.use = cols,comparison = c(4,5))
gg2
library(eoffice)
getwd()
setwd("G:\\P5\\P5.5R25数量变化热图\\R25-30VSR2510/")
topptx(gg2,"P5.5R25数量变化热图_R2530VSR2510.pptx",width =5.89,height =5.58)

### G:\P5\P5.6单次CLU11####
setwd("G:\\P5\\P5.6单次CLU11")
cols <- c("#1A1A1A","#0099CC","#336699")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, color.use = cols,idents.use = "11",comparison = c(2,3))
gg1
library(eoffice)
getwd()
setwd("G:\\P5\\P5.6单次CLU11")
topptx(gg1,"P5.6单次CLU11.pptx",width =5.59,height =4.05)

### G:\P5\P5.65黏着斑和ECM受体比较####
setwd("G:\\P5\\P5.65黏着斑和ECM受体比较")
#G:\\yan\\24.1.22\\3050suoyoudeSSGSEA\\画图只有R13050.csv

az <- read.csv("画图只有R13050.csv")
library(ggplot2)
library(ggpubr)
colnames(az)
p4=ggplot(az,aes(x=SAMPLE,y=KEGG_FOCAL_ADHESION,fill=SAMPLE)) +
  geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓???)
  #"trim"如果为TRUE(默认???),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部???
  geom_boxplot(width=0.15,position=position_dodge(0.9))+ #绘制箱线???
  scale_fill_manual(values= c("#0099CC","#336699"))+ #设置填充的颜???
  labs(title="KEGG_FOCAL_ADHESION", x="", y = " ") +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",), # 将图表标题居???
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black"), #设置x轴刻度标签的字体显示倾斜角度???45度，并向下调???1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(),#设置x轴标题的字体属???
        axis.title.y=element_text(), #设置y轴标题的字体属???
        
        panel.grid.major = element_blank(), #不显示网格线
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("R1-30","R1-50"),c("R25-30","R25-50"))
        
  ) #不显示网格线
#绘图
p6 <- p4+ geom_signif(comparisons = list(c("R1-30","R1-50")),
                      test = "t.test",
                      step_increase = 0.1,
                      tip_length = 0,
                      textsize = 8,
                      map_signif_level = TRUE)+NoLegend()
p6

library(eoffice)
getwd()
setwd("G:\\P5\\P5.65黏着斑和ECM受体比较")
topptx(p6,"黏着斑_R13050.pptx",width =4.5,height =5.5)


### G:\P5\P5.7多次热图两个####
setwd("G:\\P5\\P5.7多次热图两个")

library(ComplexHeatmap)

pathway.union <- union(object.list[[5]]@netP$pathways, object.list[[6]]@netP$pathways)
object.list
##all outgoing incoming
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
library(RColorBrewer)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[5]], color.use = cols, pattern = "all", signaling = pathway.union, title = names(object.list)[5], width = 5, height = 6,color.heatmap = "Reds")
ht1
ht2 = netAnalysis_signalingRole_heatmap(object.list[[6]],  color.use = cols,pattern = "all", signaling = pathway.union, title = names(object.list)[6], width = 5, height = 6,color.heatmap = "Reds")
pp1 <- ht1+ht2
pp1
library(eoffice)
getwd()
setwd("G:\\P5\\P5.7多次热图两个")
topptx(pp1,"P5.7多次热图两个2.pptx",width =9.89,height =5.58)

### G:\P5\P5.8多次信息流####
setwd("G:\\P5\\P5.8多次信息流")
cols <- c("#FF6666","#993333")
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, color.use = cols, do.stat = TRUE,comparison = c(5,6))
gg2
library(eoffice)
getwd()
setwd("G:\\P5\\P5.8多次信息流")
topptx(gg2,"P5.8多次信息流.pptx",width =6,height =5.5)
#PPT改了字体的颜色



###G:\P5\P5.85两种变化的三种类型饼图 ####
#G:\yan\23.6.5
setwd("G:\\P5\\P5.85两种变化的三种类型饼图")
#显示提取了一个去重的表格，然后根据表格中这三个通路类型的数字绘制饼图 两个
#保存一下点图
#R25

##只能两两比


object.list <- list(R2530 = cellchat_R2530,R2550 = cellchat_R2550)

#run netAnalysis_computeCentrality
object.list = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})
pos.dataset = "R2550"
net.up <- read.csv("shangtiaodejiyin_R2550VS2530.csv")
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_R2550)


pairLR.use.up = net.up[,"interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12),
                        comparison = c(1, 2),
                        max.dataset = 2,
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

gg1

setwd("G:\\P5\\P5.85两种变化的三种类型饼图/")
library(eoffice)
getwd()
topptx(gg1,"R25所有方法三的dotplot.pptx",width =19.19,height =10.16)
#R1

##只能两两比


object.list <- list(R130 = cellchat_R130,R150 = cellchat_R150)

#run netAnalysis_computeCentrality
object.list = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})
pos.dataset = "R150"
net.up <- read.csv("shangtiaodejiyin_R150VS30.csv")
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_R150)


pairLR.use.up = net.up[,"interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12),
                        comparison = c(1, 2),
                        max.dataset = 2,
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

gg1

setwd("G:\\P5\\P5.85两种变化的三种类型饼图/")
library(eoffice)
getwd()
topptx(gg1,"R1所有方法三的dotplot.pptx",width =19.19,height =10.16)

#只画单次的那三个的



object.list <- list(R130 = cellchat_R130,R150 = cellchat_R150)

#run netAnalysis_computeCentrality
object.list = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})

# merge data
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

gg1 <- netVisual_bubble(cellchat,
                        sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12),
                        comparison = c(1, 2),
                        max.dataset = 2,
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

gg1
gg2 <-gg1$data 
write.csv(gg2,file="单次需要整合的那三个的1.CSV")

h1 <- read.csv("单次需要整合的那三个的_WANCHENG.CSV.csv")


hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p4 <- ggplot(h1,aes(y=name,x=x))+
  geom_point(aes(size=log2(prob.original),color=log2(prob.original)))+# 淇敼鐐圭殑澶у皬
  scale_color_gradientn(colours = c('#5A55A6','#83CDA4','#ECEE93','#F47044','#9F0242'))+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p4

setwd("G:\\P5\\P5.85两种变化的三种类型饼图/")
library(eoffice)
getwd()
topptx(p4,"R1那三条的.pptx",width =7.5,height =3.5)
#要个bar值
p4 <- ggplot(h1,aes(y=name,x=x))+
  geom_point(aes(size=prob.original,color=prob.original))+# 淇敼鐐圭殑澶у皬
  scale_color_gradientn(colours = c('#5A55A6','#83CDA4','#ECEE93','#F47044','#9F0242'))+
  theme_bw()+
  theme(text=element_text(colour="black",face="bold",family = 'Times New Roman'))
p4

library(eoffice)
getwd()
topptx(p4,"R1那三条的——bar.pptx",width =7.5,height =3.5)




#把多次的这对加到liana的表格中


object.list <- list(R2530 = cellchat_R2530,R2550 = cellchat_R2550)

#run netAnalysis_computeCentrality
object.list = lapply(object.list, function(x){
  x = netAnalysis_computeCentrality(x)
})

# merge data
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

gg1 <- netVisual_bubble(cellchat,
                        sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12),
                        comparison = c(1, 2),
                        max.dataset = 2,
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

gg1
gg2 <-gg1$data 
write.csv(gg2,file="多需要整合的那1个的1.CSV")
#就不画dotplot了，直接加在liana里面就好了

#用方法二重新统计饼图
#G:\P5\P5.8.5两种变化的三种类型饼图\用方法二重新统\R1所有增加的
#直接统计就画好了，30代的受配体对50代其实都有，直接用50代就好。



### G:\P5\P5.9贡献####
setwd("G:\\P5\\P5.9贡献")

pp1 <- netAnalysis_contribution(cellchat_R2550, signaling = "ANGPTL")
pp1
library(eoffice)
getwd()
setwd("G:\\P5\\P5.9贡献")
topptx(pp1,"P5.9贡献.pptx",width =5.23,height =4.11)
#PPT砍掉了空白 改了颜色 



###G:\P6\P6.1ANGPTL4-SDC4在自己样本的表达 ####
library(Seurat)


load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)

DefaultAssay(Rz) <- "RNA"
Rz

DimPlot(Rz1,reduction = "umap",ncol=2,pt.size = 0.1,label = T)


newname <- c("High_CNV cell","High_CNV cell","Low_CNV cell","Low_CNV cell","Low_CNV cell","Low_CNV cell","High_CNV cell","High_CNV cell","Low_CNV cell","High_CNV cell","High_CNV cell","High_CNV cell")
names(newname) <- levels(Rz1)
Rz1 <- RenameIdents(Rz1,newname)
DimPlot(Rz1,reduction = "umap",ncol=2,pt.size = 0.1,label = T)
Rz1@meta.data$cnv <- Rz1@active.ident

p1 <- DotPlot(Rz1, features ="ANGPTL4", cols = c("white","red"), dot.scale = 8
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.1ANGPTL4-SDC4在自己样本的表达")
topptx(p2,"1-P6.1ANGPTL4-SDC4在自己样本的表达.pptx",width =6.4,height =5.4)


p1 <- DotPlot(Rz1, features ="SDC4", cols = c("white","red"), dot.scale = 8
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.1ANGPTL4-SDC4在自己样本的表达")
topptx(p2,"2-P6.1ANGPTL4-SDC4在自己样本的表达.pptx",width =6.4,height =5.4)
###G:\P6\P6.2ANGPTL4-SDC4在三个多次样本的表达####
setwd("G:\\P6\\P6.2ANGPTL4-SDC4在三个多次样本的表达")
Rz2 <- Rz1[,Rz1@active.ident %in% c("High_CNV cell") ]

Rz3 <- Rz2[,Rz2@meta.data$sample %in% c("R25-10","R25-30","R25-50") ]
Rz3
table(Rz3$sample)


p1 <- DotPlot(Rz3, features ="ANGPTL4", cols = c("white","red"), dot.scale = 8,group.by = "sample"
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()

setwd("G:\\P6\\P6.2ANGPTL4-SDC4在三个多次样本的表达")
topptx(p2,"1-P6.2ANGPTL4-SDC4在三个多次样本的表达.pptx",width =6.4,height =5.4)


p1 <- DotPlot(Rz3, features ="SDC4", cols = c("white","red"), dot.scale = 8,group.by = "sample"
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.2ANGPTL4-SDC4在三个多次样本的表达")
topptx(p2,"2-P6.2ANGPTL4-SDC4在三个多次样本的表达.pptx",width =6.4,height =5.4)


### G:\P6\P6.4预后####

setwd("G:\\P6\\P6.4预后")



#setwd("G:\\yan\\24.2.26\\TCGA")

library(survival)
library(survminer)

a <- read.csv("G:\\yan\\24.2.26\\TCGA\\所有的LUAD+LUSC(1+2).csv")
b <- read.csv("G:\\yan\\24.2.26\\TCGA\\LUAD+LUSC生存.csv")
c <- merge(a,b,gene.by="sample")
table(c$sample_type)




table(c$sample_type)

colnames(surdata)

library(survival)
library(survminer)
surdata <- c




surdata <- c
surdata$OS.time <- surdata$OS.time/365
surdata$level <- ifelse(surdata[,"ANGPTL4"]>median(surdata[,"ANGPTL4"]),'High','Low')
table(surdata$level)
fit <- survfit(Surv(OS.time,OS)~level,data = surdata)
surv_pvalue(fit)$pval
surv_pvalue(fit)$hr
p1 <- ggsurvplot(fit,pval = T,pval.method = T,risk.table = T)
p1
p2 <- p1$plot

library(eoffice)
getwd()
setwd("G:\\P6\\P6.4预后")
topptx(p2,"angptg4_P64预后.pptx",width =4.63,height =4.68)


library(survival)
data.survdiff <- survdiff(Surv(OS.time,OS) ~ level,data = surdata)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 

#SDC4
library(survival)
library(survminer)

a <- read.csv("G:\\yan\\24.2.26\\TCGA\\所有的LUAD+LUSC(1+2).csv")
b <- read.csv("G:\\yan\\24.2.26\\TCGA\\LUAD+LUSC生存.csv")
c <- merge(a,b,gene.by="sample")
table(c$sample_type)




table(c$sample_type)

colnames(surdata)

library(survival)
library(survminer)
surdata <- c




surdata <- c
surdata$OS.time <- surdata$OS.time/365
surdata$level <- ifelse(surdata[,"SDC4"]>median(surdata[,"SDC4"]),'High','Low')
table(surdata$level)
fit <- survfit(Surv(OS.time,OS)~level,data = surdata)
surv_pvalue(fit)$pval
surv_pvalue(fit)$hr
p1 <- ggsurvplot(fit,pval = T,pval.method = T,risk.table = T)
p1
p2 <- p1$plot

library(eoffice)
getwd()
setwd("G:\\P6\\P6.4预后")
topptx(p2,"SDC4_P64预后.pptx",width =4.63,height =4.68)

table(surdata$level)
library(survival)
data.survdiff <- survdiff(Surv(OS.time,OS) ~ level,data = surdata)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 



#双基因

setwd("G:\\yan\\24.2.26\\TCGA")
a <- read.csv("3月5日的双基因的_LUAD+LUSC1+@.CSV.csv")
s <- a
table(s$level1)
s$level <- ifelse(s[,"zonghe"]>median(s[,"zonghe"]),'High','Low')
table(s$level)
fit <- survfit(Surv(OS.time,OS)~level,data = s)
surv_pvalue(fit)$pval
p1 <- ggsurvplot(fit,pval = T,pval.method = T,risk.table = T)
p1
p2 <- p1$plot

library(eoffice)
getwd()
setwd("G:\\P6\\P6.4预后")
topptx(p2,"双基因_P64预后.pptx",width =4.63,height =4.68)

table(s$level)
library(survival)
data.survdiff <- survdiff(Surv(OS.time,OS) ~ level,data = s)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 

###G:\P6\P6.5肺腺癌注释####
setwd("G:\\P6\\P6.5肺腺癌注释")
#G:\\yan\\24.1.22\\MAIOKEYAN\\自己注释玩.RData
load("G:\\yan\\24.1.22\\MAIOKEYAN\\自己注释玩.RData")
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))
cols <- c("#63b2ee","#76da91","#f8cb7f","#f89588","#7cd6cf","#9192ab","#7898E1","#efa666","#eddd86","#9987ce","#87CEEB","#76da91")
p1 <- DimPlot(Rz2,reduction = "umap",pt.size = 0.1,label = T,cols = cols)+NoLegend()
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.5肺腺癌注释/")
topptx(p2,"zhushi.pptx",width =5.6,height =5.5)
###G:\P6\P6.6肺腺癌dotplot####
setwd("G:\\P6\\P6.6肺腺癌dotplot")



p1 <- DotPlot(Rz2,  features = c("CD3D","NKG7","CD68","EPCAM","CD79A","KIT","COL1A1","VWF"), cols = c("white","red"), dot.scale = 8,
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 

library(eoffice)
getwd()
setwd("G:\\P6\\P6.6肺腺癌dotplot/")
topptx(p2,"DOTPLOT.pptx",width =8.5,height =6)


###G:\P6\P6.7ANGPTL4-SDC4在TUMOENORMAL的表达####
setwd("G:\\P6\\P6.7ANGPTL4-SDC4在TUMOENORMAL的表达/")

table(Rz2$group)
p1 <- DotPlot(Rz2, features ="ANGPTL4", cols = c("white","red"), dot.scale = 8,group.by = "group"
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.7ANGPTL4-SDC4在TUMOENORMAL的表达")
topptx(p2,"P6.7ANGPTL4在TUMOENORMAL的表达.pptx",width =6.4,height =5.4)

#sdc4
setwd("G:\\P6\\P6.7ANGPTL4-SDC4在TUMOENORMAL的表达/")

table(Rz2$group)
p1 <- DotPlot(Rz2, features ="SDC4", cols = c("white","red"), dot.scale = 8,group.by = "group"
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.7ANGPTL4-SDC4在TUMOENORMAL的表达")
topptx(p2,"P6.7SDC4在TUMOENORMAL的表达.pptx",width =6.4,height =5.4)

###G:\P6\P6.8ANGPTL4-SDC4在细胞类型的表达####
setwd("G:\\P6\\P6.8ANGPTL4-SDC4在细胞类型的表达")


table(Rz2@active.ident)
p1 <- DotPlot(Rz2, features ="ANGPTL4", cols = c("white","red"), dot.scale = 8,
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.8ANGPTL4-SDC4在细胞类型的表达")
topptx(p2,"P6.8ANGPTL4在细胞类型的表达.pptx",width =6.4,height =5.4)

#sdc4
setwd("G:\\P6\\P6.8ANGPTL4-SDC4在细胞类型的表达")

table(Rz2$group)
p1 <- DotPlot(Rz2, features ="SDC4", cols = c("white","red"), dot.scale = 8,
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.8ANGPTL4-SDC4在细胞类型的表达")
topptx(p2,"P6.8SDC4在细胞类型的表达达.pptx",width =6.4,height =5.4)


###G:\P6\P6.9互作的验证####
load("G:\\yan\\24.2.26\\缪可言的marker_肿瘤上皮\\cellchat_type=truncatedMean_trim_005.RData")
setwd("G:\\P6\\P6.9互作的验证")


pairLR.use <- extractEnrichedLR(cellchat, signaling = c("ANGPTL"))

pairLR.use  = pairLR.use[c(8),,drop=F]
P1 <- netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(1:2), 
                       pairLR.use =  pairLR.use, remove.isolate = FALSE)
P1
library(eoffice)
getwd()
setwd("G:\\P6\\P6.9互作的验证")
topptx(P1,"P6.9互作的验证.pptx",width =5,height =4)



###G:\P6\P6.10angptl4在RNA-SEQ的表达####
setwd("G:\\P6\\P6.10angptl4在RNA-SEQ的表达")


cx <- read.csv("G:\\yan\\24.2.26\\肺纤维化\\ANGPTL4_所有IPF和对照的fpkm.csv")

a <- cx
library(ggplot2)
library(ggpubr)
table(a$TYPE)
p4=ggplot(a,aes(x=TYPE,y=log2(ANGPTL4+1),fill=TYPE)) +
  
  #geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓???)
  #"trim"如果为TRUE(默认???),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部???
  geom_boxplot(width=0.6,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线???
  # scale_fill_manual(values= c("#FF6666","#993333"))+ #设置填充的颜???
  labs(title="ANGPTL4 expression", x="", y = "") +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",), # 将图表标题居???
        axis.text.x=element_text(size=15), #设置x轴刻度标签的字体显示倾斜角度???45度，并向下调???1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(),#设置x轴标题的字体属???
        axis.title.y=element_text(), #设置y轴标题的字体属???   
        
        panel.grid.major = element_blank(), #不显示网格线
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("Ctrl",c("R25-30","R25-50")))
        
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("R1-10","R1-30"),c("R1-30","R1-50"),c("R25-10","R25-30"),c("R25-30","R25-50"),c("R1-50","R25-50"))
        
  ) #不显示网格线
p4
#绘图
library(Seurat)
p5 <- p4+stat_summary(fun.data = "mean_se", geom = "point", show.legend = F)+NoLegend()
p6 <- p5+ geom_signif(comparisons = list(c("control","ipf")),
                      test = "t.test",
                      step_increase = 0.1,
                      tip_length = 0,
                      textsize = 8,
                      map_signif_level = TRUE)
p6
library(eoffice)
topptx(p6,"G:\\P6\\P6.10angptl4在RNA-SEQ的表达\\P6.10angptl4在RNA-SEQ的表达1.pptx",width = 6,height = 6.8)


###G:\P6\P6.11IPF分数####
setwd("G:\\P6\\P6.11IPF分数")




a <- read.csv("最后画图的low_R253050_res_肺纤维化分数_daiyangben.csv")

library(ggplot2)
library(ggpubr)
table(a$r3050.sample)
p4=ggplot(a,aes(x=r3050.sample,y=FEIXIANWEIHUA,fill=r3050.sample)) +
  
  #geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓???)
  #"trim"如果为TRUE(默认???),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部???
  geom_boxplot(width=0.6,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线???
   scale_fill_manual(values= c("#FF6666","#993333"))+ #设置填充的颜???
  labs(title="IPF signature", x="", y = " ") +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",), # 将图表标题居???
        axis.text.x=element_text(size=15), #设置x轴刻度标签的字体显示倾斜角度???45度，并向下调???1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(),#设置x轴标题的字体属???
        axis.title.y=element_text(), #设置y轴标题的字体属???   
        
        panel.grid.major = element_blank(), #不显示网格线
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("Ctrl",c("R25-30","R25-50")))
        
        # panel.grid.minor = element_blank())+stat_compare_means(method = "t.test",comparisons = list(c("R1-10","R1-30"),c("R1-30","R1-50"),c("R25-10","R25-30"),c("R25-30","R25-50"),c("R1-50","R25-50"))
        
  ) #不显示网格线
#绘图
library(Seurat)
p5 <- p4+stat_summary(fun.data = "mean_se", geom = "point", show.legend = F)+NoLegend()
p6 <- p5+ geom_signif(comparisons = list(c("R25-30","R25-50")),
                      test = "t.test",
                      step_increase = 0.1,
                      tip_length = 0,
                      textsize = 8,
                      map_signif_level = TRUE)
p6
library(eoffice)
topptx(p6,"P6.11IPF分1数.pptx",width = 6,height = 6.8)


###G:\P6\P6.12肺纤维化表型补充FGF和IFN等表达####

setwd("G:\\P6\\P6.12肺纤维化表型补充FGF和IFN等表达")

library(Seurat)
library(CellChat)


load("G:\\yan\\22.11.14\\Rz_0.5.RData")
library(Seurat)


Rz3 <- Rz1[,Rz1@meta.data$sample %in% c("R25-10","R25-30","R25-50") ]
Rz3

p1 <- DotPlot(Rz3, features ="FGFR1", cols = c("white","red"), dot.scale = 8,group.by="sample",
) +scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))



p1
p2 <- p1 + theme(
  axis.text.y=element_text(colour="black",family = 'B'),)+
  theme(axis.title.x = element_text(colour="black",family = "B"),)+
  theme(axis.title.y = element_text(colour="black",family = "B"),)+
  theme(axis.text.x = element_text(colour="black",family = "B"),)+
  theme(legend.text = element_text(colour="black",family = 'B'),
        legend.title = element_text(colour="black",family = 'B'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2 
library(eoffice)
getwd()
setwd("G:\\P6\\P6.12肺纤维化表型补充FGF和IFN等表达")
topptx(p2,"FGFR1.pptx",width =6.4,height =5.4)
