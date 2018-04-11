###TCGA Provisional , PMN-SIGNATURE
source("function_heatmap_cluster_combn_after_jianjun.R")
EDA=function(x){
  windows()
  par(mfrow=c(3,2))
  hist(x)
  dotchart(x)
  boxplot(x,horizontal = T)
  qqnorm(x);qqline(x)
  mtext("title",outer = TRUE)
  par(mfrow=c(1,1))
}
##定义参数列表
para_List=c()
ori_low=67;ori_middle=260;ori_high=171##CXCL17
ori_low=118;ori_middle=206;ori_high=174##cxcl15
ori_low=37;ori_middle=46;ori_high=34##cxcl15
getGene0="CXCL17"
getGene0="CXCL5"
#getGene0="PML"
ref_Gene="ZBTB7A"
ref_Gene1="PTEN"
k_n0=3

##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
##path="D:\\min数据\\tcga\\prad_su2c_2015\\prad_su2c_2015\\"
out_put_dir=paste(path,getGene0,"\\",sep="")
label_Func_readInRnaSeqFile=function(){print("读入RNAseq表达数据！")}
file=paste(path,list.files(path = path,pattern = "data.*?RNA_Seq.*?expression.*?median.*?txt"),sep="")
RDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=file)
useSymbol = FALSE
#file=paste(path,"data_RNA_Seq_v2_mRNA_median_Zscores.RDS",sep="")
if(file.exists(RDSfile)){
  rna_seq_data=readRDS(file=RDSfile)
}else{
  rna_data=read.csv(file = file,sep = "\t")
  saveRDS(rna_data,file=RDSfile)
  rna_seq_data=readRDS(file=RDSfile)
}
##rnaSEQ文件探针合并.max,average??
##这里我们使用data.table的.SD方法
library(data.table)
rna_seq_data=data.table(rna_seq_data)
receive=message(cat("发现重复基因",as.character(nrow(rna_seq_data)!=length(unique(rna_seq_data$Hugo_Symbol)))))
rna_seq_data=rna_seq_data[,lapply(.SD,max),by=.(Hugo_Symbol),.SDcols=colnames(rna_seq_data)[c(2:ncol(rna_seq_data))]]
rna_seq_data=na.omit(as.data.frame(rna_seq_data))
rownames(rna_seq_data)=rna_seq_data$Hugo_Symbol


if(FALSE){##在这里不读取CNA文件了
  label_Func_FileCNA=function(){print("读取CNA文件！")}
  CNAfile=paste(path,list.files(path = path,pattern = "data_CNA.*txt"),sep="")
  CNARDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=CNAfile)
  if(file.exists(CNARDSfile)){
    data_CNA=readRDS(file=CNARDSfile)
  }else{
    data_CNA=read.csv(file = CNAfile,sep = "\t")
    saveRDS(data_CNA,file=CNARDSfile)
    data_CNA=readRDS(file=CNARDSfile)
  }
  fix(data_CNA)
}


label_Func_OutputDir=function(){print("打开输出文件文件夹！")}
if(!dir.exists(paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep=""))){
  dir.create(recursive=TRUE,paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep=""))
}
out_put_dir=paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep="")
shell.exec(path)
shell.exec(paste(out_put_dir,sep=""))


###开始进行聚类
fix(rna_seq_data)
geneSignatures=c("ANPEP","CD38","CSF3","CXCR2","CD80","IDO","SLEB2","CD274","PTPRC","TLR3","CEACAM8","STAT1","CXCR4","LOC116196","CSF1R","ITGAM","ITGAX","CD14","TGFB1","CXCL12","LOC199828","ENTPD1","FUT4","STAT3","IL4R","STAT5A","TLR4","CSF2","CXCL8","S100A8","S100A9","TNF")
# geneSignatures0=c("ARG1","CCL2","CD40","IL4R","FCGR2A","CD68","CD33","CSF1R","ITGAM","HLA-DRB1","FOXP3","PDCDILG2","CD163","LY75","PTPRC","CCR2","CD200R1","IL10","FCER2")
# geneSignatures=c("ARG1","CCL2","CD40","IL4R","FCGR2A","CD68","LOC116196","CSF1R","ITGAM","HLA-DRB1","JM2","PDCD1LG2","CD163","LY75","PTPRC","LOC90262","CD200R1","IL10","FCER2")
data2use=rna_seq_data[geneSignatures,];
rownames(data2use)=geneSignatures;fix(data2use)
###完成了数据的提取，下面进行数据的预处理，分为两部分，复制出一部分进行data_plot
data_plot=data2use
##去除无用的列表
  if("Entrez_Gene_Id" %in% colnames(data_plot)){
    data_plot=data_plot[,-which(colnames(data_plot)=="Entrez_Gene_Id")]  
  }
  if("Hugo_Symbol" %in% colnames(data_plot)){
    data_plot=data_plot[,-which(colnames(data_plot)=="Hugo_Symbol")]
  };fix(data_plot)
##进行数据的处理，例如log2,normalize
library(lattice)
#win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()
{"进行log2处理！"}
data_plot=log2(data_plot+1)
win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()
{"进行scale处理"}
data_plot=data.frame(t(scale(t(data_plot)+1)))
#win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()

####开始进行热图的绘制
##计算break
goHeatmap=function(d1=data_plot,distmethod="euclidean",clustmethod="average"){
  PERCENTILE=0.01;lowQ=as.numeric(quantile(unlist(data_plot),PERCENTILE,na.rm = TRUE));highQ=as.numeric(quantile(unlist(data_plot),1-PERCENTILE,na.rm = TRUE))
  BREAKS=c(min(data_plot)-0.01,seq(lowQ,highQ,0.005),max(data_plot)+0.01)
  library(factoextra)
   print("自定义聚类函数")
  myClust<-function(x,aclustmethod=clustmethod){
    hclust(x,method=aclustmethod)
  }
  print("自定义距离函数")
  myDist<-function(x,adistmethod=distmethod){
    get_dist(x,method = adistmethod)
  }
  
  ##主要的绘图函数
  library(gplots)
  hm=heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1), distfun=myDist,hclustfun=myClust,
               scale="none",
               trace="none",
               dendrogram = "column",
               Rowv = FALSE,
               Colv = TRUE,
               sepcolor = "white",
               symkey = TRUE,
               #breaks = c(min(d1)-0.05,seq(-3,3,0.005),max(d1)+0.05)
               # breaks = c(min(d1)-0.5,seq(-3,3,0.005),max(d1)+0.5)
               breaks = BREAKS,
               main = paste("distmethod",distmethod,"clustermethod",clustmethod,sep="_")
               #,lmat=rbind(c(0,3,0),c(0,1,2),c(0,4,0)) ,lhei=c(3,10,3),lwid=c(1,9,1)
  )
  return(hm)
}



##############开始循环餐数
Cluster_Method<-c( 
  "ward.D",
  "ward.D2",
  #"single",
  "complete",
  "average"# ,
  #"mcquitty",
  #"median",
  #"centroid"
)
Dist_Methods<-  c("euclidean"
                  #, "maximum", "manhattan", 
                  #"canberra", 
                  #"binary", 
                  #"minkowski",
                  #"pearson", "spearman","kendall"
)

for(a in dev.list()){
  dev.off()
}

if(FALSE){
for(onedistmethod in Dist_Methods){
  for(oneclustmethod in Cluster_Method){
    win.graph();
    goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod)
  }
}}



goboxplot=function(d2=data_plot,getGene0_=getGene0,one_dist_method,one_clust_method){
  ###??????????????????
  library(factoextra)
    dist=get_dist(x=t(d2),method = one_dist_method)
    hc=hclust(d=dist,method = one_clust_method)
    hcc=data.frame(cutree(tree = hc,k=7));colnames(hcc)="treeIndex"
    hcc[,"Sample"]=rownames(hcc)
    data_getGene0=data.frame(t(rna_seq_data[getGene0_,])[c(-1,-2),]);colnames(data_getGene0)="Value";data_getGene0[,"Sample"]=rownames(data_getGene0)
    hcc=merge(hcc,data_getGene0,by.x = "Sample",by.y ="Sample",all = TRUE )
    for(one in unique(hcc$treeIndex)){
      ##从热图获取区块数据
      get_block_average=d2[,hcc[which(hcc$treeIndex==one),"Sample"]]
      meanV=mean(as.numeric(as.character(unlist(get_block_average))))
      print(meanV)
      hcc[which(hcc$treeIndex==one),"groupMean"]=meanV
    }
    hcc=hcc[order(hcc$groupMean),]
    ranks=1
    for(onemean in sort(unique(hcc$groupMean))){
      hcc[which(hcc$groupMean==onemean),"rank_low2high"]=paste("low2high",ranks,sep="_")
      ranks=ranks+1
    };cat(sort(unique(hcc$groupMean)),"\n")
    
    rankindex=data.frame(unique(hcc$groupMean));colnames(rankindex)="groupMean";rownames(rankindex)=rankindex$groupMean;
    hcc2=hclust(get_dist(rankindex,method = "euclidean"),method = "average");#win.graph();
    plot(hcc2)
    hcc2=data.frame(cutree(hcc2,k=3));colnames(hcc2)="rank_group_index";hcc2[,"groupMean"]=rownames(hcc2)
    hcc2$rank_group_index=paste("rank_group_index",hcc2$rank_group_index,sep = "_")
     ##合并重新分组后的
    hcc=merge(hcc,hcc2,by.x ="groupMean",by.y="groupMean")
    hcc$Value=as.numeric(hcc$Value)
    for(one in unique(hcc$rank_group_index)){
      hcc[which(hcc$rank_group_index==one),"group_count"]=length(which(hcc$rank_group_index==one))
    }
    
    fix(hcc)
   ##进行排序
    my_comparision_matrix=as.data.frame(t(combn(x=unique(hcc$rank_group_index),2)))
    colnames(my_comparision_matrix)=c("first","second")
    for(onerow in c(1:nrow(my_comparision_matrix))){
      my_comparision_matrix[onerow,"differ"]=as.numeric(stringi::stri_sub(str=my_comparision_matrix[onerow,"second"],from = stringi::stri_length(my_comparision_matrix[onerow,"second"])))-
        as.numeric(stringi::stri_sub(str=my_comparision_matrix[onerow,"first"],from = stringi::stri_length(my_comparision_matrix[onerow,"second"])))
    }
    ##按照差值进行排序
    my_comparision_matrix=my_comparision_matrix[order(my_comparision_matrix$first),]
    my_comparision_matrix=my_comparision_matrix[order(my_comparision_matrix$differ),]
    my_comparisons <- list()
    for(one in c(1:nrow(my_comparision_matrix))){
      oneterm=as.vector(c(as.character(my_comparision_matrix[one,"first"]),as.character(my_comparision_matrix[one,"second"])))
      my_comparisons[[one]]=oneterm
    }
    
    
    ###绘制图形qlot
  library(ggpubr);
  ##十分之一percent，用于绘图
    particle=max(diff(hcc$Value))/10
    #win.graph();
    fix(hcc)
    q=qplot(data = hcc,x=rank_group_index,y=hcc$Value,geom = "boxplot",outlier.colour = "black",outlier.colour="black",
                        ylab = getGene0_,main = paste("dist_method",one_dist_method,"hclust_method",one_clust_method,sep = "_"))
    q=q+ geom_jitter(aes(colour = rank_low2high))
    q=q+geom_text(data = hcc,aes(label=paste("n=",group_count,sep=""),y=min(hcc$Value)-particle))
    q=q+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),
      comparisons = my_comparisons,paired = FALSE,#label = "pb.format",
                           #hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                           label.y = c(seq(max(hcc$Value),max(hcc$Value)+max(diff(hcc$Value))/10*10,particle))[c(1:3)]
                          )
    q=q+stat_compare_means(label.y = c(seq(max(hcc$Value),max(hcc$Value)+max(diff(hcc$Value))/10*10,particle))[4],
                           )
    q=q+geom_hline(yintercept = mean(hcc$Value), linetype=2)
    q=q+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(hcc$Value)-particle/2)# Pairwise comparison against all
    return(q)
}


if(FALSE){
for(onedistmethod in Dist_Methods){
  for(oneclustmethod in Cluster_Method){
    #win.graph();
    q=goboxplot(d2=data_plot
              ,getGene0_ = getGene0
              ,one_dist_method = onedistmethod
              ,one_clust_method = oneclustmethod
              );plot(q)
  }
}
}

#######把图画在一起
plot_together=function(){
  for(onedistmethod in Dist_Methods){
    for(oneclustmethod in Cluster_Method){
      #win.graph();
      h=goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod);
      h
      #win.graph();
      q=goboxplot(d2=data_plot
                  ,getGene0_ = getGene0
                  ,one_dist_method = onedistmethod
                  ,one_clust_method = oneclustmethod
      );plot(q)
      
    }
  }
}






if(FALSE){

#提取zbtb7a
data_refGene=rna_seq_data[ref_Gene,]
if("Hugo_Symbol" %in% colnames(rna_seq_data)){
  data_refGene=data_refGene[,-which(colnames(data_refGene)=="Hugo_Symbol")] 
}
fix(data_refGene)
##进行参考基因的分组
criteria_refGene.high=as.numeric(quantile(data_refGene,1-0.205))
criteria_refGene.low=as.numeric(quantile(data_refGene,1-0.795))

Data_refGene.high=data_refGene[,which(data_refGene>=criteria_refGene.high)]
Data_refGene.low=data_refGene[,which(data_refGene<=criteria_refGene.low)]

Sample_refGene.high=colnames(Data_refGene.high)
Sample_refGene.high=Sample_refGene.high[which(Sample_refGene.high %in% colnames(rna_seq_data))]
Sample_refGene.low=colnames(Data_refGene.low)
Sample_refGene.low=Sample_refGene.low[which(Sample_refGene.low %in% colnames(rna_seq_data))]


##开始获取signature基因
data_signatire=rna_seq_data[getGene0,]

tryCatch({
  if("Entrez_Gene_Id" %in% colnames(rna_seq_data)){
    data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]  
  }
  if("Hugo_Symbol" %in% colnames(rna_seq_data)){
    data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
  }
},finally = {
  "全部去除完毕!"
})

##查看未找到的名称
rownames(data_signatire)=getGene0
fix(data_signatire)
print(summary(is.na(rownames(data_signatire))))
##寻找criteria基因的样品进行分组,绘图时加标签进行分组
lable_Func_preProcess=function(){print("对signature基因进行分组处理！")}
Sample_Alt_Normal=data_CNA[which(data_CNA$Hugo_Symbol==ref_Gene1),]
rownames(Sample_Alt_Normal)=Sample_Alt_Normal$Hugo_Symbol
if( "Hugo_Symbol"%in%colnames(Sample_Alt_Normal)){
  Sample_Alt_Normal=Sample_Alt_Normal[,-which(colnames(Sample_Alt_Normal)=="Hugo_Symbol")]
}
if("Entrez_Gene_Id"%in%colnames(Sample_Alt_Normal)){
  Sample_Alt_Normal=Sample_Alt_Normal[,-which(colnames(Sample_Alt_Normal)=="Entrez_Gene_Id")]
}

Sample_Alt=rownames(as.data.frame(t(Sample_Alt_Normal)[which(Sample_Alt_Normal==-2),]))
Sample_Alt=Sample_Alt[which(Sample_Alt %in% colnames(rna_seq_data))]
Sample_Normal=rownames(as.data.frame(t(Sample_Alt_Normal)[which(abs(Sample_Alt_Normal)<2),]))
Sample_Normal=Sample_Normal[which(Sample_Normal %in% colnames(rna_seq_data))]
##画图
data_plot=as.data.frame(t(data_signatire))
data_plot[Sample_Alt,"ALT"]=paste(ref_Gene1,"ALT",sep = "_")
data_plot[Sample_Normal,"ALT"]=paste(ref_Gene1,"Normal",sep = "_")
data_plot[Sample_refGene.low,paste(ref_Gene,"Level",sep="_")]=paste(ref_Gene,"Low",sep = "_")
data_plot[Sample_refGene.high,paste(ref_Gene,"Level",sep="_")]=paste(ref_Gene,"High",sep = "_")
library(ggpubr);library(ggrepel)
secondMin=min(data_plot[,getGene0][which(data_plot[,getGene0]!=0)])
data_plot[,getGene0]=log2(data_plot[,getGene0]+secondMin)
write.csv(x=data_plot,file = paste(out_put_dir,"tcga_signature.csv",sep=""))


#可视化DEPDC1基因表达谱
ggboxplot(data=data_plot,x="ALT",y="CXCL5", 
          color = "ALT",add = "jitter", legend="none")+ 
  
  geom_hline(yintercept = mean(data_plot[,getGene0]), linetype=2)+# Add horizontal line at base mean 
  stat_compare_means(method = "t.test") # Add global annova p-value 
win.graph();plot(q)
}