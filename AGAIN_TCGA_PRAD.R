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
#getGene0="CXCL5"
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
#geneSignatures=c("ANPEP","CD38","CSF3","CXCR2","CD80","IDO","SLEB2","CD274","PTPRC","TLR3","CEACAM8","STAT1","CXCR4","LOC116196","CSF1R","ITGAM","ITGAX","CD14","TGFB1","CXCL12","LOC199828","ENTPD1","FUT4","STAT3","IL4R","STAT5A","TLR4","CSF2","CXCL8","S100A8","S100A9","TNF")
geneSignatures0=c("ARG1","CCL2","CD40","IL4R","FCGR2A","CD68","CD33","CSF1R","ITGAM","HLA-DRB1","FOXP3","PDCDILG2","CD163","LY75","PTPRC","CCR2","CD200R1","IL10","FCER2")
geneSignatures=c("ARG1","CCL2","CD40","IL4R","FCGR2A","CD68","LOC116196","CSF1R","ITGAM","HLA-DRB1","JM2","PDCD1LG2","CD163","LY75","PTPRC","LOC90262","CD200R1","IL10","FCER2")
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
#win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()
{"进行scale处理"}
data_plot=data.frame(t(scale(t(data_plot),center =TRUE,scale = TRUE)))
#win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()
########对data_plot进行最终的可视化ggpurb
library(ggpubr)
library(reshape)
a=melt(t(data_plot))
colnames(a)=c("Sample","Gene","Value")
win.graph();ggpubr::ggdensity(data=a,x="Value",color = "Gene",add = "mean")
####开始进行热图的绘制
##计算break
goHeatmap=function(d1=data_plot,distmethod="euclidean",clustmethod="average"){
  PERCENTILE=0.01;lowQ=as.numeric(quantile(unlist(data_plot),PERCENTILE,na.rm = TRUE));highQ=as.numeric(quantile(unlist(data_plot),1-PERCENTILE,na.rm = TRUE))
  BREAKS=c(min(data_plot)-0.01,seq(lowQ,highQ,0.05),max(data_plot)+0.01)
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
  library(pheatmap)
  #设置colsidebar的颜色
  # h0<<-heatmap.2_adj(x=as.matrix(d1),Rowv = TRUE, Colv =TRUE)
  # h0.1=as.hclust(h0$colDendrogram)
  # h0.2=data.frame(cutree(h0.1,k=10));colnames(h0.2)="groupindex";h0.2[,"colname"]=rownames(h0.2)
  # a=data.frame(colnames(data_plot));colnames(a)="colname";a$colname=as.character(a$colname)
  # a_left2right=data.frame(a[h0$colInd,]);colnames(a_left2right)="colname"
  # ##merge
  # a_left2right_merged=merge(h0.2,a_left2right,by.x = "colname",by.y = "colname")
  # a_left2right_merged[which(a_left2right_merged$groupindex==1),"color"]="green"
  # a_left2right_merged[which(a_left2right_merged$groupindex==2),"color"]="red"
  # a_left2right_merged[which(a_left2right_merged$groupindex==3),"color"]="blue"
  # lastcolors=colors()[which(regexpr(pattern = ".*red.*|.*green.*|.*blue.*",text=colors())==-1)]
  # for(onecutreeindex in unique(a_left2right_merged$groupindex)){
  #   if(onecutreeindex>=4){
  #   a_left2right_merged[which(a_left2right_merged$groupindex==onecutreeindex),"color"]=lastcolors[onecutreeindex-3]
  #   }
  # }
  # colsidecolors=a_left2right_merged$color
  #     h1=heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1), distfun=myDist,hclustfun=myClust,
  #              scale="none",
  #              trace="none",
  #              dendrogram = "column",
  #              Rowv = TRUE,
  #              Colv =TRUE,
  #              sepcolor = "white",
  #              symkey = TRUE,
  #              #breaks = c(min(d1)-0.05,seq(-3,3,0.005),max(d1)+0.05)
  #              # breaks = c(min(d1)-0.5,seq(-3,3,0.005),max(d1)+0.5)
  #              breaks = BREAKS,
  #              labCol = NA,
  #              key.title = NA,
  #              main = paste("distmethod",distmethod,"clustermethod",clustmethod,sep="_")
  #           #,lmat=rbind(c(0,3,0),c(1,1,2),c(4,4,0)) ,lhei=c(3,10,3),lwid=rbind(c(3,9,1),c(3,9,1),c(3,9,1))
  #           ,lmat=rbind(c(0,0,3,0,0,0),c(0,0,4,4,4,0),c(0,0,1,1,1,0),c(0,0,2,2,2,0),c(0,5,5,5,0,0)) ,lhei=c(1.5,3,0.3,6,3),lwid=rbind(c(1,0.6,8,0.5,0.46,1))
  #           #,lmat=rbind(c(0,3,0),c(0,4,0),c(0,1,0),c(0,2,0),c(0,5,0)) ,lhei=c(5,3,3,3,3),lwid=c(2,9,2)
  #           ,na.rm = TRUE
  #           ,ColSideColors = colsidecolors
  # )
            
      
      
  hm=pheatmap(mat=as.matrix((d1)),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1)
              #,kmeans_k = 3
              ,breaks = BREAKS
              ,cluster_rows = TRUE
              ,cluster_cols = TRUE
               #,cutree_rows = 3
               ,cutree_cols = 3
              ,drows = myDist(d1)
              ,dcols = myDist(t(d1))
              ,clustering_method = clustmethod
              ,main = paste("distmethod",distmethod,"clustermethod",clustmethod,sep="_")
              ,silent = TRUE
              ##添加注释
)
  ##按照hv$tree_col$labels为顺序
  ##添加分组信息
  hccut=as.data.frame(cutree(hm$tree_col,k=3));
  hccut=as.data.frame(hccut[hm$tree_col$labels,]);colnames(hccut)="groupindex";rownames(hccut)=hm$tree_col$labels
  ##添加getgene表达值信息
  hccut[,paste(getGene0,"value",sep="_")]=scale((unlist(rna_seq_data[getGene0,rownames(hccut)])))
  
  
  lastcolors=colors()[which(regexpr(pattern = ".*red.*|.*green.*|.*blue.*",text=colors())==-1)]
  for(onecutreeindex in unique(hccut$groupindex)){
    if(onecutreeindex>=4){
      hccut[which(hccut$groupindex==onecutreeindex),"color"]=lastcolors[onecutreeindex-3]
    }
  }
  colsidecolors=hccut$color
  hm1=pheatmap(mat=as.matrix((d1)),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1)
              #,kmeans_k = 3
              ,breaks = BREAKS
              ,cluster_rows = TRUE
              ,cluster_cols = TRUE
              #,cutree_rows = 3
              ,cutree_cols = 3
              ,drows = myDist(d1)
              ,dcols = myDist(t(d1))
              ,clustering_method = clustmethod
              ,main = paste("distmethod",distmethod,"clustermethod",clustmethod,sep="_")
              ,silent = FALSE
              ##添加注释
              ,annotation_col = hccut
  )
  hm1
}



##############开始循环餐数
Cluster_Method<-c( 
  "ward.D",
  "ward.D2",
  # "single",
   "complete",
  "average" #,
   #"mcquitty",
   #"median",
   #"centroid"
)
Dist_Methods<-  c("euclidean"
                  #, "maximum"
                  , "manhattan" 
                  #,"canberra", 
                  #"binary", 
                  , "minkowski"
                  # "pearson", "spearman","kendall"
)

# for(a in dev.list()){
#   dev.off()
# }

if(TRUE){
for(onedistmethod in Dist_Methods){
  for(oneclustmethod in Cluster_Method){
    windows(width = 60,height=40);
    #tiff(filename = "a.tiff",width = 800,height = 600)
    hv=goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod)
    #dev.off();
    windows(width = 60,height=40);
    hv1=goboxplot(hv_ = hv,getGene0_ = getGene0,d2=data_plot
                  ,one_dist_method = onedistmethod
                  ,one_clust_method = oneclustmethod
                  )
    windows(width = 60,height=40);
    plot(hv1[[1]])
  }
}}



if(TRUE){
  goboxplot=function(hv_,k0=3,getGene0_=getGene0_,d2=data_plot,one_dist_method,one_clust_method){
  ###??????????????????
  library(factoextra)
    
    hcc=data.frame(cutree(tree = hv_$tree_col,k=k0));colnames(hcc)="treeIndex"
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
    #plot(hcc2)
    hcc2=data.frame(cutree(hcc2,k=3));colnames(hcc2)="rank_group_index";hcc2[,"groupMean"]=rownames(hcc2)
    hcc2$rank_group_index=paste("rank_group_index",hcc2$rank_group_index,sep = "_")
     ##合并重新分组后的
    hcc=merge(hcc,hcc2,by.x ="groupMean",by.y="groupMean")
    hcc$Value=log2(as.numeric(hcc$Value))
    
   ###返回list循环输出
    return_list=list()
  #连续型组合,aggregate ("rank_low2high")  使用复制和更改的模式
   #  hcc[,"rank_group_index"]=hcc$rank_low2high
   #  combination_index=c(1:(k0-2))+shift   ##shiftMax=k0-3
   # hcc[which(hcc$rank_group_index %in% paste("low2high_",combination_index,sep = "")),"rank_group_index"]=paste("low2high_",combination_index[3],sep="")
        
    for(one in unique(hcc$rank_group_index)){
      hcc[which(hcc$rank_group_index==one),"group_count"]=length(which(hcc$rank_group_index==one))
    }

    #fix(hcc)
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
      #  fix(hcc)
    q=qplot(data = hcc,x=rank_group_index,y=hcc$Value,geom = "boxplot",outlier.colour = "black",outlier.colour="black",
                        ylab = paste("Log2(",getGene0_,")",sep = ""),main = paste("dist_method",one_dist_method,"hclust_method",one_clust_method,sep = "_"))
    q=q+ geom_jitter(aes(colour = rank_low2high))
    q=q+geom_text(data = hcc,aes(label=paste("n=",group_count,sep=""),y=min(hcc$Value)-particle))
    q=q+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),
      comparisons = my_comparisons,paired = FALSE,#label = "pb.format",
                           #hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                           label.y = c(seq(max(hcc$Value),max(hcc$Value)+max(diff(hcc$Value))/10*10,particle))[c(1:3)]
                          ,method = "wilcox.test")
    q=q+stat_compare_means(label.y = c(seq(max(hcc$Value),max(hcc$Value)+max(diff(hcc$Value))/10*10,particle))[4],
                           )
    q=q+geom_hline(yintercept = mean(hcc$Value), linetype=2)
    q=q+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(hcc$Value)-particle/2)# Pairwise comparison against all
    return_list[[length(return_list)+1]]=q

    return(return_list)
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
}
