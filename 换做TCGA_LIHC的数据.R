print(commandArgs())
###TCGA_LIHC
source("function_heatmap_cluster_combn_after_jianjun.R")
setwd("E:/data/tcga/lihc_tcga/lihc_tcga")
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
fig_order=1
##定义参数列表
para_List=commandArgs();print(para_List[3]);

if(is.na(para_List[8])){message("请输出基因3");quit(save = "no")}
getGene0=para_List[8]
ori_low=67;ori_middle=260;ori_high=171##CXCL17
ori_low=118;ori_middle=206;ori_high=174##cxcl15
ori_low=37;ori_middle=46;ori_high=34##cxcl15
if(is.na(para_List[7])){message("请输出signature4:PMN,Mo,T");quit(save = "no")}
whichSignature=para_List[7];
#getGene0="CXCL17"
#getGene0="CXCL5"
#getGene0="PML"
ref_Gene="ZBTB7A"
ref_Gene1="PTEN"
k_n0=3
##读取数据ran_seq_median
path="E:/data/tcga/lihc_tcga/lihc_tcga/"
##path="D:\\min数据\\tcga\\prad_su2c_2015\\prad_su2c_2015\\"
#out_put_dir=paste(path,getGene0,whichSignature,"\\",sep="")
out_put_dir=paste("e:\\temp",getGene0,whichSignature,"\\",sep="\\")
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
##out_put_dir=paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep="")
setwd(out_put_dir)
shell.exec(path)
shell.exec(paste(out_put_dir,sep=""))


###开始进行聚类
if(!(grepl("Rterm",paste(commandArgs(),collapse = "_")))){
  fix(rna_seq_data);
}
PMN_geneSignatures=c("ANPEP","CD38","CSF3","CXCR2","CD80","IDO","SLEB2","CD274","PTPRC","TLR3","CEACAM8","STAT1","CXCR4","LOC116196","CSF1R","ITGAM","ITGAX","CD14","TGFB1","CXCL12","LOC199828","ENTPD1","FUT4","STAT3","IL4R","STAT5A","TLR4","CSF2","CXCL8","S100A8","S100A9","TNF")
#Mo_geneSignatures0=c("ARG1","CCL2","CD40","IL4R","FCGR2A","CD68","CD33","CSF1R","ITGAM","HLA-DRB1","FOXP3","PDCDILG2","CD163","LY75","PTPRC","CCR2","CD200R1","IL10","FCER2")
Mo_geneSignatures=c("ARG1","CCL2","CD40","IL4R","FCGR2A","CD68","LOC116196","CSF1R","ITGAM","HLA-DRB1","JM2","PDCD1LG2","CD163","LY75","PTPRC","LOC90262","CD200R1","IL10","FCER2")
T_geneSignatures=c("HLA-DOB","CXCL10","CXCL9","ICOS","CD8A","GZMK","HLA-DOA","HLA-DMB","CCL2","IRF1","HLA-DMA","CCL4","CCL3")
PI3K_AKT_mTOR_geneSignatures=c("PIK3CA","PIK3R1","PIK3R2","PTEN","PDPK1","AKT1","AKT2","FOXO1","FOXO3","FRAP2","RICTOR","TSC1","LOC124054","RHEB","AKT1S1","LOC654218","MLST8")
geneSignatures=eval(parse(text = paste(whichSignature,"_geneSignatures",sep = "")))

useCbioportal=FALSE
if(!useCbioportal){
  data2use=rna_seq_data[geneSignatures,];
}else{
  bigFile="G:\\TcgaTargetGtex_RSEM_Hugo_norm_count"
  message(paste("花费了",Sys.time()-begin_time_stamp,"s"))
}
rownames(data2use)=geneSignatures;if(!(grepl("Rterm",paste(commandArgs(),collapse = "_")))){fix(data2use)}
###完成了数据的提取，下面进行数据的预处理，分为两部分，复制出一部分进行data_plot
data_plot=data2use
##去除无用的列表
if("Entrez_Gene_Id" %in% colnames(data_plot)){
  data_plot=data_plot[,-which(colnames(data_plot)=="Entrez_Gene_Id")]  
}
if("Hugo_Symbol" %in% colnames(data_plot)){
  data_plot=data_plot[,-which(colnames(data_plot)=="Hugo_Symbol")]
};if(!(grepl("Rterm",paste(commandArgs(),collapse = "_")))){fix(data_plot)}
##进行数据的处理，例如log2,normalize
library(lattice)
library(ggpubr)
library(reshape)
#win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()
if(TRUE){
  {"进行log2处理！"}
  data_plot=log2(data_plot+1)
  ########对data_plot进行最终的可视化ggpurb
  a=melt(t(data_plot))
  colnames(a)=c("Sample","Gene","Value")
  win.graph();ggpubr::ggdensity(data=a,x="Value",color = "Gene",add = "mean")
}
if(FALSE){
  {"进行-rowmedian处理"}
  library(matrixStats)
  data_plot=data_plot-rowMedians(data.matrix(data_plot))
  ########对data_plot进行最终的可视化ggpurb
  a=melt(t(data_plot))
  colnames(a)=c("Sample","Gene","Value")
  win.graph();ggpubr::ggdensity(data=a,x="Value",color = "Gene",add = "mean")
}
if(TRUE){
  {"进行scale处理"}
  data_plot=data.frame(t(scale(t(data_plot)+1)))
  ########对data_plot进行最终的可视化ggpurb
  a=melt(t(data_plot))
  colnames(a)=c("Sample","Gene","Value")
  win.graph();ggpubr::ggdensity(data=a,x="Value",color = "Gene",add = "mean")
}
#win.graph();densityplot(unlist(data_plot));Sys.sleep(5);dev.off()
####开始进行热图的绘制
##计算break
goHeatmap=function(d1=data_plot,distmethod="euclidean",clustmethod="complete"){
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
  #设置colsidebar的颜色
  win.graph();
  h0=heatmap.2_adj(x=as.matrix(d1),Rowv = TRUE, Colv =TRUE, distfun=myDist,hclustfun=myClust)
  dev.off();
  ##主要的绘图函数
  #    使用 dendextend 包增强热图
  #   软件包 dendextend 可以用于增强其他软件包的功能
  ################
  library(dendextend)# order for rows
  library(magrittr)
  hcc<<-as.data.frame(cutree(h0$colDendrogram,k=5))
  colnames(hcc) <- "group_index_cut"  ##此时样本的顺序已经一致
  hcc[,"ori_rowindex"]<-c(1:ncol(data_plot))
  hcc[,"Sample"]=rownames(hcc)
  ##添加getgene0的表达值
  Value_getgene0_s=as.data.frame(t(rna_seq_data[getGene0,][,-c(1,2)]))
  Value_getgene0_s[,"Sample"]=rownames(Value_getgene0_s)
  hcc=merge(hcc,Value_getgene0_s,by.x = "Sample",by.y = "Sample",sort = FALSE)
  rownames(hcc)=hcc$Sample
  ##不断恢复顺序
  hcc=hcc[colnames(data_plot)[h0$colInd],]##变换为热图上的顺序
  
  
  ##分配cut_group_mean(0.25-0.75mean),热图区域的平均值
  for(onegroupindex in unique(hcc$group_index_cut)){
    sub_samples=rownames(hcc)[which(hcc$group_index_cut==onegroupindex)]
    sub_samples_data=data_plot[,sub_samples]##获取区块数据，使用作图时的数据计算平均值
    sub_samples_data_list=unlist(sub_samples_data)
    QS=quantile(sub_samples_data_list)##使用interquantile的数据的平均值
    sub_samples_data_interquantile_mean=mean(sub_samples_data_list[intersect(which(sub_samples_data_list>QS[2]),which(sub_samples_data_list<QS[4]))])##0.25~0.75作为平均值
    block_weight=sub_samples_data_interquantile_mean*10000##放大一万倍作为权重
    hcc[sub_samples,"interquantile_mean"]=sub_samples_data_interquantile_mean
    hcc[sub_samples,"block_weight"]=block_weight
  }
  ##hcc分流分配到作图
  hcc=hcc[colnames(data_plot),]##变换为原始数据的顺序，为datadriver不要重排
  go_wts=(hcc$block_weight)
  
  library(gplots);
  tiff(filename = paste(fig_order,paste("Heatmap","Distmethod",distmethod,"Clustermethod",clustmethod,sep="_"),".tiff",sep=""),width = 1000,height = 800)
  fig_order<<-fig_order+1
  hm<<-heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1), distfun=myDist,hclustfun=myClust,
                 scale="none",
                 trace="none",
                 dendrogram = "column",
                 Rowv = FALSE ,
                 ##########################################################这里一定要设置为mean不要设置为默认的sum
                 Colv =rev(reorder(as.dendrogram(h0$colDendrogram),wts = go_wts,agglo.FUN = (function(x){mean(x)}))%>%branches_color(k=5)),
                 #Colv = TRUE,
                 sepcolor = "white",
                 symkey = TRUE,
                 #breaks = c(min(d1)-0.05,seq(-3,3,0.005),max(d1)+0.05)
                 # breaks = c(min(d1)-0.5,seq(-3,3,0.005),max(d1)+0.5)
                 breaks = BREAKS,
                 main = paste("distmethod",distmethod,"clustermethod",clustmethod,sep="_")
                 #此参数仅仅为win.graph()设计，markdown不适用 
                 ####,lmat=rbind(c(0,0,3,0,0,0),c(0,0,4,4,4,0),c(0,0,1,1,1,0),c(0,0,2,2,2,0),c(0,5,5,5,0,0)) ,lhei=c(1.5,3,0.3,6,3),lwid=rbind(c(1,0.6,8,0.5,0.46,1))
                 #### ,margins = rep(4,2)
                 #,ColSideColors = colsidecolors
  )
  
  dev.off()   
  ###
  hcc=hcc[colnames(data_plot)[hm$colInd],]##变换为热图上的顺序
  data_boxplot_return=data.frame(matrix(nrow = 0,ncol = 3));colnames(data_boxplot_return)=list("group_index_cut","Value","plot_prio")
  prio=1;
  for(onegroupindex in unique(hcc$group_index_cut)){
    sample_onegroupindex=rownames(hcc)[which(hcc$group_index_cut==onegroupindex)]
    data_sample_onegroupindex=data_plot[,sample_onegroupindex]
    batchreturn_data_sample_onegroupindex=data.frame(group_index_cut=onegroupindex,Value=as.numeric(unlist(data_sample_onegroupindex)),plot_prio=paste(prio,onegroupindex,sep="_"))
    data_boxplot_return=rbind(data_boxplot_return,batchreturn_data_sample_onegroupindex)
    prio=prio+1
  }
  ##添加groupcount,不使用merge的方法
  for(onegroupindex in  unique(data_boxplot_return$group_index_cut)){
    data_boxplot_return[which(data_boxplot_return$group_index_cut==onegroupindex),"group_count"]=
      length(which(data_boxplot_return$group_index_cut==onegroupindex))/nrow(data_plot)
  }
  ##添加颜色
  cols=data.frame(cols=unique(dendextend::get_leaves_branches_col(hm$colDendrogram)),
                  group_index_cut=unique(data_boxplot_return$group_index_cut))
  data_boxplot_return=data_boxplot_return[order(data_boxplot_return$plot_prio),]
  data_boxplot_return=merge(data_boxplot_return,cols,by.x = "group_index_cut", by.y = "group_index_cut")
  
  ##merge后的恢复plot顺序
  data_boxplot_return=data_boxplot_return[order(data_boxplot_return$plot_prio),]
  #反正heatmap.2的图形输出不受控制这样就直接返回boxplot吧
  ###进行集合的比较
  data_boxplot_return$group_index_cut=as.character(data_boxplot_return$group_index_cut)
  my_comparision_matrix=as.data.frame(matrix(ncol=2,t(combn(x=unique(data_boxplot_return$plot_prio),2))))
  colnames(my_comparision_matrix)=c("first","second")
  for(onerow in c(1:nrow(my_comparision_matrix))){
    my_comparision_matrix[onerow,"differ"]=as.numeric(stringi::stri_sub(str=my_comparision_matrix[onerow,"second"],from = stringi::stri_length(my_comparision_matrix[onerow,"second"])))-
      as.numeric(stringi::stri_sub(str=my_comparision_matrix[onerow,"first"],from = stringi::stri_length(my_comparision_matrix[onerow,"second"])))
  }
  data_boxplot_return$plot_prio=as.factor(data_boxplot_return$plot_prio)
  ##按照差值进行排序
  my_comparision_matrix=my_comparision_matrix[order(my_comparision_matrix$first),]
  my_comparision_matrix=my_comparision_matrix[order(my_comparision_matrix$differ),]
  my_comparisons <- list()
  for(one in c(1:nrow(my_comparision_matrix))){
    oneterm=as.vector(c(as.character(my_comparision_matrix[one,"first"]),as.character(my_comparision_matrix[one,"second"])))
    my_comparisons[[one]]=oneterm
  }
  ###进行图形的绘制
  ###绘制图形qlot
  library(ggpubr);
  ##十分之一percent，用于绘图
  particle=max(diff(data_boxplot_return$Value))/9.5
  q=qplot(data = data_boxplot_return,x=data_boxplot_return$plot_prio,y=Value,geom = "boxplot",outlier.colour = "black",outlier.colour="black",
          ylab = "Value",main = paste("dist_method",distmethod,"hclust_method",clustmethod,sep = "_"))
  q= q+geom_jitter(colour = data_boxplot_return$cols)
  q=q+geom_text(data =data_boxplot_return,aes(label=paste("n=",group_count,sep=""),y=min(data_boxplot_return$Value)-particle))
  q=q+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),method = "t.test",
                         comparisons = my_comparisons,paired = FALSE,#label = "pb.format",
                         #hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                         label.y = c(seq(max(data_boxplot_return$Value),max(data_boxplot_return$Value)+max(diff(data_boxplot_return$Value))/10*10,particle))[c(1:length(my_comparisons))]
  )
  q=q+stat_compare_means(label.y = c(seq(max(data_boxplot_return$Value),max(data_boxplot_return$Value)+max(diff(data_boxplot_return$Value))/10*(length(my_comparisons)+1),particle))[length(my_comparisons)+1] )
  q=q+geom_hline(yintercept = mean(data_boxplot_return$Value), linetype=2)
  q=q+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(data_boxplot_return$Value)-particle/2)# Pairwise comparison against all
  #由于点数过多交换1,2的位置
  temp=q$layers[[1]];q$layers[[1]]=q$layers[[2]];q$layers[[2]]=temp
  tiff(filename = paste(fig_order,paste("Heatmap_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,sep="_"),".tiff",sep = ""),width = 1000,height = 800)         
  fig_order<<-fig_order+1
  plot(q)
  dev.off()
  
  if(TRUE){
    ######
    ######get_gene0      boxplot
    data_getgene0_boxplot=hcc[,c(getGene0,"group_index_cut")]
    data_getgene0_boxplot[,getGene0]=log2(data_getgene0_boxplot[,getGene0]+1)
    ##恢复顺序
    data_getgene0_boxplot=data_getgene0_boxplot[colnames(data_plot)[hm$colInd],]##变换为热图上的顺序
    ##添加plot顺序和groupcount
    prio=1;
    for(onegroupindex in unique(data_getgene0_boxplot$group_index_cut)){
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"plot_prio"]=
        paste(prio,onegroupindex,sep="_")
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"group_count"]=
        length(which(data_getgene0_boxplot$group_index_cut==onegroupindex))
      prio=prio+1
    }
    ##添加颜色
    cols=data.frame(cols=unique(dendextend::get_leaves_branches_col(hm$colDendrogram)),
                    group_index_cut=unique(data_getgene0_boxplot$group_index_cut))
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    data_getgene0_boxplot=merge(data_getgene0_boxplot,cols,by.x = "group_index_cut", by.y = "group_index_cut")
    ##merge后的恢复plot顺序
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    if(TRUE){
      ####进行组融合
      groups=unique(data_getgene0_boxplot$plot_prio)
      merges=list()
      for(onemerge in merges){
        print(onemerge)
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"plot_prio"]=
          groups[onemerge[[1]][1]]
        ##重新计算group_count
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"group_count"]=
          length(which(data_getgene0_boxplot$plot_prio==groups[onemerge[[1]][1]]))
      }
      ##重新计算my_comparision组合
      my_comparision_matrix=as.data.frame(matrix(ncol=2,t(combn(x=unique(data_getgene0_boxplot$plot_prio),2))))
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
    }
    ###开始绘图boxplot
    particle=max(diff(data_getgene0_boxplot[,getGene0]))/9.5
    q1=qplot(data = data_getgene0_boxplot,x=data_getgene0_boxplot$plot_prio,y=data_getgene0_boxplot[,getGene0],geom = "boxplot",outlier.colour = "black",outlier.colour="black",
             ylab = "Value",main = paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_")
    )
    q1= q1+geom_jitter(colour = data_getgene0_boxplot$cols)
    q1=q1+geom_text(data =data_getgene0_boxplot,aes(label=paste("n=",group_count,sep=""),y=min(data_getgene0_boxplot[,getGene0])-7.5*particle))
    q1=q1+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),method = "t.test",
                             comparisons = my_comparisons,paired = FALSE,label = "pb.format",
                             hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*10,particle))[c(1:length(my_comparisons))]
    )
    q1=q1+stat_compare_means(label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*(length(my_comparisons)+1),particle))[length(my_comparisons)+1] )
    q1=q1+geom_hline(yintercept = mean(data_getgene0_boxplot[,getGene0]), linetype=2)
    q1=q1+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(data_getgene0_boxplot[,getGene0])-particle/2)# Pairwise comparison against all
    #由于点数过多交换1,2的位置
    temp=q1$layers[[1]];q1$layers[[1]]=q1$layers[[2]];q1$layers[[2]]=temp
    tiff(filename = paste(fig_order,paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_"),".tiff",sep = ""),width = 1000,height = 800)         
    fig_order<<-fig_order+1
    plot(q1)
    dev.off()
  }###endif
  
  
  if(TRUE){
    ######
    ######get_gene0      boxplot
    data_getgene0_boxplot=hcc[,c(getGene0,"group_index_cut")]
    data_getgene0_boxplot[,getGene0]=log2(data_getgene0_boxplot[,getGene0]+1)
    ##恢复顺序
    data_getgene0_boxplot=data_getgene0_boxplot[colnames(data_plot)[hm$colInd],]##变换为热图上的顺序
    ##添加plot顺序和groupcount
    prio=1;
    for(onegroupindex in unique(data_getgene0_boxplot$group_index_cut)){
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"plot_prio"]=
        paste(prio,onegroupindex,sep="_")
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"group_count"]=
        length(which(data_getgene0_boxplot$group_index_cut==onegroupindex))
      prio=prio+1
    }
    ##添加颜色
    cols=data.frame(cols=unique(dendextend::get_leaves_branches_col(hm$colDendrogram)),
                    group_index_cut=unique(data_getgene0_boxplot$group_index_cut))
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    data_getgene0_boxplot=merge(data_getgene0_boxplot,cols,by.x = "group_index_cut", by.y = "group_index_cut")
    ##merge后的恢复plot顺序
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    if(TRUE){
      ####进行组融合
      groups=unique(data_getgene0_boxplot$plot_prio)
      merges=list(c(1,2),c(4,5))
      for(onemerge in merges){
        print(onemerge)
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"plot_prio"]=
          groups[onemerge[[1]][1]]
        ##重新计算group_count
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"group_count"]=
          length(which(data_getgene0_boxplot$plot_prio==groups[onemerge[[1]][1]]))
      }
      ##重新计算my_comparision组合
      my_comparision_matrix=as.data.frame(matrix(ncol=2,t(combn(x=unique(data_getgene0_boxplot$plot_prio),2))))
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
    }
    ###开始绘图boxplot
    particle=max(diff(data_getgene0_boxplot[,getGene0]))/9.5
    q1=qplot(data = data_getgene0_boxplot,x=data_getgene0_boxplot$plot_prio,y=data_getgene0_boxplot[,getGene0],geom = "boxplot",outlier.colour = "black",outlier.colour="black",
             ylab = "Value",main = paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_")
    )
    q1= q1+geom_jitter(colour = data_getgene0_boxplot$cols)
    q1=q1+geom_text(data =data_getgene0_boxplot,aes(label=paste("n=",group_count,sep=""),y=min(data_getgene0_boxplot[,getGene0])-7.5*particle))
    q1=q1+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),method = "t.test",
                             comparisons = my_comparisons,paired = FALSE,label = "pb.format",
                             hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*10,particle))[c(1:length(my_comparisons))]
    )
    q1=q1+stat_compare_means(label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*(length(my_comparisons)+1),particle))[length(my_comparisons)+1] )
    q1=q1+geom_hline(yintercept = mean(data_getgene0_boxplot[,getGene0]), linetype=2)
    q1=q1+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(data_getgene0_boxplot[,getGene0])-particle/2)# Pairwise comparison against all
    #由于点数过多交换1,2的位置
    temp=q1$layers[[1]];q1$layers[[1]]=q1$layers[[2]];q1$layers[[2]]=temp
    tiff(filename = paste(fig_order,paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_"),".tiff",sep = ""),width = 1000,height = 800)         
    fig_order<<-fig_order+1
    plot(q1)
    dev.off()
  }###1245endif
  
  if(TRUE){
    ######
    ######get_gene0      boxplot
    data_getgene0_boxplot=hcc[,c(getGene0,"group_index_cut")]
    data_getgene0_boxplot[,getGene0]=log2(data_getgene0_boxplot[,getGene0]+1)
    ##恢复顺序
    data_getgene0_boxplot=data_getgene0_boxplot[colnames(data_plot)[hm$colInd],]##变换为热图上的顺序
    ##添加plot顺序和groupcount
    prio=1;
    for(onegroupindex in unique(data_getgene0_boxplot$group_index_cut)){
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"plot_prio"]=
        paste(prio,onegroupindex,sep="_")
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"group_count"]=
        length(which(data_getgene0_boxplot$group_index_cut==onegroupindex))
      prio=prio+1
    }
    ##添加颜色
    cols=data.frame(cols=unique(dendextend::get_leaves_branches_col(hm$colDendrogram)),
                    group_index_cut=unique(data_getgene0_boxplot$group_index_cut))
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    data_getgene0_boxplot=merge(data_getgene0_boxplot,cols,by.x = "group_index_cut", by.y = "group_index_cut")
    ##merge后的恢复plot顺序
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    if(TRUE){
      ####进行组融合
      groups=unique(data_getgene0_boxplot$plot_prio)
      merges=list(c(1,2),c(3,4))
      for(onemerge in merges){
        print(onemerge)
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"plot_prio"]=
          groups[onemerge[[1]][1]]
        ##重新计算group_count
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"group_count"]=
          length(which(data_getgene0_boxplot$plot_prio==groups[onemerge[[1]][1]]))
      }
      ##重新计算my_comparision组合
      my_comparision_matrix=as.data.frame(matrix(ncol=2,t(combn(x=unique(data_getgene0_boxplot$plot_prio),2))))
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
    }
    ###开始绘图boxplot
    particle=max(diff(data_getgene0_boxplot[,getGene0]))/9.5
    q1=qplot(data = data_getgene0_boxplot,x=data_getgene0_boxplot$plot_prio,y=data_getgene0_boxplot[,getGene0],geom = "boxplot",outlier.colour = "black",outlier.colour="black",
             ylab = "Value",main = paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_")
    )
    q1= q1+geom_jitter(colour = data_getgene0_boxplot$cols)
    q1=q1+geom_text(data =data_getgene0_boxplot,aes(label=paste("n=",group_count,sep=""),y=min(data_getgene0_boxplot[,getGene0])-7.5*particle))
    q1=q1+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),method = "t.test",
                             comparisons = my_comparisons,paired = FALSE,label = "pb.format",
                             hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*10,particle))[c(1:length(my_comparisons))]
    )
    q1=q1+stat_compare_means(label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*(length(my_comparisons)+1),particle))[length(my_comparisons)+1] )
    q1=q1+geom_hline(yintercept = mean(data_getgene0_boxplot[,getGene0]), linetype=2)
    q1=q1+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(data_getgene0_boxplot[,getGene0])-particle/2)# Pairwise comparison against all
    #由于点数过多交换1,2的位置
    temp=q1$layers[[1]];q1$layers[[1]]=q1$layers[[2]];q1$layers[[2]]=temp
    tiff(filename = paste(fig_order,paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_"),".tiff",sep = ""),width = 1000,height = 800)         
    fig_order<<-fig_order+1
    plot(q1)
    dev.off()
  }###1234endif
  
  if(TRUE){
    ######
    ######get_gene0      boxplot
    data_getgene0_boxplot=hcc[,c(getGene0,"group_index_cut")]
    data_getgene0_boxplot[,getGene0]=log2(data_getgene0_boxplot[,getGene0]+1)
    ##恢复顺序
    data_getgene0_boxplot=data_getgene0_boxplot[colnames(data_plot)[hm$colInd],]##变换为热图上的顺序
    ##添加plot顺序和groupcount
    prio=1;
    for(onegroupindex in unique(data_getgene0_boxplot$group_index_cut)){
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"plot_prio"]=
        paste(prio,onegroupindex,sep="_")
      data_getgene0_boxplot[which(data_getgene0_boxplot$group_index_cut==onegroupindex),"group_count"]=
        length(which(data_getgene0_boxplot$group_index_cut==onegroupindex))
      prio=prio+1
    }
    ##添加颜色
    cols=data.frame(cols=unique(dendextend::get_leaves_branches_col(hm$colDendrogram)),
                    group_index_cut=unique(data_getgene0_boxplot$group_index_cut))
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    data_getgene0_boxplot=merge(data_getgene0_boxplot,cols,by.x = "group_index_cut", by.y = "group_index_cut")
    ##merge后的恢复plot顺序
    data_getgene0_boxplot=data_getgene0_boxplot[order(data_getgene0_boxplot$plot_prio),]
    if(TRUE){
      ####进行组融合
      groups=unique(data_getgene0_boxplot$plot_prio)
      merges=list(c(2,3,4))
      for(onemerge in merges){
        print(onemerge)
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"plot_prio"]=
          groups[onemerge[[1]][1]]
        ##重新计算group_count
        data_getgene0_boxplot[which(data_getgene0_boxplot$plot_prio %in% groups[onemerge]),"group_count"]=
          length(which(data_getgene0_boxplot$plot_prio==groups[onemerge[[1]][1]]))
      }
      ##重新计算my_comparision组合
      my_comparision_matrix=as.data.frame(matrix(ncol=2,t(combn(x=unique(data_getgene0_boxplot$plot_prio),2))))
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
    }
    ###开始绘图boxplot
    particle=max(diff(data_getgene0_boxplot[,getGene0]))/9.5
    q1=qplot(data = data_getgene0_boxplot,x=data_getgene0_boxplot$plot_prio,y=data_getgene0_boxplot[,getGene0],geom = "boxplot",outlier.colour = "black",outlier.colour="black",
             ylab = "Value",main =paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_")
    )
    q1= q1+geom_jitter(colour = data_getgene0_boxplot$cols)
    q1=q1+geom_text(data =data_getgene0_boxplot,aes(label=paste("n=",group_count,sep=""),y=min(data_getgene0_boxplot[,getGene0])-7.5*particle))
    q1=q1+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),method = "t.test",
                             comparisons = my_comparisons,paired = FALSE,label = "pb.format",
                             hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*10,particle))[c(1:length(my_comparisons))]
    )
    q1=q1+stat_compare_means(label.y = c(seq(max(data_getgene0_boxplot[,getGene0]),max(data_getgene0_boxplot[,getGene0])+max(diff(data_getgene0_boxplot[,getGene0]))/10*(length(my_comparisons)+1),particle))[length(my_comparisons)+1] )
    q1=q1+geom_hline(yintercept = mean(data_getgene0_boxplot[,getGene0]), linetype=2)
    q1=q1+stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y =min(data_getgene0_boxplot[,getGene0])-particle/2)# Pairwise comparison against all
    #由于点数过多交换1,2的位置
    temp=q1$layers[[1]];q1$layers[[1]]=q1$layers[[2]];q1$layers[[2]]=temp
    tiff(filename = paste(fig_order,paste(getGene0,"_blockplot","Distmethod",distmethod,"Clustermethod",clustmethod,"merged",paste(unlist(merges),collapse="_"),sep="_"),".tiff",sep = ""),width = 1000,height = 800)         
    fig_order<<-fig_order+1
    plot(q1)
    dev.off()
  }###endif
  
  
  return(hm)
}


##############开始循环餐数
Cluster_Method<-c( 
  "ward.D",
  "ward.D2",
  "single",
  "complete",
  "average",
  "mcquitty",
  "median",
  "centroid"
)
Dist_Methods<-  c("euclidean"
                  #, "maximum"
                  ,"manhattan"#, 
                  #"canberra", 
                  #"binary", 
                  #"minkowski",
                  #"pearson", "spearman","kendall"
)

for(a in dev.list()){
  if(class(dev.list())!="NULL"){
    dev.off()
  }
}



if(FALSE){
  for(onedistmethod in Dist_Methods){
    for(oneclustmethod in Cluster_Method){
      win.graph();
      goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod)
    }
  }}



goboxplot=function(d2=data_plot,heatmap_result,getGene0_=getGene0,one_dist_method,one_clust_method){
  ###??????????????????
  library(factoextra)
  #dist=get_dist(x=t(d2),method = one_dist_method)
  #hc=hclust(d=dist,method = one_clust_method)
  hcc=data.frame(cutree(tree = heatmap_result$colDendrogram,k=3));colnames(hcc)="treeIndex"
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
  for(one in unique(hcc$rank_group_index)){
    hcc[which(hcc$rank_group_index==one),"group_count"]=length(which(hcc$rank_group_index==one))
  }
  
  #fix(hcc);
  color_hcc=data.frame(rank_group_index=unique(hcc$rank_group_index),color=c("blue","white","red"))
  hcc=merge(hcc,color_hcc,by.x = "rank_group_index",by.y = "rank_group_index")
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
  #fix(hcc)
  q=qplot(data = hcc,x=rank_group_index,y=hcc$Value,geom = "boxplot",outlier.colour = "black",outlier.colour="black",
          ylab = paste("Log2(",getGene0_,")",sep=""),main = paste("dist_method",one_dist_method,"hclust_method",one_clust_method,sep = "_"))
  q=q+ geom_jitter(aes(colour = rank_low2high))
  q=q+geom_text(data = hcc,aes(label=paste("n=",group_count,sep=""),y=min(hcc$Value)-particle))
  q=q+stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)),method = "t.test",
                         comparisons = my_comparisons,paired = FALSE,#label = "pb.format",
                         #hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                         label.y = c(seq(max(hcc$Value),max(hcc$Value)+max(diff(hcc$Value))/10*10,particle))[c(1:3)]
  )
  q=q+stat_compare_means(label.y = c(seq(max(hcc$Value),max(hcc$Value)+max(diff(hcc$Value))/10*10,particle))[4] )
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
plot_together_rmd=function(){
  for(onedistmethod in Dist_Methods){
    for(oneclustmethod in Cluster_Method){
      h=goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod);
      #h自动输出热图，不可控
      q1<-goboxplot(d2=data_plot
                    ,getGene0_ = getGene0
                    ,one_dist_method = onedistmethod
                    ,one_clust_method = oneclustmethod
      );
      #print(q1)
      library(plotly)
      plotlyOutput(renderPlotly(ggplotly(q1)), width = "100%", height = "800px", inline = FALSE)
    }
  }
}



#######把图画在一起
plot_together=function(){
  for(onedistmethod in Dist_Methods){
    for(oneclustmethod in Cluster_Method){
      #win.graph();
      tryCatch({
        h=goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod);
      },error=function(e){cat("ERROR :",conditionMessage(e),"\n")});
      #h自动输出热图，不可控
      tryCatch({
        q=goboxplot(d2=data_plot
                    ,heatmap_result=h
                    ,getGene0_ = getGene0
                    ,one_dist_method = onedistmethod
                    ,one_clust_method = oneclustmethod
        );
        plot(q)
      },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
      #win.graph();
    }
  }
}
setwd(out_put_dir)
#plot_together()


par_fun_plot_together<-function(input_x){
  print(input_x);
  onedistmethod=input_x[1];oneclustmethod=input_x[2];
  #win.graph();
  tryCatch({
    h=goHeatmap(d1=data_plot,distmethod = onedistmethod,clustmethod = oneclustmethod);
  },error=function(e){cat("ERROR :",conditionMessage(e),"\n")});
  #h自动输出热图，不可控
  tryCatch({
    q=goboxplot(d2=data_plot
                ,heatmap_result=h
                ,getGene0_ = getGene0
                ,one_dist_method = onedistmethod
                ,one_clust_method = oneclustmethod
    );
    plot(q)
  },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  #win.graph();
}



  # combinations=merge(Dist_Methods,Cluster_Method)
  # apply(combinations,MARGIN = 1, FUN = par_fun_plot_together)
plot_together();  



source("E:/src/r/AGAIN_TCGA_PRAD/sendMail.R")
goMail(paste(getGene0,whichSignature,sep="_"))

