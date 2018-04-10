print("加载自定义函数包")

print("自定义聚类函数")
myclust<-function(x,clustmethod){
  hclust(x,method=clustmethod)
}


print("自定义距离函数")
myDist<-function(x,distmethod){
  get_dist(x,method = distmethod)
}


# Function to average the expression of probesets which map to same gene
probeset2genelevel = function (x){
  return(tapply(x, factor(x[,"Hugo_Symbol"], mean)))
}


print("热图函数")
library("RColorBrewer")
library(gplots)
goHeatmap=function(d1,clustmethod,distmethod){
  
  ##计算break
  PERCENTILE=0.005;lowQ=as.numeric(quantile(unlist(d1),PERCENTILE));highQ=as.numeric(quantile(unlist(d1),1-PERCENTILE))
  BREAKS=c(min(d1)-0.5,seq(lowQ,highQ,0.005),max(d1)+0.5)
  ##主要的绘图函数
  heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1), distfun=myDist,hclustfun=myclust,scale="none",
            trace="none",
            dendrogram = "column",
            Rowv = FALSE,
            Colv = TRUE,
            sepcolor = "white",
            symkey = TRUE,
            #breaks = c(min(d1)-0.05,seq(-3,3,0.005),max(d1)+0.05)
            # breaks = c(min(d1)-0.5,seq(-3,3,0.005),max(d1)+0.5)
            breaks = BREAKS
  )
}