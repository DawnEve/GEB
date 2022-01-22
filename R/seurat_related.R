# Part I Seurat and ext ==================
# functions enhencing or relating to Seurat


#' convert human genes to mouse genes.
#'
#' @param s a vector of human gene symbols, like IFIT1
#'
#' @return a vector of mouse gene symbols, like Ifit1
#' @export
#'
#' @examples
#' convertH2M( c("IFIT1") )
convertH2M=function(s){
  paste0(substr(s, 1,1), tolower(substring(s, 2)) )
}

if(0){
  Mcc.gene =list(
    s.genes=convertH2M(Seurat::cc.genes.updated.2019$s.genes),
    g2m.genes=convertH2M(Seurat::cc.genes.updated.2019$g2m.genes)
  )
}


#' shortcut: HVG, Scale(all), PCA
#'
#' @param object Seurat obj
#' @param nfeatures number of HVG
#'
#' @return
#' @export
#'
#' @examples
#' Feature_scale_PCA_all(sce, 3000)
Feature_scale_PCA_all=function(object, nfeatures){
  object <- FindVariableFeatures(object = object, nfeatures=nfeatures)
  object <- ScaleData(object = object, features = rownames(object))
  object <- RunPCA(object=object, features=VariableFeatures(object))
  return(object)
}


#' shortcut: HVG, Scale(HVG), PCA
#'
#' @param object Seurat obj
#' @param nfeatures numbre of HVG
#'
#' @return
#' @export
#'
#' @examples
#' Feature_scale_PCA_variable(sce, 3000)
Feature_scale_PCA_variable = function(object, nfeatures){
  object <- FindVariableFeatures(object = object,nfeatures=nfeatures)
  object <- ScaleData(object = object,features = VariableFeatures(object))
  object <- RunPCA(object=object,features=VariableFeatures(object))
  return(object)
}




#' shortcut: cell Cluster, UMAP
#'
#' @param object Seurat obj
#' @param dims dimentions to use of PC
#' @param resolution a number between 0.1-2, sometimes evern 50, 200
#'
#' @return
#' @export
#'
#' @examples
#' Neighbor_cluster_umap_tsne(sce, 1:30, 0.8)
Neighbor_cluster_umap_tsne = function(object,dims,resolution){
  object <- FindNeighbors(object = object,dims=dims)
  object <- FindClusters(object = object, resolution = resolution)
  object <- RunUMAP(object = object, dims=dims)
  return(object)
}



#' FindMarkersDIY, quick version
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
#' FindMarkersDIY(sce)
FindMarkersDIY = function(object = object){
  l <- length(levels(object))
  idents <- levels(object)
  temp <- c()
  for (i in 1:l) {
    temp[[i]] <- apply(
      X = object@assays$RNA@data[, WhichCells(object, idents = idents[i]), drop = FALSE],
      MARGIN = 1,
      FUN = function(x){log(x = mean(x = expm1(x = x)) + 1)}
    )
  }
  temp.exp <- as.data.frame(temp, col.names = idents)
  for (i in 1:l) {
    temp[[i]] <- temp.exp[,i] - log(x = rowMeans(x = expm1(x = temp.exp[, -c(i)])) + 1)
  }
  temp.exp.diff <- as.data.frame(temp, col.names = idents)
  temp.exp.diff$gene <- rownames(temp.exp.diff)
  temp.exp.diff <- reshape2::melt(temp.exp.diff, id.vars = "gene", variable.name = "cluster", value.name = "avg_logFC")
  return(temp.exp.diff)
}



#' 2D FACS like scatter plot
#'
#' @param object Seurat ojb
#' @param feature1 gene symbol 1
#' @param feature2 gene symbol 2
#' @param cells cid
#' @param group.by group by
#' @param cols colors
#' @param pt.size dot size
#' @param shape.by shape by
#' @param span span?
#' @param smooth smooth?
#' @param combine combine?
#' @param slot slot
#' @param plot.cor plot.correlation?
#' @param raster is raster
#'
#' @return
#' @export
#'
#' @examples
#' nFeatureScatter(sce, "CD3D", "CD4")
nFeatureScatter=function(
  object,
  feature1,
  feature2,
  cells = NULL,
  group.by = NULL,
  cols = NULL,
  pt.size = 1,
  shape.by = NULL,
  span = NULL,
  smooth = FALSE,
  combine = TRUE,
  slot = 'data',
  plot.cor = TRUE,
  raster = NULL
) {
  cells <- cells %||% colnames(x = object)
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  data <-  FetchData(
    object = object,
    vars = c(feature1, feature2, group.by),
    cells = cells,
    slot = slot
  )
  if (!grepl(pattern = feature1, x = colnames(x = data)[1])) {
    stop("Feature 1 (", feature1, ") not found.", call. = FALSE)
  }
  if (!grepl(pattern = feature2, x = colnames(x = data)[2])) {
    stop("Feature 2 (", feature2, ") not found.", call. = FALSE)
  }
  data <- as.data.frame(x = data)
  feature1 <-  colnames(x = data)[1]
  feature2 <-  colnames(x = data)[2]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      SingleCorPlot(
        data = data[,c(feature1, feature2)],
        col.by = data[, x],
        cols = cols,
        pt.size = pt.size,
        smooth = smooth,
        legend.title = 'Identity',
        span = span,
        plot.cor = plot.cor,
        raster = raster
      )
    }
  )
  if (isTRUE(x = length(x = plots) == 1)) {
    return(plots[[1]])
  }
  if (isTRUE(x = combine)) {
    plots <- wrap_plots(plots, ncol = length(x = group.by))
  }
  return(plots)
}





#' SingleCorPlot
#'
#' @param data matrix?
#' @param col.by color by?
#' @param cols colors?
#' @param pt.size dot size
#' @param smooth smooth?
#' @param rows.highlight which row?
#' @param legend.title legend title
#' @param na.value value of NA
#' @param span span?
#' @param raster is raster?
#' @param plot.cor plot correlation
#'
#' @return
#' @export
#'
#' @examples
#' SingleCorPlot(sce)
SingleCorPlot = function(
  data,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  smooth = FALSE,
  rows.highlight = NULL,
  legend.title = NULL,
  na.value = 'grey50',
  span = NULL,
  raster = NULL,
  plot.cor = TRUE
) {
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  if ((nrow(x = data) > 1e5) & !isFALSE(raster)){
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e5)
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  names.plot <- colnames(x = data) <- gsub(
    pattern = ':',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  if (ncol(x = data) < 2) {
    msg <- "Too few variables passed"
    if (ncol(x = data) == 1) {
      msg <- paste0(msg, ', only have ', colnames(x = data)[1])
    }
    stop(msg, call. = FALSE)
  }
  plot.cor <- if (isTRUE(x = plot.cor)) {
    round(x = cor(x = data[, 1], y = data[, 2]), digits = 2)
  }
  else(
    ""
  )
  if (!is.null(x = rows.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = rows.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size,
      cols.highlight = 'red',
      col.base = 'black',
      pt.size = pt.size
    )
    cols <- highlight.info$color
    col.by <- factor(
      x = highlight.info$highlight,
      levels = rev(x = highlight.info$plot.order)
    )
    plot.order <- order(col.by)
    data <- data[plot.order, ]
    col.by <- col.by[plot.order]
  }
  if (!is.null(x = col.by)) {
    data$colors <- col.by
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = names.plot[1], y = names.plot[2])
  ) +
    labs(
      x = orig.names[1],
      y = orig.names[2],
      title = plot.cor,
      color = legend.title
    )
  if (smooth) {
    # density <- kde2d(x = data[, names.plot[1]], y = data[, names.plot[2]], h = Bandwidth(data = data[, names.plot]), n = 200)
    # density <- data.frame(
    #   expand.grid(
    #     x = density$x,
    #     y = density$y
    #   ),
    #   density = as.vector(x = density$z)
    # )
    plot <- plot + stat_density2d(
      mapping = aes(fill = ..density.. ^ 0.25),
      geom = 'tile',
      contour = FALSE,
      n = 200,
      h = Bandwidth(data = data[, names.plot])
    ) +
      # geom_tile(
      #   mapping = aes_string(
      #     x = 'x',
      #     y = 'y',
      #     fill = 'density'
      #   ),
      #   data = density
      # ) +
      scale_fill_continuous(low = 'white', high = 'dodgerblue4') +
      guides(fill = FALSE)
  }
  if (!is.null(x = col.by)) {
    if (raster) {
      plot <- plot + geom_scattermore(
        mapping = aes_string(color = 'colors'),
        position = 'jitter',
        pointsize = pt.size
      )
    } else {
      plot <- plot + geom_point(
        mapping = aes_string(color = 'colors'),
        position = 'jitter',
        size = pt.size
      )
    }
  } else {
    if (raster) {
      plot <- plot + geom_scattermore(position = 'jitter', pointsize = pt.size)
    } else {
      plot <- plot + geom_point(position = 'jitter', size = pt.size)
    }
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (length(x = cols) == 1 && cols %in% rownames(x = brewer.pal.info)) {
      scale_color_brewer(palette = cols)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + cols.scale
    if (!is.null(x = rows.highlight)) {
      plot <- plot + guides(color = FALSE)
    }
  }
  plot <- plot + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(x = span)) {
    plot <- plot + geom_smooth(
      mapping = aes_string(x = names.plot[1], y = names.plot[2]),
      method = 'loess',
      span = span
    )
  }
  return(plot)
}



#' Show nGene hist plot in each cluster
#'
#' v2.0 order output
#'
#' @param object Seurat obj
#' @param nth a number
#' @param draw wether output text of draw picture
#'
#' @return
#' @export
#'
#' @examples
#' rs=meanNFeatureByCluster(sce, 3)
meanNFeatureByCluster = function(object, nth, draw=T ){
  df0=sapply( split(object@meta.data$nFeature_RNA, object@meta.data$seurat_clusters), mean)
  df1=data.frame( mean=as.numeric(df0), cluster=as.character(names(df0)) )

  if(T==draw){
    oPar=par(no.readonly=T)
    par(mfrow=c(2,1), mar=c(4,4,1,1))
    plot(sort(df1$mean), ylab="mean of nFeature", mgp=c(2, 0.5, 0))
    abline(v=nth, col='red', lty=2)

    plot( log10(sort(df1$mean)), ylab="mean of nFeature(log10)", mgp=c(2, 0.5, 0) )
    abline(v=nth, col='red', lty=2)
    par(oPar)
  }

  cat( "mean(nFeature), nth, (n+1)th, median: ",
       c(sort(df1$mean)[nth], sort(df1$mean)[nth+1], median(df1$mean)), "\n" )
  df1=df1[order(df1$mean),]
  cat( "remove:",  df1[which(df1$mean<= df1$mean[nth] ), "cluster"], "\n")
  cat( "keep:",  df1[which(df1$mean> df1$mean[nth] ), "cluster"], "\n")
  #cat( "left:", setdiff( 0:(nrow(df1)-1), df1[which(df1$mean< df1$mean[nth] ), "cluster"]), "\n" )
  return(df1)
}





#' draw barplot, from a data.frame generated by table(para1, para2), colored by 1st parameter of table()
#'
#' v2.1
#' v2.2 第一参数的空值跳过
#'
#' @param tbl1 data.frame, draw by each column
#' @param colors colors of each row(default NULL, auto-color)
#' @param scale whether scale to 1 or not(default T)
#' @param title the main title of the figure,(is main in the function)
#' @param legendY the position y of legend, adjust as needed, default -0.25
#' @param omit whether to omit some columns
#' @param ...
#'
#' @return NULL
#' @export
#'
#' @examples
#' table2barplot(table(sce$seurat_clusters, sce$sample) )
table2barplot=function(tbl1, colors=NULL,levels=NULL, scale=T, title="",
                       omit=NULL, xlab="", ylab="", legendTitle=NULL){
  tbl1= tbl1[, which(colSums(tbl1)>0)] #remove all 0 columns
  tbl1= tbl1[which(rowSums(tbl1)>0), ] #remove all 0 rows
  # remove some columns by column names
  if(!is.null(omit)){
    tbl1=tbl1[, setdiff( colnames(tbl1), as.character( omit ) ) ]
  }
  # table to data.frame(wide to long)
  df2=as.data.frame(tbl1)
  if(!is.null(levels)){
    df2$Var1=factor(df2$Var1, levels =levels )#change order
  }
  if(""==ylab){ ylab=ifelse(scale, "Freq", "Count") }
  if(""==xlab){ xlab="Index" }

  # draw
  g1=ggplot(df2, aes(x=Var2, y=Freq, fill=Var1))+
    geom_bar(stat="identity", position=ifelse(scale, "fill","stack") )+
    labs(x=xlab, y=ylab, title=title)+
    theme_classic(base_size = 14)+
    theme(axis.text.x=element_text(angle=60, hjust=1,size=rel(1.2)) )
  if(is.null(colors)){
    return (g1 + scale_fill_discrete( ifelse(is.null(legendTitle), "", legendTitle) ) )
  }else{
    return( g1+scale_fill_manual(legendTitle, values = colors) )
  }
}
# debug
if(0){
  temp=table( mtcars$gear, substring( rownames(mtcars), 1,3))
  table2barplot( temp )
  table2barplot( temp, scale=F )
  table2barplot( temp, xlab="Marker", title="xxx" )
  table2barplot( temp, xlab="Marker", legendTitle = "Gear" )
  table2barplot( temp, colors=c("#00AFBB", "#E7B800", "#FC4E07") )
  #
  table2barplot( temp, colors=c("purple", "cyan", "yellow") )
  table2barplot( temp, colors=c("purple", "cyan", "yellow"), levels=c(5,4,3) )
  table2barplot( temp, colors=c("purple", "cyan", "yellow"), levels=c(5,4,3), title="XX" )
}





#' hist plot of nFeature in each cluster
#'
#' @param object Seurat obj
#' @param cluster which cluster
#' @param xlim limits of x axis
#'
#' @return
#' @export
#'
#' @examples
#' nFeature.hist(sce, 2, c(0,200))
nFeature.hist=function(object, cluster, xlim=c(0, 500)){
  params=list(
    subset(object, idents = cluster)@meta.data$nFeature_RNA,
    n=200,
    mgp=c(2,0.5,0),
    main=paste0("cluster=",cluster),
    xlab="nFeature_RNA"
  )
  if(!is.null(xlim))
    params$xlim=xlim;

  do.call(hist, params)
}






#' Show the position of a cluster on umap
#'
#' @param scObj Seurat obj
#' @param clusterId which cluster
#' @param slot default seurat_clusters
#' @param color highlight color, default darkred
#'
#' @return
#' @export
#'
#' @examples
#' showCluster(scObj, 17)
#' showCluster(scObj, "BL1Y", slot="sample")
#' showCluster(scObj, "BL1Y", slot="sample", color="purple")
showCluster=function(scObj, clusterId, slot="seurat_clusters", color="darkred"){
  df2 = as.data.frame(scObj@meta.data);
  DimPlot(scObj, label = F, #group.by = "sample",
          #cells.highlight = df2[which(df2[, slot] == clusterId), "cell"],
          cells.highlight = rownames( df2[which(df2[, slot] == clusterId), ]),
          cols.highlight = color, cols = "grey")+ #c("darkred", "darkblue")
    labs(title=paste0(slot, ": ", clusterId))+
    theme(
      legend.position = "none"
    )
}




#' get 3col-df (exp value/ gene symbol /  cluster id) for stack violin plot
#'
#' @param scRNA Seurat obj
#' @param genenames gene symbol list
#' @param group.by use which cluster cretia
#'
#' @return
#' @export
#'
#' @examples
#' df2= getLongDf(scObj_known2, c("Cd3e", "Cd19"), group.by='sample')
getLongDf=function(scRNA, genenames, group.by="seurat_clusters"){
  df1=NULL;
  meta=scRNA@meta.data
  level_list = unique(as.data.frame(scRNA@meta.data)[, group.by]);
  for( cluster in level_list){
    message(cluster)
    cid=rownames( meta[which(meta[, group.by]== cluster),] )
    mtx=scRNA@assays$RNA@data[, cid]
    for( gene in genenames){
      df2=data.frame(value=as.numeric( mtx[gene, ]) )
      df2$gene=gene
      df2$cluster=cluster
      df1=rbind(df1, df2)
    }
  }
  return(df1);
}





#' Plot stacked violin plot using 3col-df(value/ gene/  cluster)
#'
#' @param df1 df of 3 columns
#' @param type h or v
#' @param colors color list
#' @param fill_by like group by
#'
#' @return
#' @export
#'
#' @examples
#' VlnPlot_stack(df1, "h")
VlnPlot_stack=function(df1, type="h", colors=NULL, fill_by="cluster"){
  # 分类变量 to factor
  df1$cluster=factor(df1$cluster)
  if( is.na(match(fill_by, colnames(df1))) ){
    stop("Error: fill_by must be a column name of df1!")
  }

  if(is.null(colors)){
    print(colors)
    colors=rainbow( length( unique( df1[,fill_by])) )
  }

  if("h"==type){
    # VlnPlot 2: x=gene, y=ident, violin=horizotal
    p1 = ggplot(df1,aes(x=cluster, y=value, fill=df1[,fill_by]))+
      geom_violin(scale = "width",colour="white",alpha=0.85,width=1) +
      coord_flip() + guides(fill="none")+
      facet_wrap(gene~.,nrow = 1, strip.position = "bottom", scale="free_x") +
      theme_bw(base_size = 14)+labs(x="", y="")+
      scale_fill_manual(values=colors)+ #指定颜色
      theme(
        panel.grid = element_blank(), #不要背景网格
        axis.text.x = element_blank(), #不要x坐标轴刻度文字

        #axis.text.y = element_text(size=15), #y坐标刻度字号
        #axis.title.y = element_text(size=15), #y标题字号
        axis.ticks.x =  element_blank(), #不要x轴刻度线

        legend.position = "none", #不要图例
        panel.spacing=unit(0,"cm"), #分面的间距

        strip.placement = "outside", #分面标签位置
        strip.text.x = element_text(angle=-90,vjust=0.5, hjust = 0), #分面标签 文字倾斜
        strip.background = element_blank() #分面标签 不要背景
      )
  }else{
    # VlnPlot 1： x=gene, y=ident, violin=vertical
    p1=ggplot(df1, aes(gene, value, fill=df1[,fill_by]))+
      geom_violin(scale="width", color="black")+
      facet_wrap(cluster~. ,strip.position = "left", ncol=1, scales = "free_y") +
      theme_bw(base_size = 14)+labs(x="", y="")+
      scale_fill_manual(values=colors)+ #指定颜色
      theme(panel.grid = element_blank(), #不要背景
            axis.ticks.y=element_blank(), #不要y坐标刻度
            axis.text.x = element_text(angle=45,vjust=1,hjust = 1), #x文字倾斜45度
            axis.text.y = element_blank(), #不要y轴文字
            legend.position = "none", #不要图例
            panel.spacing=unit(0,"cm"), #分面间距为0

            strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5), #很神奇，没见过的
            strip.background = element_blank())
  }
  return(p1)
}






#' Plot stack violin plot from Seurat object
#'
#' @param object Seurat obj
#' @param features gene symbols
#' @param cols color list
#' @param group.by group by
#' @param type which type, h or v
#' @param fill_by like group by
#'
#' @return
#' @export
#'
#' @examples
#' genes2=c("Cd79a", "Igha", "Cd3d", "Itgax", "Ly6g");
#' colors2=c("red", "blue", "red3", "orange", "purple");
#' DoVlnPlot_stack(scObj_known2, features = genes2, cols =colors2,
#'   group.by = "ident", fill_by = "gene", type="v")
DoVlnPlot_stack=function(object, features, cols = NULL, group.by = NULL, type="h", fill_by="cluster"){
  params=list(object, features)
  if(!is.null(group.by)){
    params$slot=group.by;
  }
  df2= do.call(getLongDf, params)
  df2$gene=factor(df2$gene, levels = features)

  params=list(df2)
  return( VlnPlot_stack(df2, type=type, colors=cols, fill_by=fill_by) )
}







#' Draw UMAP, my method 2
#'
#' @param scObj Seurat obj
#' @param group.by group by
#'
#' @return
#' @export
#'
#' @examples
#' DimPlot_2(scObj_known2, "ident")
DimPlot_my2 = function(scObj, group.by="seurat_clusters", cols=NULL){
  # (1)为了调用ggplot2我们把UMAP的坐标放到metadata中：
  pbmc<-AddMetaData(scObj,scObj@reductions$umap@cell.embeddings,col.name = colnames(scObj@reductions$umap@cell.embeddings))
  #head(pbmc@meta.data)

  # (2)读入一套我珍藏多年的颜色列表：
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

  # (3) 开始画图
  df1= as.data.frame(pbmc@meta.data) #[, c("UMAP_1", "UMAP_2", group.by)]
  df1$cluster=as.data.frame(pbmc@meta.data)[, group.by]
  # colnames(df1)=c("UMAP_1", "UMAP_2","cluster")

  class_avg <-  df1%>%
    group_by(cluster) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  #print(head(class_avg))
  umap <-  ggplot(df1, aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=cluster))+
    scale_color_manual(values = allcolour)+
    geom_text(aes(label = cluster), data = class_avg)+
    theme(text=element_text(family="Arial",size=18)) +
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(), axis.title = element_text(color='black',
                                                                family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          #axis.ticks.margin = unit(0.6,"lines"),
          #`axis.ticks.margin` is deprecated. Please set `margin` property of `axis.text` instead

          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=18),
          axis.title.y=element_text(colour='black', size=18),
          axis.text=element_text(colour='black',size=18),
          legend.title=element_blank(),
          legend.text=element_text(family="Arial", size=18),
          legend.key=element_blank())+
    #theme(plot.title = element_text(size=22,colour = "black",face = "bold"))  +
    guides(colour = guide_legend(override.aes = list(size=5)))
  return(umap)
}







#' 快速桑基图，预览几种分类方法的来源和去向
#'
#' @param scObj0 Seurat obj
#' @param col_list column list of meta data
#' @param downsample sample size of each cluster, to get a small subset
#' @param color_res color by resolution, default 1
#'
#' @return
#' @export
#'
#' @examples
#' Q_sankey(pbmc)
#' Q_sankey(pbmc, c(5, 7))
Q_sankey=function(scObj0, col_list=NULL, downsample=NULL, color_res=1){
  library(ggforce)
  if(!is.null(downsample)){
    scObj=subset( scObj0, downsample=100)
  }else{
    scObj=scObj0;
  }

  if(is.null(col_list)){
    col_list=grep("RNA_snn", colnames(scObj@meta.data))
  }

  scObj@meta.data %>%
    gather_set_data( col_list ) %>% #选取要看的列 RNA_snn_res....
    ggplot(aes(x, id = id, split = y, value = 1))  +
    geom_parallel_sets(aes_string(fill = paste0("RNA_snn_res.", color_res )), show.legend = FALSE, alpha = 0.3) +
    geom_parallel_sets_axes(axis.width = 0.1, color = "lightgrey", fill = "white") +
    geom_parallel_sets_labels(angle = 0) +
    theme_no_axes()
}





#' 每个亚群中的细胞数量饼图
#'
#' @param scObj Seurat obj
#' @param group.by group by
#' @param colors color list
#'
#' @return
#' @export
#'
#' @examples
#' Q_cellNumber(sce)
#' Q_cellNumber(scObj_known2, colors=allcolour)
#' Q_cellNumber(scObj_known2, "ident", colors=allcolour)
#'
Q_cellNumber=function(scObj, group.by="seurat_clusters", colors=NULL){
  library(ggforce)
  g1=data.frame(cluster= as.data.frame(scObj@meta.data)[, group.by]) %>%
    count( cluster ) %>%
    #arrange(-n)  %>%
    ggplot() +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1,
                     amount = n, fill = cluster), alpha = 1, stat = "pie")+
    theme_no_axes()
  if(!is.null(colors))
    g1=g1+scale_fill_manual(values=allcolour)
  return(g1)
}






#' plot gene expression level with color, on nGene x nUMI plot
#'
#' @param scObj Seurat obj
#' @param gene_symbol gene symbol
#' @param keyword sample info on title
#'
#' @return
#' @export
#'
#' @examples
#' QC_nGene_nUMI_expr(sce)
#' QC_nGene_nUMI_expr(sce, "CD4")
QC_nGene_nUMI_expr=function(scObj, gene_symbol=NULL, keyword=""){
  metadata=scObj@meta.data
  if(""==keyword){
    keyword=scObj@project.name
  }
  # no gene symbol
  if(is.null(gene_symbol)){
    g1= ggplot(metadata, aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point(size=0.3, alpha=0.3) +
      labs(title=keyword)
  }else{
    # with gene symbol
    #gene_symbol="Ly6g"
    if(!(gene_symbol %in% rownames(scObj))){
      return(NULL)
    }
    g1=metadata %>%
      cbind( FetchData(scObj, gene_symbol) ) %>%
      ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=get(gene_symbol))) +
      geom_point(size=0.3) +
      scale_colour_gradient(name= gene_symbol, low = "#FF000000", high = "#FF0000FF") +
      # stat_smooth(method=lm, linetype=2, color="#aaaaaa55") +
      labs(title= ifelse(""==keyword, gene_symbol, sprintf("%s in sample %s", gene_symbol, keyword) ) )
  }

  # return
  g1+ scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 200, linetype=2, color="grey") +
    geom_hline(yintercept = c(100,200), linetype=2, color="grey")
}









#' heatmap of average expression, by cluster
#'
#' @param sce Seurat obj
#' @param features gene list
#' @param slot which expression type
#' @param assay assay
#' @param main title
#'
#' @return
#' @export
#'
#' @examples
#' DoHeatmap_my(scObj_colon, gene_203307, group.by = 'sample')
#' DoHeatmap_my(scObj_colon, gene_203307, group.by = 'sample', rotate=T)
DoHeatmap_my=function(sce, features, slot="scale.data", assay="RNA", main="",
                      group.by="seurat_clusters", rotate=F){
  #features=gene_203307
  #sce=scObj_colon
  #cut_dims=c(-0.5, 0.5)
  # 1. get df from Seurat: scaled gene expression
  .data = GetAssayData(	object = sce, assay = assay, slot = slot)[features,]
  .group=as.character( as.data.frame(sce@meta.data)[, group.by] )
  #dim(.data)
  #2. get mean by sample
  .data2=sapply( split( as.data.frame(t(.data)), .group), function(x){
    colMeans(x)
  })
  # remove all 0 rows
  .data2=.data2[apply(.data2, 1, sd)>0,]

  #3. clean large outliers
  #.data2[.data2<cut_dims[1]] <-cut_dims[1]
  #.data2[.data2>cut_dims[2]] <-cut_dims[2]
  if(0){ # how is this fig?
    tmp=AverageExpression(sce, slot="data", group.by = group.by)$RNA[features, ]
    dim(tmp)
    pheatmap(log(tmp+1), border_color = NA,
             clustering_method = "ward.D2",
             scale = "row",
             main="try4")
  }

  params2=list(mat=.data2,
               border_color = NA,
               clustering_method = "ward.D2",
               scale = "row",
               main=main)
  if(rotate){
    params2$mat=t(.data2)
    params2$scale="column"
  }
  library(pheatmap)
  do.call(pheatmap::pheatmap, params2)
}

