


#' enhenced FeaturePlot, with one unified legend color-bar.
#'
#' @param sce Seurat obj
#' @param features gene symbol list
#' @param ncol num of columns
#' @param cut_dims the min,max of scaled data
#' @param colors color list, default will be ok
#'
#' @return
#' @export
#'
#' @examples
#' FeaturePlot_my(scObj_colon, features =c("Itgam", "Ly6g", "Cxcr2", "Cxcr4", "Cxcl10"), ncol=5 )
FeaturePlot_my=function(sce, features, ncol=5, cut_dims=c(-1,1.5), colors=brewer.pal(9, "YlOrRd")){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  # 1. get df from Seurat, umap coord and scaled gene expression

  umap_data = GetAssayData(	object = sce, assay = "RNA", slot = "scale.data")[features,]

  umap_data = sce@reductions$umap@cell.embeddings %>%  #UMAP coord
    as.data.frame() %>%
    #cbind( t(as.data.frame( sce@assays$RNA@scale.data)[features,]) ) # slow
    cbind( t(GetAssayData(	object = sce, assay = "RNA", slot = "scale.data")[features,]) ) # fast
  # 2. wide to long
  umap_data_long=umap_data %>% pivot_longer(cols=all_of(features) )
  hist(umap_data_long$value, n=100)
  # 3. cut if too high and too low
  if(!(is.null(cut_dims))){
    umap_data_long$value=ifelse(umap_data_long$value < cut_dims[1], cut_dims[1], umap_data_long$value)
    umap_data_long$value=ifelse(umap_data_long$value > cut_dims[2], cut_dims[2], umap_data_long$value)
  }
  # 4. to factor
  umap_data_long$name=factor(umap_data_long$name, levels = features)

  # 5. plot
  g1=ggplot(umap_data_long, aes(UMAP_1, UMAP_2, color=value) )+
    geom_point(size=0.3)+
    facet_wrap(name~.)+
    scale_color_gradientn(name="Z Score", colors=colors)+
    theme_void(base_size = 12)+
    theme(
      legend.position = c("left"),
      legend.text=element_text(size=8),
      legend.title=element_text(size=8),
      legend.key.width = unit(4, "mm"),
      legend.key.height = unit(4, "mm"),
      legend.margin=margin(t = 0, r = 0, b = 60, l = 0, unit = "pt"),
    );
  return(g1)
}






#' Calc gene score
#'
#' @param sce Seurat obj
#' @param gene_list gene list
#' @param group.by group by
#'
#' @return
#' @export
#'
#' @examples
#' getGeneScoreValue(sce, gene_list)
getGeneScoreValue=function(sce, gene_list, group.by=NULL){
  if(is.null(group.by)){
    group.by="seurat_clusters"
  }
  #1. calc score
  sce = AddModuleScore(sce, features = list("g1"=gene_list), name="score")
  return(data.frame(
    cluster=as.data.frame(sce@meta.data)[, group.by],
    score=sce$score1
  ))
}


#' plot violin plot of gene score among groups(ggplot2)
#'
#' @param sce Seurat obj
#' @param gene_list gene list
#' @param group.by group by
#'
#' @return
#' @export
#'
#' @examples
#' plotGeneScore(sce, gene_list)
plotGeneScore=function(sce, gene_list, group.by=NULL){
  df1=getGeneScoreValue(sce, gene_list, group.by )
  # head(df1)
  g1=ggplot(df1, aes(cluster, score, fill=cluster))+
    geom_violin()+
    theme_classic(base_size = 14)+
    theme(
      legend.position = "none"
    )+labs(x="")
  return(g1)
}





#' Draw UMAP or tSNE plot with circle backgrounded text
#'
#' @param sce Seurat obj
#' @param group.by group by
#' @param reduction reduction method
#' @param isRaster is raster?
#' @param show.legend default not show legend
#' @param theme.use which theme
#' @param par.geom_point geom_point parameters
#' @param par.geom_text geom_text parameters
#' @param par.circle circle parameters
#' @param title main title of the figure
#' @param colSet color list
#'
#' @return
#' @export
#'
#' @examples
#' DimPlot_my(sce, title="pbmc_3k") #default umap
#' DimPlot_my(sce, reduction = "tsne", title="pbmc_3k tSNE") #tsne
#' DimPlot_my(sce, reduction = "pca", title="pbmc_3k PCA") #pca
DimPlot_my=function(sce,
                    group.by="seurat_clusters",
                    reduction="umap",

                    isRaster=T,
                    show.legend=F,
                    theme.use=theme_classic,
                    par.geom_point=list(size=0.3),
                    par.geom_text = list(size=4, color="black"),
                    par.circle=list(
                      size=10, #radisu of ring
                      stroke=1.2, #width of ring
                      #color : for ring
                      #fill="white", alpha=0.8,# alpha may affect lines, so set alpha in fill
                      fill="#FFFFFF66" #fill color, last 2 digits is alpha(0.6->#99)
                    ),
                    title="UMAP plot",
                    colSet=NULL)
{
  # get coord and cluster
  umap_data = cbind(
    Embeddings(sce, reduction =reduction)[,1:2],#coord
    data.frame(cluster=as.data.frame(sce@meta.data)[, group.by] ) # cluster
  )
  colnames(umap_data)=c("Dim1", "Dim2", "cluster")
  # head(umap_data)

  # get coord for labels
  dat.plot.label=as.data.frame(t(sapply(
    split(umap_data[,1:2], umap_data[,3]),
    function(x){
      apply(x, 2, median)
    }
  )))
  dat.plot.label$cluster=rownames(dat.plot.label)
  # label circle r
  dat.plot.label$r=1
  dat.plot.label$notion=paste0("c", ifelse( as.numeric(dat.plot.label$cluster)>=10,
                                            dat.plot.label$cluster,
                                            paste0("0", dat.plot.label$cluster)))
  # plot settings
  if(isRaster){
    my.ggPoint <- geom_point_rast
  }else{
    my.ggPoint <- geom_point
  }

  # plot
  g1=ggplot(umap_data, aes(x=Dim1, y=Dim2, color=cluster))+
    do.call(my.ggPoint, c( list(show.legend=show.legend), par.geom_point) )+
    theme.use(base_size = 14)+
    do.call(my.ggPoint, c(list(aes(Dim1, Dim2, color=cluster), data=dat.plot.label,
                               show.legend = F,
                               shape=21), par.circle ) # shape=21 is ring
    )+
    do.call(geom_text, c(list(aes(Dim1, Dim2, label= notion ), data=dat.plot.label,show.legend = F),
                         par.geom_text))+
    labs(title=title)+
    theme( plot.title = element_text(hjust = 0.5) )
  if(!is.null(colSet)){ #set colors
    g1=g1+scale_color_manual(values=colSet)
  }
  return(g1)
}




