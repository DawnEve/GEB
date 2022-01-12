# Part II other functions==================




#' Title: Make Sankey plot data from a table(paraLeft, paraRight)
#' ref: https://observablehq.com/ at d3/sankey
#' version: 2.0 define flow color as source//done
#' v2.1 supply 3 color method: static, source, target//done
#'
#' @param tbl1 a table produced by table()
#' @param omit omit which rows/paraLeft
#' @param prefix_left add a prefix to left nodes
#' @param prefix_right add a prefix to right nodes
#' @param colorMethod "static", "source", "target"
#' @param omitR omitR rownames
#' @param omitC omitC column names
#'
#' @return
#' @export
#'
#' @examples
#' getSankeyData( table(mtcars$gear, mtcars$carb))
getSankeyData=function(tbl1,
                       prefix_left="", prefix_right="",
                       colorMethod="source",
                       omitR=NULL, omitC=NULL, rm.na=T){
  preL=prefix_left
  preR=prefix_right
  #
  if(!is.null(omitR)){
    tbl1=tbl1[setdiff( rownames(tbl1), omitR), ]
  }
  if(!is.null(omitC)){
    tbl1=tbl1[ , setdiff( colnames(tbl1), omitC)]
  }
  # default, remove all 0 rows and columns
  tbl1=tbl1[rowSums(tbl1)>0, ]
  tbl1=tbl1[ , colSums(tbl1)>0]
  # table to data.frame
  df2=as.data.frame(tbl1)
  df2=df2[which(df2$Freq!=0),]
  #
  df2$Var1=paste0(preL, df2$Var1)
  df2$Var2=paste0(preR, df2$Var2)
  #
  nodes=data.frame(node=unique( c( paste0(preL, rownames(tbl1)),
                                   paste0(preR, colnames(tbl1))) ))
  nodes$index=0:(nrow(nodes)-1) #0-based index, as in JS
  df3=data.frame(
    source=nodes$index[match(df2$Var1, nodes$node)],
    target=nodes$index[match(df2$Var2, nodes$node)],
    value=df2$Freq
  )
  # use color as node in source
  if("source"==colorMethod){
    df3$group= nodes$node[(1+df3$source)]
  }else if("target"==colorMethod){
    df3$group= nodes$node[(1+df3$target)]
  }else if("static"==colorMethod){
    df3$group="grey";
  }

  return(list(link=df3, node=nodes))
}


#' Sankey plot from sankeyDt, produced by getSankeyData()
#' version: 2.0 define flow color as source//done
#'
#' @param sankeyDt dataset produced by getSankeyData
#' @param colors a color list for sankeyDt$node$node
#' depend getSankeyData()
#'
#' @return
#' @export
#'
#' @examples
#' getSankeyData( table(mtcars$gear, mtcars$carb))
#' sankeyDt2=getSankeyData(table(mtcars$gear, mtcars$carb),"Gear_", "Carb_")
#' sankeyPlot(sankeyDt2 )
#' sankeyPlot( getSankeyData( table(mtcars$gear, mtcars$carb), "Gear_", "Carb_", colorMethod="source") )
#' sankeyPlot( getSankeyData( table(mtcars$gear, mtcars$carb), "Gear_", "Carb_", colorMethod="target") )
#' sankeyPlot( getSankeyData( table(mtcars$gear, mtcars$carb), "Gear_", "Carb_", colorMethod="static") )
#'
#' library(RColorBrewer);
#' mypalette_1=brewer.pal(8,"Set1")
#' sankeyPlot(sankeyDt2, colors=mypalette_1 )
sankeyPlot=function(sankeyDt, colors=NULL){
  library(RColorBrewer)
  library(networkD3)
  if(is.null(colors)){
    #define the color of each nodes(from/source, to/target)
    colors = brewer.pal(8,"Set2")
    #c("#FF0000", "#00FF00", "#0000FF", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#FF00004D", "#0000FF4D")
  }
  #R array to json
  colorsJson=  paste0('d3.scaleOrdinal().domain(', jsonlite::toJSON(sankeyDt$node$node), ').range(',  jsonlite::toJSON(colors) , ')');
  #colorsJson <- 'd3.scaleOrdinal().domain(["pH", "C", "N", "P", "K", "1", "-1"]).range(["#377EB8", "#4DAF4A", "#984EA3"])'
  params=list(Links = sankeyDt$link, Nodes = sankeyDt$node,
              Source = "source", Target = "target", Value = "value",
              NodeID = "node",
              colourScale = colorsJson, #JS("d3.scaleOrdinal(d3.schemeCategory20);"),
              fontSize = 12, nodeWidth = 30)
  # if only one type of links, then it must be static grey!
  if( length(unique(sankeyDt$link$group))>1 ){
    params$LinkGroup = "group" #color of flow
  }
  do.call(sankeyNetwork, params)
}
# test
#temp=table(scObj_known2@meta.data$oCluster, scObj_known2@meta.data$seurat_clusters)
# rownames(temp)[rowSums(temp)<10] #check to remove rows with very small numbers
#sankeyPlot( getSankeyData(temp, "old_") )
#sankeyPlot( getSankeyData(temp, "old_", omitR=c("10", "22", "24", "25", "28", "29") ) )
#sankeyPlot( getSankeyData(temp, "old_", omitR=c("10", "22", "24", "25", "28", "29") ), colors= brewer.pal(8,"Set2") )

