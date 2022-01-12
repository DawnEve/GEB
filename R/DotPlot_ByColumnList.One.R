#' facet_grid DotPlot by column, only accept named marker gene list.
#' v0.3
#' 模仿张泽民实验室的图(https://pubmed.ncbi.nlm.nih.gov/34914499/  Fig.2A)，竖排，分面方框
#'
#' @param object Seurat obj
#' @param features a list of named marker genes
#' @param scale scale?
#' @param scale.by scale method
#' @param idents idents to use
#' @param cluster.idents whether cluster idents
#' @param base_size font size
#' @param cols colors of gene expression
#' @param direction color RColorBrewer direction
#' @param col.min min?
#' @param col.max max?
#' @param dot.min dot size min
#' @param dot.scale dot max size
#' @param scale.min min?
#' @param scale.max max?
#' @param group.by group by
#' @param split.by split by
#'
#' @return
#' @export
#'
#' @examples
#' DotPlot_ByColumnList(scObj_colon, features = common_markers)
DotPlot_ByColumnList=function(object, features, scale = TRUE, scale.by = "radius", idents=NULL,
                              cluster.idents = FALSE,
                              base_size=list(),
                              cols=NULL,
                              direction=1,
                              col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                              scale.min = NA, scale.max = NA,
                              group.by=NULL, split.by=NULL){
  if(!is.list(features)){
    stop("Error: features must be a named list")
  }

  # 返回高于某个阈值的细胞百分比
  PercentAbove=function (x, threshold)
  {
    return(length(x = x[x > threshold])/length(x = x))
  }

  # default values
  if(!("x" %in% names(base_size))){ base_size$x=10 }
  if(!("y" %in% names(base_size))){ base_size$y=10 }
  if(!("right" %in% names(base_size))){ base_size$right=10 }
  if(!("legend" %in% names(base_size))){ base_size$legend=6 }



  scale.func <- switch(EXPR = scale.by,
                       size = scale_size,
                       radius = scale_radius,
                       stop("'scale.by' must be either 'size' or 'radius'"))

  feature.groups <- NULL
  # 如果是list，则解开
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features),
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x],
                                                     each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE,
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }

  # 获取 idents 的细胞
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))


  ##==--  获取表达数据
  # 获取 给定 cells 的 features 值(基因表达矩阵，默认是 slot="data")也就是 normalize 后的数据
  data.features <- FetchData(object = object, vars = features, cells = cells)

  # 生成 id 列
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE] # 如有 无 group.by
  }else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE] #如果有 group.by
  }


  # 如果id列不是factor，变 factor
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  # id 列的 unique 值
  id.levels <- levels(x = data.features$id)
  # id 列 变回 vector
  data.features$id <- as.vector(x = data.features$id)


  # 如果有分割~分面
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits,
                              sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }


  ##==--  获取平均表达量和表达百分比，按类 id 列
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), #去掉最后一列 id
                              drop = FALSE] # drop 啥意思？

    # 求平均，按列。求的是取对数前-1，也就是 counts的平均。
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x))) #expm1 是啥？ expm1(x) computes exp(x) - 1 accurately also for |x| << 1.
    })

    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })

  names(x = data.plot) <- unique(x = data.features$id)



  ##==-- 如果对 cluster 进行聚类（可选）
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, FUN = unlist)) #解开list
    mat <- scale(x = mat) #按列scale
    id.levels <- id.levels[hclust(d = dist(x = mat))$order] #按聚类结果对cluster排序
  }

  # data.plot 是个 list，name就是各个id
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x

    # browser() #1----- //在这里调试点
    #data.use
    #    avg.exp pct.exp features.plot          id
    #    CD3D  9.48778156    0.84          CD3D Naive CD4 T
    #    CD79A 0.08742276    0.02         CD79A Naive CD4 T
    return(data.use)
  })
  # 合并list为df
  data.plot <- do.call(what = "rbind", args = data.plot)

  # id列转为 factor
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }


  ##==--  scale 各类(id列)的 mean 值
  # 如果只有一类，则不进行 scale 。怎么 scale 呢？
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE)
  }

  # 按 基因 取出其各个类的mean，并 scale。也就是一对基因必有一个avg.exp 的最值，该假设比较合理。
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == x, "avg.exp"] #取 平均表达值的一列
                             if (scale) { #如果大于1个类，则scale
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use) #log1p(x) computes log(1+x) accurately also for |x| << 1.
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled)) #一列一个基因，变成了按行拼接成的向量

  data.plot$avg.exp.scaled <- avg.exp.scaled # 记录 scale 后的mean
  data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features) #变因子

  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100


  # 最小
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  # 最大
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }

  #
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot],
                                       levels = unique(x = feature.groups))
  }

  # plot using data.frame data.plot
  plot_g1=ggplot(data = data.plot, mapping = aes_string(y = "features.plot", x = "id")) +
    geom_point(mapping = aes_string(size = "pct.exp", color ="avg.exp.scaled")) +
    scale_radius(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    # scale_x_continuous(expand=c(0,0))+
    facet_grid(facets = feature.groups~.,
               #switch = "y", #y的label由默认的右侧转为左侧
               space = "free_y", #自由分配x轴的空间
               scales = "free_y")+ #x坐标范围自由
    theme_grey(base_size = 12)+
    theme(
      text=element_text(size=base_size$legend),
      # 每个画板
      panel.background = element_blank(), #去背景灰
      panel.border = element_rect(color="black", size=1, fill="#00112200"), #要边框,填充为透明
      panel.grid = element_line(colour = "grey92", size = rel(0.5)),# 方格线
      panel.spacing.y = unit(1, "mm"), #各个画板的间距

      # 分面标题
      strip.background = element_rect(fill = "#00112200", colour = "black"),

      strip.placement = "outside", #分面标签和主图分离
      strip.text.y = element_text(margin = margin(l=0.3, unit="cm"), angle = 90), # 分面标题框内边距
      strip.text = element_text(size=base_size$right, face="bold"),

      # 图例背景
      legend.key = element_rect(fill = "#00112200"),

      #坐标轴
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), #xlab, ylab
      axis.text.x.bottom = element_text(hjust = 1, vjust = 1, size=base_size$x,  angle = 60 ), #text.x
      axis.text.y=element_text(size=base_size$y, face="bold", margin = margin(r=-1, unit="pt")), #text.y
    )+
    guides(size = guide_legend(title = "Percent Expressed"),
           color = guide_colorbar(title = "Average Expression"))

  # 设置渐变颜色
  if(is.null(cols)){
    plot_g1 <- plot_g1 + scale_color_distiller(palette = "YlOrRd", direction =direction)
  }
  else if (length(x = cols) == 1) {
    plot_g1 <- plot_g1 + scale_color_distiller(palette = cols, direction =direction)
  }
  else if (length(x = cols) == 2) {
    plot_g1 <- plot_g1 + scale_color_gradient(low = cols[1], high = cols[2])
  }
  else if (length(x = cols) == 3) {
    plot_g1 <- plot_g1 + scale_color_gradient2(low = cols[1], mid=cols[2], high = cols[3])
  }
  else if (length(x = cols) >3) {
    plot_g1 <- plot_g1 + scale_color_gradientn( colors  = cols )
  }


  library(grid)
  lg = linesGrob(x=unit(c(0,0)+0.2,"npc"), y=unit(c(0,1),"npc"),
                 gp=gpar(col="black", lwd=3))

  #grid.newpage();
  # grid.draw(lg) #预览
  q = ggplotGrob(plot_g1) #ggplot2 to grid Grob
  # q$layout$name #先预览

  # 替换
  for (k in grep("strip-r",q$layout$name)) {
    q$grobs[[k]]$grobs[[1]]$children[[1]] = lg
  }
  # 画图
  # grid.draw(q)
  # message("Draw the plot: > grid::grid.draw(obj)")
  # gave an class to draw when print
  class(q) <- c("GEB_lite", class(q))
  return(q)
}

#' 回车直接画 grid 图形
#'
#' @param x an Grob object with class GEB_lite
#' @param newpage whether newpage?
#'
#' @return
#' @export
#'
#' @examples
#' print(sth)
print.GEB_lite <- function(x, newpage = TRUE) {
  if (newpage) {
    grid::grid.newpage()
  }
  grid::grid.draw(x)
  invisible(x)
}


