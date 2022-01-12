#' My enhanced DotPlot
#' version:0.1
#'
#' @param object Seurat obj
#' @param assay which assay
#' @param features gene list
#' @param cols colors for expression
#' @param col.min min?
#' @param col.max max?
#' @param dot.min dot min size
#' @param dot.scale dot scale method
#' @param idents which idents to show
#' @param group.by group by
#' @param split.by split by
#' @param cluster.idents is cluster idents?
#' @param scale scale?
#' @param scale.by scale by?
#' @param scale.min min?
#' @param scale.max max?
#' @param base_size font size list
#' @param rel_widths relative widths
#' @param color_list color list for each cluster
#'
#' @return
#' @export
#'
#' @examples
#' DotPlot_my(sce, features = c("CD3D", "CD79A"))
#' DotPlot_my(sce, features = c("CD3D", "CD79A"), cols = c("white", "red2"), rel_widths = c(1,4))
#' colors=c("#96C3D8", "#5F9BBE", "#F5B375", "#C0937E", "#67A59B", "#A5D38F", "#4A9D47", "#F19294", "#E45D61", "#3377A9",
#' "#BDA7CB", "#684797", "#8D75AF", "#CD9C9B", "#D62E2D", "#DA8F6F", "#F47D2F")
#' DotPlot_my(pbmc, features = split(top11$gene, top11$cluster), cols = c("white", "red3"), color_list=colors, rel_widths=c(1,45))
#' DotPlot_my(scObj_colon, features = common_markers, color_list = colors,
#'           rel_widths = c(1,35), base_size = list(legend=10, x=12))
DotPlot_my=function (object, assay = NULL, features,
                     cols = colorRampPalette(colors=c("lightyellow", "orangered", "red4"),
                                             bias=0.5,
                                             interpolate = c("spline"))(100),
                     col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                     idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE,
                     scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA,

                     circle_size=6, #左侧圆圈大小
                     base_size=list(x=14, y=14, top=14, legend=14), #文字字号
                     rel_widths=c(1, 8), #左右比例
                     color_list=NULL) #颜色列表
{
  library(cowplot)
  # 返回高于某个阈值的细胞百分比
  PercentAbove=function (x, threshold)
  {
    return(length(x = x[x > threshold])/length(x = x))
  }


  # default values
  if(!("x" %in% names(base_size))){ base_size$x=14 }
  if(!("y" %in% names(base_size))){ base_size$y=14 }
  if(!("top" %in% names(base_size))){ base_size$top=14 }
  if(!("legend" %in% names(base_size))){ base_size$legend=14 }


  # 哪个 assay，默认是 RNA
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay

  # 是否分割，必须定义 split.by 且 cols 是画板中的值
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))

  # scale.by 是啥？按面积或者按半径？
  #switch(3, "a", "B", "C") #C
  #switch("size", size = scale_size, radius = scale_radius)
  scale.func <- switch(EXPR = scale.by,
                       size = scale_size,
                       radius = scale_radius,
                       stop("'scale.by' must be either 'size' or 'radius'"))


  feature.groups <- NULL
  # 如果是list，则解开
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features),
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
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

  ################
  # 获取表达数据，
  ################
  # 获取 给定 cells 的 features 值(基因表达矩阵，默认是 slot="data")也就是 normalize 后的数据
  data.features <- FetchData(object = object, vars = features, cells = cells)

  # 生成 id 列
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE] # 如有 无 group.by
  }
  else {
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


  ###############################
  # 获取平均表达量和表达百分比，按类 id 列
  ###############################
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


  # 如果对 cluster 进行聚类（可选）
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


  ##############
  # scale 各类(id列)的 mean 值
  ##############
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
  # browser() #2----- //在这里调试点
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled)) #一列一个基因，变成了按行拼接成的向量


  # 如果有分开的颜色
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }

  # browser(text = "3") #3----- //在这里调试点

  data.plot$avg.exp.scaled <- avg.exp.scaled # 记录 scale 后的mean
  data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features) #变因子

  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100

  #
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id),
                         FUN = gsub, FUN.VALUE = character(length = 1L),
                         pattern = paste0("^((", paste(sort(x = levels(x = object),
                                                            decreasing = TRUE), collapse = "|"), ")_)"),
                         replacement = "", USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")

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


  #################
  ## ggplot2 画主图 plot 对象
  #################
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) +
    geom_point(mapping = aes_string(size = "pct.exp", fill = color.by),
               shape=21, stroke=0.5) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    # scale_x_continuous(expand=c(0,0))+
    theme(
      text=element_text(size=base_size$legend), #文字字号
      plot.margin =margin(l=0), #no margin on the left

      # 每个画板
      panel.background = element_blank(), #去背景灰
      panel.border = element_rect(color="black", size=1, fill="#00112200"), #要边框,填充为透明
      panel.grid = element_line(colour = "grey92", size = rel(0.5)),# 方格线
      panel.spacing.x = unit(1, "mm"), #各个画板的间距

      # 分面顶标题
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(size=base_size$top),
      # 图例背景
      legend.key = element_rect(fill = "#00112200"),

      # 坐标轴不要刻度
      axis.ticks = element_blank(),
      # 坐标轴不要label
      axis.title = element_blank(), # axis.title.x = element_blank(), axis.title.y = element_blank(),

      #x坐标文字 gene 旋转60度
      axis.text.x.bottom = element_text(hjust = 1, vjust = 1, size=base_size$x,  angle = 60 ),

      # 不要y坐标文字
      axis.text.y = element_blank()
    )+
    guides(size = guide_legend(title = "Percent Expressed"))
  #labs(x = "Features", y = ifelse(test = is.null(x = split.by), yes = "Identity", no = "Split Identity")) +



  # 如果gene有分组，就分面
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups,
                              switch = "y",
                              space = "free_x",
                              scales = "free_x")+
      theme(
        # 每个画板
        panel.spacing = unit(x = 1, units = "lines"),

        # 不要顶部分面标签
        #strip.text = element_blank(),
        #strip.background = element_blank(),
      )
  }

  # 如果有分组颜色
  if (split.colors) {
    plot <- plot + scale_fill_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_fill_distiller(palette = cols)
  }
  else if (length(x = cols) == 2) {
    plot <- plot + scale_fill_gradient(low = cols[1], high = cols[2])
  }
  else if (length(x = cols) == 3) {
    plot <- plot + scale_fill_gradient2(low = cols[1], mid=cols[2], high = cols[3])
  }
  else if (length(x = cols) >3) {
    plot <- plot + scale_fill_gradientn( colors  = cols )
  }

  # 如果没分组颜色，添加填充图例
  if (!split.colors) {
    plot <- plot + guides(fill = guide_colorbar(title = "Average Expression"))
  }

  ############################
  # left: color circle list
  ############################
  len=length(unique(data.plot$id));
  if(is.null(color_list)){
    color_list=rainbow(len)
  }
  g1_left=ggplot(data.frame(x=0, y= unique(data.plot$id) ), aes(x,y, color=factor(y) ))+
    geom_point(size=circle_size, show.legend = F)+
    scale_color_manual(values=color_list)+theme_classic()+
    scale_x_continuous(expand=c(0,0))+
    theme(
      plot.margin =margin(r=0), #no margin on the right
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text.y=element_text(size=base_size$y)
    );


  #library(cowplot)
  # https://wilkelab.org/cowplot/articles/aligning_plots.html
  # we can align both the bottom and the top axis (axis = "bt").
  plot2=plot_grid(g1_left, plot, align ="h", axis="bt", rel_widths = rel_widths)
  return(plot2)
}
