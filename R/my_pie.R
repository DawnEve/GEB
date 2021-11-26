#' Using a dataframe specific column to plot a pie figure
#'
#' @param data a dataframe or tibble
#' @param ncol specific column, factor
#'
#' @return a figure as ggplot2 object
#' @export
#'
#' @examples
#' my_pie(mtcars,ncol=2)
#' my_pie(diamonds,ncol=3)
#' my_pie(diamonds,ncol=3)+guides(fill="none")
#'
my_pie <- function(data, ncol){
  plotdat <- as.data.frame(table(data[,ncol])) %>% dplyr::arrange(-Freq)
  plotdat$Label <- paste(plotdat$Var1, paste0("(",round(((plotdat$Freq/sum(plotdat$Freq))*100),2),"%)"))

  p <- ggplot(plotdat, aes (x="", y = Freq, fill = factor(Var1))) +
    geom_col(position = 'stack', width = 1) +
    geom_text_repel(aes(label = Label, x = 1.3),size=5,
                    position = position_stack(vjust = 0.5)) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(fill = "Category",x = NULL,y = NULL) +
    coord_polar("y")
  return(p)
}
