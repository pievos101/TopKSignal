#' elbowPlot
#'
#' The elbow plot permits the identification of subsets of objects, e.g. top-$k$ or bottom-$q$ objects. On the x-axis all objects are ordered according to their rank positions. On the y-axis the corresponding estimated signal values are displayed. The idea of the elbow plot is to scan for 'jumps' in the sequence of ordered objects ? i.e. find signal estimates next to each other that are visually much distant - in an exploratory manner. The elbowPlot function requires the estimation results from the estimateTheta function.
#'
#' @param estimation Results from the estimateTheta() function
#' @param title A title for the plot  
#' @keywords elbowtPlot
#' @return A elbow plot 
#' @examples
#' data(estimatedSignal)
#' elbowPlot(estimatedSignal)
#' @export
elbowPlot <- function(estimation, title = "") {
    estimation <- estimation$estimation
    ord <- order(estimation$signal.estimate, decreasing = T)
    signalOrdered <- estimation$signal.estimate[ord]
    position <- 1:length(ord)

    elbowPlotRes <- ggplot(data = data.frame(), aes(x = position, 
        y = signalOrdered)) + geom_point() + geom_text(aes(label = estimation$id[ord], 
        hjust = -0.1, vjust = -0.5)) + ggtitle(title) + 
        xlab("Rank position") + ylab("Estimated value") + 
        scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, 
            (max(x) + 1) * 1.1)))))

    return(elbowPlotRes)
}
