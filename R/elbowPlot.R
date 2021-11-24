#' elbowPlot
#'
#' The Elbow plot permits to identify a subset of top-$k$ objects. Objects are ordered according to their rank position shown on the x-axis. On the y-axis we have the corresponding estimated signal value of each object. The idea of the elbow plot is to scan for 'jumps' - neighbouring signal estimates which are visually much distant - in an exploratory manner. The elbowPlot function requires the estimation results from the estimateTheta function.
#' @param estimation Results from the estimateTheta() function
#' @param title A title for the plot  
#' @keywords elbowtPlot
#' @return A elbow plot 
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
