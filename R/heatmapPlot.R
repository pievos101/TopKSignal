#' Heatmap noise matrix plot
#'
#' The heatmap plot allows us to control for specific error patterns associated with the assessors. The heatmap plot displays information about the noises involved in the estimation process. The rows of the noise matrix are ordered by the estimated ranks of the consensus signal values. The columns are ordered by the column error sums. In the plot, the column with the lowest sum is positioned on the left side and the column with the highest sum is positioned on the right side. Hence, assessors positioned on the left show substantial consensus and thus are more reliable than those positioned to the far right. The heatmap plot is also an exploratory tool for the search for a subset of top-ranked objects (notion of top-$k$ objects ? see the package TopKLists on CRAN for details and functions). Please note, beyond exploratory tasks, the noise matrix can serve as input for various inferential purposes such as testing for assessor group differences. The heatmapPlot function requires the estimation results obtained from the estimateTheta function.
#' @param estimation The bootstrap estimation obtained from the estimateTheta function
#' @param type The type of method used: Two options are available, 'full' or 'reduced'
#' @param title The title of the plot  
#' @keywords heatmapPlot
#' @return A list with:
#' \itemize{
#'   \item plot - A heatmap plot with the noise matrix (ordered values).
#'   \item matrixNoiseOrdered - The matrix noise ordered by the columns. The objects are ordered by the estimated value.
#'   \item estimateThetaOrdered - The theta vector ordered by their importance (from the highest value to the lowest).
#' }
#' @examples
#' data(estimatedSignal)
#' heatmapPlot(estimatedSignal)
#' @export
heatmapPlot <- function(estimation, type = "full", 
    title = "") {
    Var2 <- NULL
    Var1 <- NULL
    value <- NULL
    ord2 <- order(colSums(estimation$estimatedMatrixNoise), 
        decreasing = F)
    ord <- order(estimation$estimation$signal.estimate, 
        decreasing = T)

    if (type == "full") {
        matNoise <- estimation$estimatedMatrixNoise[ord, 
            ord2]
    }
    if (type == "reduced") {
        matNoise <- estimation$estimatedMatrixNoise[, 
            ord2]
    }

    matrixNoiseMelted <- melt(matNoise)

    toRet <- ggplot(data = matrixNoiseMelted, aes(x = Var2, 
        y = Var1, fill = value)) + geom_tile() + ggtitle(title) + 
        xlab("Assessor") + ylab("Position") + scale_fill_gradient(low = "white", 
        high = "black") + theme(legend.position = "bottom") + 
        scale_y_reverse(breaks = function(x) unique(floor(pretty(seq(1, 
            (max(x) + 1) * 1.1)))))


    temp <- estimation$estimation[ord, ]
    rownames(temp) <- 1:nrow(estimation$estimation)

    return(list(plot = toRet, matrixNoiseOrdered = matNoise, 
        estimateThetaOrdered = temp))
}
