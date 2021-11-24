#' Heatmap noise matrix plot
#'
#' The heatmap plot allows us to check whether specific patterns are present in the data. The heatmapPlot contains the information about the involved noises in the estimation process. The rows of the noise matrix are ordered by the estimated rank of the signal value. The columns are ordered by the column sum. The column with the lowest sum is positioned on the left side and the column with the highest sum is positioned on the right side. Hence, the assessors positioned on the left should be more reliable than those on the right. When some assessors of low ranking ability are suspected the heatmap is very helpful. This plot is also useful when an informative set of top-ranked objects is likely or in general when patterns inside the data might be present. The heatmapPlot function requires the estimation results obtained from the estimateTheta function. 
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
