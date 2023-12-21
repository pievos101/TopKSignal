#' violinPlot
#'
#The violin plot displays the bootstrap distribution of the estimated signals along with its means. The deviations from the mean values +/-2 standard errors SE and are shown in the plot.  Analyzing the shape of the distribution and the standard error of the signal of each object, it is possible to evaluate its rank stability with respect to all other objects. The violinPlot function requires (1) the result obtained by the estimation procedure and (2) the 'true' (simulated) signals or ground truth (when available).
#' @param estimation The estimation list from the 'estimateTheta' function
#' @param trueSignal The true signal (if available)
#' @param title The title of the plot  
#' @keywords violintPlot
#' @return A violint plot with the estimated distribution of each object 
#' @examples
#' data(estimatedSignal)
#' violinPlot(estimatedSignal)
#' @export
violinPlot <- function(estimation, trueSignal = NULL, 
    title = NULL) {
    Var2 <- NULL
    SE <- NULL
    value <- NULL
    id <- NULL
    signal.estimate <- NULL
    Type <- NULL
    bootstrapEstimations <- estimation$allBootstraps
    pointEstimations <- estimation$estimation

    allEstimationsMelted <- melt(bootstrapEstimations)
    allEstimationsMelted$Var2 <- as.factor(allEstimationsMelted$Var2)

    if (!is.null(trueSignal)) {
        pointEstimations <- cbind(pointEstimations, 
            signal.true = scaleFun(trueSignal))
    }

    EstimateMelted <- melt(pointEstimations, id.vars = c("id", 
        "SE"))
    colnames(EstimateMelted) <- c("id", "SE", "Type", 
        "value")

    plotEstimation <- ggplot(allEstimationsMelted, 
        aes(x = Var2, y = value)) + geom_violin() + 
        geom_errorbar(data = pointEstimations, aes(x = id, 
            ymin = signal.estimate - 2 * SE, ymax = signal.estimate + 
                2 * SE), inherit.aes = F) + geom_line(data = EstimateMelted, 
        aes(x = id, y = value, group = Type, colour = Type, 
            linetype = Type), size = 1) + xlab("Objects") + 
        ylab("Signal value") + ggtitle(title) + theme(legend.position = "bottom")

    return(plotEstimation)
}
