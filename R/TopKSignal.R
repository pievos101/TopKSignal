#' TopKSignal: A convex optimization tool for signal reconstruction from multiple ranked lists.
#' 
#' A convex optimisation tool to estimate the underlying signal (latent variable) of the global rank order using quadratic or linear convex optimisation. 
#' The goal of consensus across the assessors is achieved by an indirect inference using a sets of order constraints.
#' TopKSignal implements a set of different functions. They permit to construct artificial ranked lists, to build sets of constraints from a rank matrix, to perform bootstrap estimation (standard and Poisson bootstrap), to run convex optimisation (with linear and quadratic objective functions), and to obtain numerical and graphical output.
#' For signal estimation two different methods are available: the restricted and the full method with two different penalization methods: quadratic and linear.
#' Two different boostrap sample schemes are implemented: the classic bootstrap and the Poisson bootstrap. 
#'
#' @section estimateTheta function:
#' The main function for the estimation of the underlying signals is called estimateTheta(). The parameters required are: (1) A rank matrix, (2) the number of bootstrap samples (500 is recommended), and (3) a constant for the support variables \(b>0\), default is 0.1, (4)  the type of model with four different options fullLinear, fullQuadratic, restrictedLinear and restrictedQuadratic, (5) the type of  bootstrap sampling scheme: poisson.bootstrap and classic.bootstrap, and (6) the number of cores for parallel computation. Each bootstrap sample is executed on a dedicated CPU core.
#'
#' @section generate.rank.matrix function:
#' The generate.rank.matrix() function requires the user to specify the number of objects (items), called p, and the number of assessors (rankers), called n. The function simulates full ranked lists (i.e. no missing assignments) without ties.
#'
#' @section violintPlot function:
#' The violint plot displays the bootstrap distribution of the estimated signals along with the mean estimates. The deviations from the mean values are equal to 2+-SE and are displayed in the violin plot.  With the analysis of the distribution and the standard error of the signal of each object, it is possible to verify its rank stability compared to all other objects. The violinPlot function requires (1) the result obtained by the estimation procedure and (2) the 'true' (simulated) signals or ground truth (when available).
#'
#' @section heatmapPlot function:
#' The heatmap plot allows us to check whether specific patterns are present in the data. The heatmapPlot contains the information about the involved noises in the estimation process. The rows of the noise matrix are ordered by the estimated rank of the signal value. The columns are ordered by the column sum. The column with the lowest sum is positioned on the left side and the column with the highest sum is positioned on the right side. Hence, the assessors positioned on the left should be more reliable than those on the right. When some assessors of low ranking ability are suspected the heatmap is very helpful. This plot is also useful when an informative set of top-ranked objects is likely or in general when patterns inside the data might be present. The heatmapPlot function requires the estimation results obtained from the estimateTheta function.
#'  
#' @section elbowPlot function:
#' The Elbow plot permits to identify a subset of top-$k$ objects. Objects are ordered according to their rank position (x-axis). On the y-axis we have the corresponding estimated signal value of each object. The idea of the elbow plot is to scan for 'jumps' - neighbouring signal estimates which are visually much distant - in an exploratory manner. The elbowPlot function requires the estimation results from the estimateTheta function.
#' 
#' @docType package
#' @name TopKSignal
NULL
