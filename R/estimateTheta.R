#' Estimation of the underlying signal.
#'
#' The main function for the estimation of the signals informing the ranks is called estimateTheta(). The required parameters are: (1) a rank matrix, (2) the number of bootstrap samples (500 is recommended), (3) a constant for the support variables \(b>0\), default is 0.1, (4)  the type of optimization technique: fullLinear, fullQuadratic, restrictedLinear, and restrictedQuadratic (the latter two recommended), (5) the type of  bootstrap sampling scheme: classic.bootstrap and poisson.bootstrap (recommended), and (6) the number of cores for parallel computation. Each bootstrap sample is executed on a dedicated CPU core.
#' 
#' @param R.input A matrix where the rows represent the objects and
#' the columns the assessors (rankers).
#' @param b The penalization term. The suggested value is 0.1.
#' @param num.boot The number of boostrap samples created from the input ranked matrix. A positive number is expected.
#' @param solver A string that indicates which solver to use. Two options are available, 'gurobi' and 'nloptr'.
#' We recommend to use gurobi for faster computation. Note, a licence is required. Check the corresponding documentation on how to install gurobi.
#' @param type A string that indicates which model to use: four approaches are available: 'restrictedQuadratic', 'fullQuadratic', 'restrictedLinear' and 'fullLinear'.
#' @param bootstrap.type  A string that indicates which bootstrap method to use: 'classic.bootstrap' or 'poisson.bootstrap'.
#' @param nCore The number of cores used for computation. Each core is used to calculate the signals from a bootstrap sample. Default number is detectCores() - 1.
#' @return A list with the estimation information obtained:
#'  \itemize{
#'   \item estimation - A data frame with the signal estimation and the standard error computed by the bootstrap for each object
#'   \item estimatedMatrixNoise - The estimated matrix noise 
#'   \item time - The execution time of the procedure
#'   \item allBootstraps - The signal estimates from all bootstrap iterations
#' }
#' @keywords estimateTheta
#' @examples
#' library(TopKSignal)
#' set.seed(1421)
#' p = 8
#' n = 10
#' input <- generate.rank.matrix(p, n)
#' rownames(input$R.input) <- c("a","b","c","d","e","f","g","h")
#' # For the following code Gurobi needs to be installed
#' \dontrun{
#' estimatedSignal <- estimateTheta(R.input = input$R.input, num.boot = 50, b = 0.1, 
#' solver = "gurobi", type = "restrictedQuadratic", bootstrap.type = "poisson.bootstrap",nCore = 1)   
#' }
#' data(estimatedSignal)
#' estimatedSignal
#' @export
estimateTheta <- function(R.input, b, num.boot, solver, 
    type, bootstrap.type, nCore = ((detectCores() - 
        1))) {
    if (requireNamespace("gurobi", quietly = TRUE)) {
        requireNamespace("gurobi")
        if (solver == "gurobi") {
            res <- bootstrap_gurobi(R.input, b, num.boot, 
                type = type, bootstrap.type = bootstrap.type, 
                nCore = nCore)
        }
    } else {
        stop("GUROBI is required <http://www.gurobi.com/>. Is possible to use nloptr a free optimisation tool.")
    }

    if (solver == "nloptr") {
        res <- bootstrap_sqp(R.input, b, num.boot, 
            type = type, bootstrap.type = bootstrap.type, 
            nCore = nCore)
    }

    return(res)
}
