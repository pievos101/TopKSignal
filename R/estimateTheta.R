#' Estimation of the underlying signal.
#'
#' The main function for the estimation of the underlying signals is called estimateTheta(). The parameters required are: (1) A rank matrix as described previously, (2) the number of bootstrap samples (it is recommend 500), and (3) a constant for the support variables \(b>0\), recommended is 0.1, (4)  the type of model with four different options fullLinear, fullQuadratic, restrictedLinear and restrictedQuadratic, (5) the type of  bootstrap sampling scheme with two options poisson.bootstrap and classic.bootstrap, and (6) the number of cores that the function can use. Each bootstrap sample is executed on a dedicated CPU core.
#'
#' @param R.input A matrix where the number of rows is equal to the number of the objects and
#' the number of columns equals the number of assessors (rankers).
#' @param b the penalization term. The suggested value is 0.1.
#' @param num.boot the number of boostrap samples created from the input ranked matrix. A positive number is expected.
#' @param solver A string that indicate which solver to use. Two options are available, 'gurobi' and 'nloptr'.
#' We recommend to use gurobi because it is faster. Note, a licence is required. Check the corresponding documentation on how to install gurobi.
#' @param type A string that indicate which model to use: four possibilites are available: 'restrictedQuadratic', 'fullQuadratic', 'restrictedLinear' and 'fullLinear'
#' @param bootstrap.type A string that indicate which bootstrap method to use: two possibilites are available: 'classic.bootstrap' or 'poisson.bootstrap'.
#' @param nCore the number of cores used for computation. Each core is used to solve a bootstrap sample. Default number is detectCores() - 1.
#' @return A list with the estimation information obtained:
#'  \itemize{
#'   \item estimation - A data frame with the signal estimation and the standard error evaluating by the bootstrap for each object
#'   \item estimatedMatrixNoise - The estimated matrix noise 
#'   \item time - the execution time of the procedure
#'   \item allBootstraps - The single bootstrap estimations of all objects for each bootstrap sample
#' }
#' @keywords estimateTheta
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
        print("GUROBI is required <http://www.gurobi.com/>. Is possible to use nloptr a free optimisation tool.")
    }

    if (solver == "nloptr") {
        res <- bootstrap_sqp(R.input, b, num.boot, 
            type = type, bootstrap.type = bootstrap.type, 
            nCore = nCore)
    }

    return(res)
}
