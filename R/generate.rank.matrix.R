#' generate.rank.matrix
#'
#' The generate.rank.matrix() function requires the user to specify the number of objects (items), called p, and the number of assessors, called n. The function simulates full ranked lists (i.e. no missing assignments) without ties.
#'
#' @param p The number of objects.
#' @param n The number of assessors.
#' @param percentageMissing The percentage of the missing values. Note, missing data should be resolved by the rank() function before calling estimateTheta().
#' @return A list with simulated data
#' \itemize{
#'   \item R.input - The rank matrix 
#'   \item thea.true - The true underlying signals from the assessments 
#'   \item sigmas - The standard error of the noise added for each assessor
#'   \item matrixNoise - The noise added to the true signals in order to get the final rank matrix
#' }
#' @keywords generate.rank.matrix
#' @examples
#' p = 8
#' n = 10
#' input <- generate.rank.matrix(p, n)
#' rownames(input$R.input) <- c("a","b","c","d","e","f","g","h")
#' @export
generate.rank.matrix <- function(p, n, percentageMissing = 0) {
    # p - number of objects n- number of
    # assessors/rankers
    sigmas = runif(n, min = 0.4, max = 0.6)
    X.input = matrix(nrow = p, ncol = n)  # the matrix of observed attributes (X=theta+Z)
    X.noise <- matrix(nrow = p, ncol = n)
    theta.true <- abs(rnorm(p, 0, 1))  # true signal
    for (i in 1:n) {
        noise <- abs(rnorm(p, 0, sigmas[i]))
        X.noise[, i] <- noise
        X.input[, i] <- theta.true + noise
    }
    nElements <- p * n
    nMiss <- round(nElements * percentageMissing/100)

    if (percentageMissing > 0) {
        repeat {
            pR <- sample(1:p, nMiss, replace = T)
            nR <- sample(1:n, nMiss, replace = T)

            if (nrow(unique(cbind(pR, nR))) == nMiss) 
                break
        }

        for (i in 1:nMiss) {
            X.input[pR[i], nR[i]] <- NA
        }
    }

    R.input = matrix(nrow = p, ncol = n)  # the observed rankings

    R.input = apply(X.input, 2, function(x) rank(-x, 
        na.last = "keep", ties.method = "random"))
    rownames(R.input) <- 1:p
    colnames(R.input) <- paste0("R", 1:n)

    return(list(R.input = R.input, theta.true = theta.true, 
        sigmas = sigmas, matrixNoise = X.noise))
}
