# return the indices of the noise with row
# i-esima and column j-esima
getNoiseIJ <- function(i, j, p) {
    return(j * p + i)
}

# get theta
getTheta <- function(res, p) {
    return(res[1:p])
}

# get the matrix noise estimated
getNoise <- function(res, p) {
    matrixNoise <- matrix(res[(p + 1):length(res)], 
        nrow = p, byrow = F)
    return(matrixNoise)
}

st.error.estim <- function(x) {
    sqrt(sum((x - mean(x))^2/(length(x) - 1)))
}  # standard error estimation


generate.bootstrap.samples <- function(R, num.boot) {
    # R - the input matrix num.boot - number of
    # bootstrap matrices
    n <- ncol(R)
    p <- nrow(R)
    boots <- list()  # list of bootstrap matrices
    boots[[1]] = R

    indices <- list()
    indices[[1]] <- 1:ncol(R)

    for (b in 2:(num.boot + 1)) {
        ind = sample.int(n, n, replace = TRUE)
        boots[[b]] = R[, ind]
        indices[[b]] <- ind
    }
    return(list(boots = boots, indices = indices))
}


generate.poisson.bootstrap.samples <- function(R, num.boot) {
    n <- ncol(R)
    p <- nrow(R)

    m <- floor((1 - exp(1)^(-1)) * n)

    boots <- list()
    weights <- list()
    indices <- list()
    for (i in 1:num.boot) {
        ind <- sample(x = n, size = m, replace = F)
        weightsTemp <- qpois(runif(m, dpois(0, 1), 
            1), 1)

        boots[[i]] <- R[, ind]
        weights[[i]] <- weightsTemp
        indices[[i]] <- ind
    }
    return(list(boots = boots, weights = weights, indices = indices))
}


generate.poisson.online.bootstrap.samples <- function(R, 
    num.boot) {
    n <- ncol(R)
    p <- nrow(R)

    boots <- list()
    weights <- list()
    indices <- list()
    for (i in 1:num.boot) {
        weightsTemp <- rpois(n, 1)
        ind <- which(weightsTemp != 0)
        weightsTemp <- weightsTemp[which(weightsTemp != 
            0)]

        boots[[i]] <- R[, ind]
        weights[[i]] <- weightsTemp
        indices[[i]] <- ind
    }
    return(list(boots = boots, weights = weights, indices = indices))
}

# get results from bootstrap
gather.results <- function(res) {
    num.boot = length(res) - 1
    p <- length(res[[1]])

    main.nrm.estimates = matrix(nrow = p, ncol = (num.boot + 
        1))

    for (b in 1:(num.boot + 1)) {
        main.nrm.estimates[, b] = res[[b]]
    }

    est.plus.SE = data.frame(signal.estimate = apply(main.nrm.estimates, 
        1, mean), SE = apply(main.nrm.estimates, 1, 
        st.error.estim))
    return(est.plus.SE)
}


gather.matrix.results <- function(matrix.noise.list, 
    bootstrap.indices, bootstrap.weights, R.input) {
    # matrix.noise.list <- matrix_noie_all
    # bootstrap.indices <- indices
    # bootstrap.weights <- weights
    if (!is.null(bootstrap.weights)) {
        for (i in 1:length(matrix.noise.list)) {
            for (j in 1:ncol(matrix.noise.list[[i]])) {
                matrix.noise.list[[i]][, j] <- matrix.noise.list[[i]][, 
                  j] * bootstrap.weights[[i]][j]
            }
        }

    }

    matInfo <- matrix(0, ncol = ncol(R.input), nrow = nrow(matrix.noise.list[[1]]))
    for (i in 1:length(matrix.noise.list)) {
        for (j in 1:ncol(matrix.noise.list[[i]])) {
            matInfo[, bootstrap.indices[[i]][j]] <- matInfo[, 
                bootstrap.indices[[i]][j]] + matrix.noise.list[[i]][, 
                j]
        }
    }

    if (!is.null(bootstrap.weights)) {
        for (i in 1:length(bootstrap.indices)) {
            temp = c()
            for (j in 1:length(bootstrap.indices[[i]])) {
                temp <- c(temp, rep(x = bootstrap.indices[[i]][j], 
                  each = bootstrap.weights[[i]][[j]]))
            }

            bootstrap.indices[[i]] <- temp
        }
    }
    bootstrap.indices[[length(bootstrap.indices) + 
        1]] <- 1:ncol(R.input)
    weightsBootSelected <- table(unlist(bootstrap.indices))
    weightsBootSelected <- weightsBootSelected - 1

    for (i in 1:ncol(matInfo)) {
        matInfo[, i] <- matInfo[, i]/weightsBootSelected[[i]]
    }

    return(matInfo)
}

## Evaluate the correlation between theta true
## and theta estimated
getCorrelations <- function(theta.true.list, theta.estimated.list) {
    pearsonCorrelations <- c()
    spearmanCorrelations <- c()
    kendallCorrelations <- c()
    for (i in 1:length(theta.true.list)) {
        pearsonCorrelations <- c(pearsonCorrelations, 
            cor(theta.true.list[[i]], theta.estimated.list[[i]], 
                method = "pearson"))
        spearmanCorrelations <- c(spearmanCorrelations, 
            cor(theta.true.list[[i]], theta.estimated.list[[i]], 
                method = "spearman"))
        kendallCorrelations <- c(kendallCorrelations, 
            cor(theta.true.list[[i]], theta.estimated.list[[i]], 
                method = "kendall"))
    }

    allCorrelation <- data.frame(pearson = pearsonCorrelations, 
        spearman = spearmanCorrelations, kendall = kendallCorrelations)
    return(allCorrelation)
}

# Get the minimum, maximum, median and mean times
# of the algorithm
getTimes <- function(listData) {
    allTimes <- c()
    for (i in 1:length(listData)) {
        allTimes <- c(allTimes, listData[[i]]$time)
    }
    minimumValue <- min(allTimes)
    maximumValue <- max(allTimes)
    meanValue <- mean(allTimes)
    medianValue <- median(allTimes)
    toReturn <- c(min = minimumValue, max = maximumValue, 
        mean = meanValue, median = medianValue)
    return(toReturn)
}

# normalization function
scaleFun <- function(x) {
    return(scale(x)[1:length(x)])
}
