buildSparseA <- function(rankMat, p, n) {
    nConstraints <- n * ((p - 1) * p/2)
    nVariables <- p + p * n

    NonZeroElement <- 4 * nConstraints
    iToSave <- rep(0, NonZeroElement)
    jToSave <- rep(0, NonZeroElement)
    xToSave <- rep(0, NonZeroElement)

    matrixIndices <- 1
    indexSave <- 1

    # index is the j-esima column of the noise
    for (index in 1:ncol(rankMat)) {
        rankI <- rankMat[, index]
        for (pos in 1:length(rankI)) {
            # print(paste('pos',pos, 'obj: ',
            # rankI[pos]))
            NlocalConstraints <- length(pos:p) - 1
            if (NlocalConstraints == 0) 
                break

            posObjpos <- which(rankI == pos)
            posZIJpos <- getNoiseIJ(which(rankI == 
                pos), index, p)

            for (a in (pos + 1):p) {
                posObjless <- which(rankI == a)
                posZIJless <- getNoiseIJ(which(rankI == 
                  a), index, p)

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posObjpos
                xToSave[indexSave] <- 1

                indexSave <- indexSave + 1

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posZIJpos
                xToSave[indexSave] <- 1

                indexSave <- indexSave + 1

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posObjless
                xToSave[indexSave] <- -1

                indexSave <- indexSave + 1

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posZIJless
                xToSave[indexSave] <- -1

                indexSave <- indexSave + 1

                matrixIndices = matrixIndices + 1
            }  ## For local constraints
        }  ## for each position
    }  ## Each column rank constraints

    A <- Matrix::spMatrix(nrow = nConstraints, ncol = nVariables, 
        i = iToSave, j = jToSave, x = xToSave)

    return(A)
}



buildSparseQ <- function(p, n, weights) {
    nElement <- p + p * n
    indices <- (p + 1):nElement
    if (is.null(weights)) {
        values <- rep(1, times = length(indices))
    } else {
        values <- rep(weights, each = p)
    }
    Q <- Matrix::spMatrix(nrow = nElement, ncol = nElement, 
        i = indices, j = indices, x = values)
    return(Q)
}


buildSparseAll <- function(rankMat, b, weights) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors

    Q <- buildSparseQ(p, n, weights)  ## The Q matrix for quadratic optmization

    A <- buildSparseA(rankMat, p, n)  ## The A matrix for linear constraints

    model <- list()
    model$A <- A
    model$Q <- Q
    model$rhs <- rep(b, (dim(A)[1]))
    model$sense <- rep(">", (dim(A)[1]))
    model$modelsense <- "min"
    model$vtype <- rep("C", dim(A)[2])
    model$lb <- 0
    model$up <- Inf

    return(model)
}

buildLinearSparseAll <- function(rankMat, b, weights = NULL) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors

    A <- buildSparseA(rankMat, p, n)  ## The A matrix for linear constraints
    if (is.null(weights)) {
        obj <- rep(1, ncol(A))
        obj[1:p] <- 0
    } else {
        obj <- c(rep(0, p), rep(weights, each = p))
    }

    model <- list()
    model$A <- A
    model$obj <- obj
    model$rhs <- rep(b, (dim(A)[1]))
    model$sense <- rep(">", (dim(A)[1]))
    model$modelsense <- "min"
    model$vtype <- rep("C", dim(A)[2])
    model$lb <- 0
    model$up <- Inf

    return(model)
}

GUROBISparsetest <- function(R.input, b, method, weights = NULL) {
    model <- buildSparseAll(R.input, b, weights)
    params <- list()
    params$Method <- method

    start_time <- Sys.time()
    res <- gurobi::gurobi(model, params = params)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")

    return(list(result = res, time = runningTime))
}

GUROBILinearSparseTest <- function(R.input, b, method, 
    weights = NULL) {
    model <- buildLinearSparseAll(R.input, b, weights)
    params <- list()
    params$Method <- method

    start_time <- Sys.time()
    res <- gurobi::gurobi(model, params = params)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")

    return(list(result = res, time = runningTime))
}


####### UTILS MEDIUM METHOD
buildMediumSparseA <- function(rankMat, nLevel) {
    p <- nrow(rankMat)
    n <- ncol(rankMat)

    tempVar <- 1:nLevel
    nLocalCon <- rev(c(tempVar, rep(nLevel, p - 1 - 
        length(tempVar))))

    nConstraints <- n * sum(nLocalCon)
    nVariables <- p + p * n

    NonZeroElement <- 4 * nConstraints
    iToSave <- rep(0, NonZeroElement)
    jToSave <- rep(0, NonZeroElement)
    xToSave <- rep(0, NonZeroElement)

    matrixIndices <- 1
    indexSave <- 1
    nLocalCon <- c(nLocalCon, 0)
    # index is the j-esima column of the noise
    for (index in 1:n) {
        rankI <- rankMat[, index]
        for (pos in 1:p) {
            # print(paste('pos',pos, 'obj: ',
            # rankI[pos]))
            NlocalConstraints <- nLocalCon[pos]
            if (NlocalConstraints == 0) 
                break

            posObjpos <- which(rankI == pos)
            posZIJpos <- getNoiseIJ(which(rankI == 
                pos), index, p)

            for (a in (pos + 1):(pos + nLocalCon[pos])) {
                posObjless <- which(rankI == a)
                posZIJless <- getNoiseIJ(which(rankI == 
                  a), index, p)

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posObjpos
                xToSave[indexSave] <- 1

                indexSave <- indexSave + 1

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posZIJpos
                xToSave[indexSave] <- 1

                indexSave <- indexSave + 1

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posObjless
                xToSave[indexSave] <- -1

                indexSave <- indexSave + 1

                iToSave[indexSave] <- matrixIndices
                jToSave[indexSave] <- posZIJless
                xToSave[indexSave] <- -1

                indexSave <- indexSave + 1

                matrixIndices <- matrixIndices + 1
            }  ## For local constraints
        }  ## for each position
    }  ## Each column rank constraints

    A <- Matrix::spMatrix(nrow = nConstraints, ncol = nVariables, 
        i = iToSave, j = jToSave, x = xToSave)

    return(A)
}


# buildMediumSparseQ <- function(p,n,weights){
# library('Matrix') nElement <- p + n*3*(p-2)
# indices <- (p+1):nElement if (is.null(weights))
# { values <- rep(1, times = length(indices)) }
# else { values <- rep(weights, each = p) } Q <-
# spMatrix(nrow = nElement,ncol=nElement,i =
# indices,j=indices,x=values) return(Q) }

buildMediumSparseAll <- function(rankMat, b, nLevel, 
    weights) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors

    Q <- buildSparseQ(p, n, weights)  ## The Q matrix for quadratic optmization

    A <- buildMediumSparseA(rankMat, nLevel)  ## The A matrix for linear constraints

    model <- list()
    model$A <- A
    model$Q <- Q
    model$rhs <- rep(b, (dim(A)[1]))
    model$sense <- rep(">", (dim(A)[1]))
    model$modelsense <- "min"
    model$vtype <- rep("C", dim(A)[2])
    model$lb <- 0
    model$up <- Inf

    return(model)
}


buildLineardMediumSparseAll <- function(rankMat, b, 
    nLevel, weights) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors

    A <- buildMediumSparseA(rankMat, nLevel)  ## The A matrix for linear constraints

    if (is.null(weights)) {
        obj <- rep(1, ncol(A))
        obj[1:p] <- 0
    } else {
        obj <- c(rep(0, p), rep(weights, each = p))
    }

    model <- list()
    model$A <- A
    model$obj <- obj
    model$rhs <- rep(b, (dim(A)[1]))
    model$sense <- rep(">", (dim(A)[1]))
    model$modelsense <- "min"
    model$vtype <- rep("C", dim(A)[2])
    model$lb <- 0
    model$up <- Inf

    return(model)
}


GUROBIMediumSparsetest <- function(R.input, b, method, 
    nLevel, weights = NULL) {
    model <- buildMediumSparseAll(R.input, b, nLevel, 
        weights)
    params <- list()
    params$Method <- method

    start_time <- Sys.time()
    res <- gurobi::gurobi(model, params = params)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")

    return(list(result = res, time = runningTime))
}

GUROBILinearMediumSparsetest <- function(R.input, b, 
    method, nLevel, weights = NULL) {
    model <- buildLineardMediumSparseAll(R.input, b, 
        nLevel, weights)
    params <- list()
    params$Method <- method

    start_time <- Sys.time()
    res <- gurobi::gurobi(model, params = params)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")

    return(list(result = res, time = runningTime))
}

