buildInequalityConstraintsFunFull <- function(rankMat, 
    b) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors

    A <- buildFullA(rankMat, p, n)  ## The A matrix for linear constraints
    IFF <- inequalityFactoryFunction(A, b)

    return(inequalityConstraintFunction = IFF)
}

buildInequalityConstraintsFunReduced <- function(rankMat, 
    b) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors

    A <- buildReducedA(rankMat)  ## The A matrix for linear constraints
    IFF <- inequalityFactoryFunction(A, b)

    return(inequalityConstraintFunction = IFF)
}

# build the starting point
startingPointFull <- function(rankMat) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors
    x0 <- rep(0, p + p * n)
    # x0[1:p] <- 1
    noiseToAdd <- rev(seq(from = 0, to = 1, length.out = p))
    index = p + 1
    for (i in 1:n) {
        col <- rankMat[, i]
        x0[index:(index + p - 1)] <- noiseToAdd[col]
        index <- index + p
    }

    return(x0)
}

startingPointReduced <- function(rankMat) {
    d <- dim(rankMat)
    p <- d[1]  #The number of objects
    n <- d[2]  #The number of assessors
    x0 <- rep(0, p + (p - 1) * n)
    # x0[1:p] <- 1
    noiseToAdd <- rev(seq(from = 0, to = 1, length.out = (p)))
    noiseToAdd <- noiseToAdd[1:(length(noiseToAdd) - 
        1)]
    index = p + 1
    for (i in 1:n) {
        col <- rankMat[, i]
        x0[index:(index + p - 2)] <- noiseToAdd
        index <- index + p - 1
    }

    return(x0)
}

# inequality factor function for
inequalityFactoryFunction <- function(A, b) {
    fun <- function(x) {
        return((A %*% x) + b)
    }
    return(fun)
}

# Build vector for objective function
quadraticObjectiveFunction <- function(p) {
    fun <- function(x) {
        return(sum(x[(p + 1):length(x)]^2))
    }
    return(fun)
}

linearObjectiveFunction <- function(p) {
    fun <- function(x) {
        return(sum(x[(p + 1):length(x)]))
    }
    return(fun)
}

buildFullA <- function(rankMat, p, n) {
    nConstraints <- n * ((p - 1) * p/2)
    nVariables <- p + p * n
    A <- matrix(0, nrow = nConstraints, ncol = nVariables)
    matrixIndices <- 1

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

                A[matrixIndices, posObjpos] = 1
                A[matrixIndices, posZIJpos] = 1
                A[matrixIndices, posObjless] = -1
                A[matrixIndices, posZIJless] = -1

                matrixIndices = matrixIndices + 1
            }  ## For local constraints
        }  ## for each position
    }  ## Each column rank constraints

    return(A)
}

buildReducedA <- function(rankMat) {
    p <- dim(rankMat)[1]
    n <- dim(rankMat)[2]

    nConstraints <- n * (p - 1)
    nVariables <- p + n * (p - 1)

    A <- matrix(0, nrow = nConstraints, ncol = p)

    matrixIndices <- 1

    for (index in 1:ncol(rankMat)) {
        rankI <- rankMat[, index]
        for (pos in 1:(p - 1)) {
            A[matrixIndices, which(rankI == pos)] = 1
            A[matrixIndices, which(rankI == (pos + 
                1))] = -1

            matrixIndices <- matrixIndices + 1
        }  ## for each position
    }  ## Each column rank constraints

    A <- cbind(A, diag(nrow = n * (p - 1), ncol = n * 
        (p - 1)))

    return(A)
}


SQPFullQuadraticTest <- function(R.input, b) {
    p <- dim(R.input)[1]
    n <- dim(R.input)[2]
    inequalityConstraints <- buildInequalityConstraintsFunFull(R.input, 
        -b)

    x0 <- startingPointFull(R.input)

    lb <- c(rep(0, p), rep(0, p * n))
    up <- c(rep(Inf, p + p * n))

    start_time <- Sys.time()
    res <- nloptr::slsqp(x0 = x0, fn = quadraticObjectiveFunction(p = p), 
        hin = inequalityConstraints, lower = lb, upper = up)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")
    return(list(result = res, time = runningTime))
}

SQPFullLinearTest <- function(R.input, b) {
    p <- dim(R.input)[1]
    n <- dim(R.input)[2]
    inequalityConstraints <- buildInequalityConstraintsFunFull(R.input, 
        -b)

    x0 <- startingPointFull(R.input)

    lb <- c(rep(0, p), rep(0, p * n))
    up <- c(rep(Inf, p + p * n))

    start_time <- Sys.time()
    res <- nloptr::slsqp(x0 = x0, fn = linearObjectiveFunction(p = p), 
        hin = inequalityConstraints, lower = lb, upper = up)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")
    return(list(result = res, time = runningTime))
}


SQPReducedLinearTest <- function(R.input, b) {
    p <- dim(R.input)[1]
    n <- dim(R.input)[2]
    inequalityConstraints <- buildInequalityConstraintsFunReduced(R.input, 
        -b)

    x0 <- startingPointReduced(R.input)

    lb <- c(rep(0, p), rep(0, (p - 1) * n))
    up <- c(rep(Inf, p + (p - 1) * n))

    start_time <- Sys.time()
    res <- nloptr::slsqp(x0 = x0, fn = linearObjectiveFunction(p = p), 
        hin = inequalityConstraints, lower = lb, upper = up)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")
    return(list(result = res, time = runningTime))
}


SQPReducedQuadraticTest <- function(R.input, b) {
    p <- dim(R.input)[1]
    n <- dim(R.input)[2]
    inequalityConstraints <- buildInequalityConstraintsFunReduced(R.input, 
        -b)

    x0 <- startingPointReduced(R.input)

    lb <- c(rep(0, p), rep(0, (p - 1) * n))
    up <- c(rep(Inf, p + (p - 1) * n))

    start_time <- Sys.time()
    res <- nloptr::slsqp(x0 = x0, fn = quadraticObjectiveFunction(p = p), 
        hin = inequalityConstraints, lower = lb, upper = up)
    end_time <- Sys.time()

    runningTime <- difftime(end_time, start_time, units = "secs")
    return(list(result = res, time = runningTime))
}
