bootstrap_sqp <- function(R.input, b, num.boot, type, 
    bootstrap.type, nCore = (detectCores() - 1)) {

    i <- NULL

    p <- nrow(R.input)
    start_time <- Sys.time()
    if (bootstrap.type == "classic.bootstrap") {
        bootsRes <- generate.bootstrap.samples(R.input, 
            num.boot - 1)
        boots <- bootsRes$boots
        indices <- bootsRes$indices
        weights <- NULL
    }
    if (bootstrap.type == "poisson.bootstrap") {
        temp <- generate.poisson.bootstrap.samples(R.input, 
            num.boot)
        boots <- temp$boots
        indices <- temp$indices
        weights <- temp$weights
    }

    cl <- makeCluster(nCore)
    registerDoParallel(cl)

    bootEstimation <- foreach(i = 1:num.boot, .packages = "nloptr") %dopar% 
        {
            if (type == "fullQuadratic") {
                res <- SQPFullQuadraticTest(boots[[i]], 
                  b)
                noiseMatrix <- getNoise(res$result$par, 
                  p)
            }
            if (type == "fullLinear") {
                res <- SQPFullLinearTest(boots[[i]], 
                  b)
                noiseMatrix <- getNoise(res$result$par, 
                  p)
            }

            theta_all <- scaleFun(getTheta(res$result$par, 
                p))

            list(theta = theta_all, noiseMatrix = noiseMatrix)
        }
    stopCluster(cl)

    thetas_norm_all <- list()
    matrix_noie_all <- list()
    for (i in 1:length(bootEstimation)) {
        thetas_norm_all[[i]] <- bootEstimation[[i]]$theta
        matrix_noie_all[[i]] <- bootEstimation[[i]]$noiseMatrix
    }

    estimate <- gather.results(thetas_norm_all)


    estimate.matrix <- gather.matrix.results(matrix.noise.list = matrix_noie_all, 
        bootstrap.indices = indices, bootstrap.weights = weights, 
        R.input = R.input)

    end_time <- Sys.time()

    estimate <- cbind(id = factor(rownames(R.input)), 
        estimate)
    colnames(estimate) <- c("id", "signal.estimate", 
        "SE")

    allEstimations <- t(data.frame(thetas_norm_all))
    rownames(allEstimations) <- paste0("boot", 1:nrow(allEstimations))
    colnames(allEstimations) <- rownames(R.input)

    colnames(estimate.matrix) <- colnames(R.input)

    time_elapsed <- difftime(end_time, start_time, 
        units = "secs")

    return(list(estimation = estimate, estimatedMatrixNoise = estimate.matrix, 
        time = time_elapsed, allBootstraps = allEstimations))
}
