
bootstrap_gurobi <- function(R.input, b, num.boot, 
    type, bootstrap.type, nCore = (detectCores() - 
        1), method = 0) {

    i <- NULL

    if (!requireNamespace("gurobi", quietly = TRUE)) {
        stop("Package gurobi needed for this function to work. Please install it.", 
            call. = FALSE)
    }


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
    # if(bootstrap.type ==
    # 'poisson.online.bootstrap'){ temp <-
    # generate.poisson.online.bootstrap.samples(R.input,
    # num.boot) boots <- temp$boots indices <-
    # temp$indices weights <- temp$weights }

    cl <- makeCluster(nCore)
    registerDoParallel(cl)

    bootEstimation <- foreach(i = 1:num.boot, .packages = c("gurobi", 
        "Matrix")) %dopar% {
        if (type == "fullQuadratic") {
            res <- GUROBISparsetest(boots[[i]], b, 
                method, weights[[i]])
            noiseMatrix <- getNoise(res$result$x, p)
        }
        if (type == "fullLinear") {
            res <- GUROBILinearSparseTest(boots[[i]], 
                b, method, weights[[i]])
            noiseMatrix <- getNoise(res$result$x, p)
        }
        if (type == "restrictedLinear") {
            res <- GUROBILinearMediumSparsetest(boots[[i]], 
                b, method, nLevel = 1, weights[[i]])
            noiseMatrix <- getNoise(res$result$x, p)
        }
        if (type == "restrictedQuadratic") {
            res <- GUROBIMediumSparsetest(boots[[i]], 
                b, method, nLevel = 1, weights[[i]])
            noiseMatrix <- getNoise(res$result$x, p)
        }

        theta_all <- scaleFun(getTheta(res$result$x, 
            p))
        # theta_all <- getTheta(res$result$x, p)
        list(theta = theta_all, noiseMatrix = noiseMatrix)
    }

    stopCluster(cl)

    thetas_norm_all <- list()
    matrix_noise_all <- list()
    for (i in 1:length(bootEstimation)) {
        thetas_norm_all[[i]] <- bootEstimation[[i]]$theta
        matrix_noise_all[[i]] <- bootEstimation[[i]]$noiseMatrix
    }

    estimate <- gather.results(thetas_norm_all)


    estimate.matrix <- gather.matrix.results(matrix.noise.list = matrix_noise_all, 
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
