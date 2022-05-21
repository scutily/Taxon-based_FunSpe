#install.packages("glmnet") 
#install.packages("msaenet")

### read data
# A ecosystem function that utilizing Glucose-1-phosphate 
EF <- read.csv("EF_G1P.csv")  

# Community omposition (ASV table) of prokaryotes.
composition_ASV <- read.csv("composition_ASV.csv")
row.names(composition_ASV) <- composition_ASV[,1]
composition_ASV <- as.matrix(composition_ASV[,-1])


#######################################################################
# Prepare the taxon-based functional specificity index using Msa-enet #
#######################################################################

## The orignal codes from library "msaenet" can not extract the minimum mse value.
## Therefore, we did a minor modification on their codes

###### begin of code modification (based on library "msaenet") ######

#' Automatic (parallel) parameter tuning for glmnet models
#'
#' @return Optimal model object, parameter set, and criterion value
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @references
#' Chen, Jiahua, and Zehua Chen. (2008).
#' Extended Bayesian information criteria for model selection with
#' large model spaces. \emph{Biometrika} 95(3), 759--771.
#'
#' @keywords internal
#' 


library(glmnet) 
library(msaenet)
library(foreach)
library(Matrix)
library(ncvreg)
library(hablar)


.is.msaenet <- function(x) "msaenet" %in% class(x)

.is.glmnet <- function(x) "glmnet" %in% class(x)

.is.ncvreg <- function(x) "ncvreg" %in% class(x)

.is.adaptive <- function(x)
  any(class(x) %in% c("msaenet.aenet", "msaenet.amnet", "msaenet.asnet"))

.is.multistep <- function(x)
  any(class(x) %in% c("msaenet.msaenet", "msaenet.msamnet", "msaenet.msasnet"))

# AIC
.aic <- function(deviance, df)
  deviance + (2L * df)

# BIC
.bic <- function(deviance, df, nobs)
  deviance + (log(nobs) * df)

# Extended BIC
.ebic <- function(deviance, df, nobs, nvar, gamma)
  deviance + (log(nobs) * df) + (2L * gamma * lchoose(nvar, df))

# deviance vector
.deviance <- function(model) {
  if (.is.glmnet(model)) return((1L - model$"dev.ratio") * model$"nulldev")
  if (.is.ncvreg(model)) return(model$"loss")
}

# degree of freedom vector
.df <- function(model) {
  if (.is.glmnet(model)) return(model$"df")
  if (.is.ncvreg(model)) {
    return(unname(colSums(
      as.matrix(abs(model$"beta"[-1L, ])) > .Machine$double.eps
    )))
  }
}

# number of observations in predictor matrix
.nobs <- function(model) {
  if (.is.glmnet(model)) return(model$"nobs")
  if (.is.ncvreg(model)) return(model$"n")
}

# dimensionality of predictor matrix
.nvar <- function(model) {
  if (.is.glmnet(model)) return(model$"dim"[[1L]])
  if (.is.ncvreg(model)) return(length(model$"penalty.factor"))
}


#####################################


msaenet.tune.glmnet.new <- function(
  x, y, family,
  alphas,
  tune,
  nfolds, rule,
  ebic.gamma,
  lower.limits, upper.limits,
  seed, parallel, ...) {
  
  if (tune == "cv") {
    if (!parallel) {
      model.list <- vector("list", length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] <- cv.glmnet(
          x = x, y = y, family = family,
          nfolds = nfolds, alpha = alphas[i],
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    } else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.glmnet(
          x = x, y = y, family = family,
          nfolds = nfolds, alpha = alphas,
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    }
    
    errors <- unlist(lapply(model.list, function(x) min(sqrt(x$"cvm"))))
    errors.min.idx <- which.min(errors) ###### Modified code
    min.errors=min(errors) ###### Modified code
    
    best.model <- model.list[[errors.min.idx]]
    
    best.alpha <- alphas[errors.min.idx]
    
    if (rule == "lambda.min") best.lambda <- best.model$"lambda.min"
    if (rule == "lambda.1se") best.lambda <- best.model$"lambda.1se"
    
    step.criterion <- errors[errors.min.idx]
  } else {
    if (!parallel) {
      model.list <- vector("list", length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] <- glmnet(
          x = x, y = y, family = family,
          alpha = alphas[i],
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    } else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        glmnet(
          x = x, y = y, family = family,
          alpha = alphas,
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    }
    
    if (tune == "aic") {
      ics.list <- mapply(
        .aic,
        deviance = lapply(model.list, .deviance),
        df = lapply(model.list, .df),
        SIMPLIFY = FALSE
      )
    }
    
    if (tune == "bic") {
      ics.list <- mapply(
        .bic,
        deviance = lapply(model.list, .deviance),
        df = lapply(model.list, .df),
        nobs = lapply(model.list, .nobs),
        SIMPLIFY = FALSE
      )
    }
    
    if (tune == "ebic") {
      ics.list <- mapply(
        .ebic,
        deviance = lapply(model.list, .deviance),
        df = lapply(model.list, .df),
        nobs = lapply(model.list, .nobs),
        nvar = lapply(model.list, .nvar),
        gamma = ebic.gamma,
        SIMPLIFY = FALSE
      )
    }
    
    ics <- sapply(ics.list, function(x) min(x))
    ics.min.idx <- which.min(ics)     ########## what is this?
    best.model <- model.list[[ics.min.idx]]
    
    best.alpha <- alphas[ics.min.idx]
    
    best.ic.min.idx <- which.min(ics.list[[ics.min.idx]])   
    best.lambda <- best.model$"lambda"[[best.ic.min.idx]]
    
    step.criterion <- ics.list[[ics.min.idx]][[best.ic.min.idx]]
  }
  
  list(
    "best.model" = best.model,
    "best.alpha" = best.alpha,
    "best.lambda" = best.lambda,
    "step.criterion" = step.criterion,
    "min.errors"=min.errors            ########## assign minimum mse as min.errors
  )
}

#' Select the number of adaptive estimation steps
#'
#' @return optimal step number
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @keywords internal

msaenet.tune.nsteps.glmnet <- function(
  model.list,
  tune.nsteps, ebic.gamma.nsteps) {
  
  nmods <- length(model.list)
  
  if (tune.nsteps == "max") {
    ics <- NULL
    best.step <- nmods
  } else {
    if (tune.nsteps == "aic") {
      ics <- .aic(
        deviance = sapply(model.list, .deviance),
        df = sapply(model.list, .df)
      )
    }
    
    if (tune.nsteps == "bic") {
      ics <- .bic(
        deviance = sapply(model.list, .deviance),
        df = sapply(model.list, .df),
        nobs = sapply(model.list, .nobs)
      )
    }
    
    if (tune.nsteps == "ebic") {
      ics <- .ebic(
        deviance = sapply(model.list, .deviance),
        df = sapply(model.list, .df),
        nobs = sapply(model.list, .nobs),
        nvar = sapply(model.list, .nvar),
        gamma = ebic.gamma.nsteps
      )
    }
    
    best.step <- which.min(ics)
  }
  
  list("best.step" = best.step, "ics" = ics)
}

#########################
#########################
#########################
msaenet.new <- function(
  x, y,
  family = c("gaussian", "binomial", "poisson", "cox"),
  init = c("enet", "ridge"),
  alphas = seq(0.05, 0.95, 0.05),
  tune = c("cv", "ebic", "bic", "aic"),
  nfolds = 5L, rule = c("lambda.min", "lambda.1se"),
  ebic.gamma = 1,
  nsteps = 2L,
  tune.nsteps = c("max", "ebic", "bic", "aic"),
  ebic.gamma.nsteps = 1,
  scale = 1,
  lower.limits = -Inf, upper.limits = Inf,
  penalty.factor.init = rep(1, ncol(x)),
  seed = 1001, parallel = FALSE, verbose = FALSE) {
  
  if (nsteps < 2L) stop("nsteps must be an integer >= 2")
  
  family <- match.arg(family)
  init <- match.arg(init)
  tune <- match.arg(tune)
  rule <- match.arg(rule)
  tune.nsteps <- match.arg(tune.nsteps)
  call <- match.call()
  
  best.alphas <- rep(NA, nsteps + 1L)
  best.lambdas <- rep(NA, nsteps + 1L)
  step.criterion <- rep(NA, nsteps + 1L)
  beta.list <- vector("list", nsteps + 1L)
  model.list <- vector("list", nsteps + 1L)
  min.errors <- rep(NA, nsteps + 1L) ################
  adapen.list <- vector("list", nsteps)
  
  if (verbose) cat("Starting step 1 ...\n")
  
  if (init == "enet") {
    model.cv <- msaenet.tune.glmnet.new(
      x = x, y = y, family = family,
      alphas = alphas,
      tune = tune,
      nfolds = nfolds, rule = rule,
      ebic.gamma = ebic.gamma,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      penalty.factor = penalty.factor.init,
      seed = seed, parallel = parallel
    )
  }
  
  if (init == "ridge") {
    model.cv <- msaenet.tune.glmnet.new(
      x = x, y = y, family = family,
      alphas = 0,
      tune = tune,
      nfolds = nfolds, rule = rule,
      ebic.gamma = ebic.gamma,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      penalty.factor = penalty.factor.init,
      seed = seed, parallel = parallel
    )
  }
  
  best.alphas[[1L]] <- model.cv$"best.alpha"
  best.lambdas[[1L]] <- model.cv$"best.lambda"
  step.criterion[[1L]] <- model.cv$"step.criterion"
  min.errors[[1L]]<- model.cv$"min.errors"    ########## assign minimum mse as min.errors
  
  
  
  model.list[[1L]] <- glmnet(
    x = x, y = y, family = family,
    alpha = best.alphas[[1L]],
    lambda = best.lambdas[[1L]],
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    penalty.factor = penalty.factor.init
  )
  
  if (.df(model.list[[1L]]) < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try a different parameter setting.")
  }
  
  bhat <- as.matrix(model.list[[1L]][["beta"]])
  if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))
  beta.list[[1L]] <- bhat
  
  # MSAEnet steps
  for (i in 1L:nsteps) {
    adpen.raw <- (pmax(abs(beta.list[[i]]), .Machine$double.eps))^(-scale)
    adapen.list[[i]] <- as.vector(adpen.raw)
    adpen.name <- rownames(adpen.raw)
    names(adapen.list[[i]]) <- adpen.name
    
    if (verbose) cat("Starting step", i + 1, "...\n")
    
    model.cv <- msaenet.tune.glmnet.new(
      x = x, y = y, family = family,
      alphas = alphas,
      tune = tune,
      nfolds = nfolds, rule = rule,
      ebic.gamma = ebic.gamma,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      seed = seed + i, parallel = parallel,
      penalty.factor = adapen.list[[i]]
    )
    
    best.alphas[[i + 1L]] <- model.cv$"best.alpha"
    best.lambdas[[i + 1L]] <- model.cv$"best.lambda"
    step.criterion[[i + 1L]] <- model.cv$"step.criterion"
    min.errors[[i+1L]]<- model.cv$"min.errors" ########## assign minimum mse as min.errors
    
    
    model.list[[i + 1L]] <- glmnet(
      x = x, y = y, family = family,
      alpha = best.alphas[[i + 1L]],
      lambda = best.lambdas[[i + 1L]],
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      penalty.factor = adapen.list[[i]]
    )
    
    if (.df(model.list[[i + 1L]]) < 0.5) {
      stop("Null model produced by the full fit (all coefficients are zero). Please try a different parameter setting.")
    }
    
    bhat <- as.matrix(model.list[[i + 1L]][["beta"]])
    if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))
    beta.list[[i + 1L]] <- bhat
  }
  
  # select optimal step
  post.ics <- msaenet.tune.nsteps.glmnet(
    model.list, tune.nsteps, ebic.gamma.nsteps
  )
  
  best.step <- post.ics$"best.step"
  post.criterion <- post.ics$"ics"
  
  msaenet.model <- list(
    "beta" = Matrix(beta.list[[best.step]], sparse = TRUE),
    "model" = model.list[[best.step]],
    "best.step" = best.step,
    "best.alphas" = best.alphas,
    "best.lambdas" = best.lambdas,
    "step.criterion" = step.criterion,
    "post.criterion" = post.criterion,
    "beta.list" = beta.list,
    "model.list" = model.list,
    "adapen.list" = adapen.list,
    "seed" = seed,
    "call" = call,
    "min.errors"=min.errors ################
  )
  
  class(msaenet.model) <- c("msaenet", "msaenet.msaenet")
  msaenet.model
}

###### end of code modification ######

###############################################
# Perform the Msa-enet based on modified code #
###############################################

#install.packages("doParallel")
library("doParallel")
registerDoParallel(detectCores()) ## speed up the computation

a <- c(1:9)/10   # Set alpha values from 0.1-0.9 for sparse regression (i.e., pure elastic net)


msaenet_list <- msaenet.new(composition_ASV,EF[,2],"gaussian","ridge" # Initial step = ridge
                           ,alpha=a,  nfolds=nrow(composition_ASV), # i.e., Leave-one-out cross-validation
                           rule="lambda.min",seed=1010,
                           nsteps = 10L, # 10 steps for variables selection
                           tune.nsteps = "ebic",parallel = TRUE) 
msaenet_mse <- msaenet_list$min.errors # Extract the mse for each step
msaenet_minmse_coefficients<- msaenet_list$beta.list[[which.min(msaenet_mse)]] # Select the step depending on the minimum mse
  
abs_faci_ASV <- length(which(msaenet_minmse_coefficients>0))  # Absolute number of facilitative ASV
abs_repr_ASV <- length(which(msaenet_minmse_coefficients<0))  # Absolute number of repressive ASV
prop_faci_ASV <- abs_faci_ASV/(abs_faci_ASV+abs_repr_ASV) # Proportion of facilitative ASV

specificity_taxon <- data.frame(abs_faci_ASV, abs_repr_ASV, prop_faci_ASV)

