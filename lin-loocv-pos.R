library(glmnet)


logp_       = seq(0,5)
plogstep    = 10**(1/5)
p_          = floor(120* plogstep**(logp_))
delta       = 1.5 # n/p
rho         = 0.2 # k/p
alpha_elnet = 1
errsd       = 1

pos = FALSE

lolims = ifelse(pos, 0, -Inf)

# Simulation



loocv_linear <- function(X, y, lambda0, alpha0){
  
  n <- nrow(X)
  
  pred_i <-function(i){
    
    X_i <- X[-i,]
    y_i <- y[-i]
    
    myfit_i <- glmnet(X_i, y_i, lambda = lambda0, alpha = alpha0, intercept = FALSE, family= "gaussian", lower.limits = lolims)
    
    beta_hat_i <- as.matrix(myfit_i$beta)[,1]
    
    dimnames(beta_hat_i) <- NULL

    df = sum(beta_hat_i != 0)/length(beta_hat_i)

    if(i %% 100 == 0){ cat(sprintf("p=%4d|n=%4d|i=%4d|df=%.2f|\n", 
                      ncol(X), n, i, df))
	}

    
    return((y[i]-sum(X[i,]*beta_hat_i))^2/2)
    
    
  }
  return(unlist(lapply(seq(1:n), pred_i)))
}


sampler = function(p, alpha_elnet){

  n = floor(p * delta)
  k = floor(p * rho)


  beta_star = rep(0,p)
  beta_star[1:k] = 1#abs(rexp(k)*(2*rbinom(k,1,0.5)-1) * sqrt(n/(2*k)))
  X = matrix(rnorm(n*p)/sqrt(n), n, p)
  
  y = X%*%beta_star + rnorm(n, sd = errsd)
  
  lambdaS = 1/n
  
  myfit <- glmnet(X, y, lambda = lambdaS, alpha = alpha_elnet, intercept = FALSE, family= "gaussian", lower.limits = lolims)
  
  beta_hat <- as.matrix(myfit$beta)[,1]
  
  dimnames(beta_hat) <- NULL
  
  lo_errs = loocv_linear(X, y, lambda0 = lambdaS, alpha0 = alpha_elnet)
  
  lo = sum(lo_errs)/n
  
  outErr = (errsd^2 + sum((beta_star-beta_hat)^2)/n)/2  
  return(c(outErr, lo))
}


sampler_all = function(){
	results = matrix(NA, 0, 3)
	for(i in 1:length(p_)){ results = rbind(results, c(p_[i], sampler(p_[i], alpha_elnet)))
		      print(results[i,])
		    }

	return(c(results))}