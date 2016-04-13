## fxns.R

library(sandwich);

## **************************************************************
## integrative gwas
## **************************************************************
## Y = n x 1 outcome
## G = n x p genes
## S = n x q SNPs
## X = n x r clinical covariates
## family = type of outcome; same as glm argument
igwas <- function(Y,G,S,X=NULL,family)
{
    ## make everything into a matrix
    Y <- as.matrix(Y); G <- as.matrix(G); S <- as.matrix(S);
    if(is.null(X[1]))
    { X <- matrix(NA,nrow=nrow(S),ncol=0); } else
    { X <- as.matrix(X); }

    ## remove missing data
    cc <- complete.cases(Y,G,S,X);
    Y <- as.matrix(Y[cc,]);
    G <- as.matrix(G[cc,]);
    S <- as.matrix(S[cc,]);
    X <- as.matrix(X[cc,]);
    n <- sum(cc); p <- ncol(G); q <- ncol(S); r <- ncol(X);

    ## fit outcome model and transcript model
    if(r==0){ fit.o <- glm(Y~G,family=family); } else
    { fit.o <- glm(Y~G+X,family=family); }
    ystar <- G%*%coef(fit.o)[2:(1+p)];
    if(r==0){ fit.t <- lm(ystar~S); } else
    { fit.t <- lm(ystar~S+X); }

    ## calculate covariance matrix of \hat{\beta}_S
    Ainv <- bread(fit.o); ## already scaled by n
    covs <- cbind(1,S,X);
    B <- t(covs)%*%G/n;
    B <- cbind(rep(0,nrow(B)),B,
               matrix(0,nrow=nrow(B),ncol=r));
    Cinv <- solve(t(covs)%*%covs/n);
    ee <- cbind(estfun(fit.o),estfun(fit.t));
    V <- t(ee)%*%ee/n;
    V11 <- V[1:(1+p+r),1:(1+p+r)];
    V12 <- V[1:(1+p+r),(2+p+r):ncol(V)];
    V21 <- V[(2+p+r):ncol(V),1:(1+p+r)];
    V22 <- V[(2+p+r):ncol(V),(2+p+r):ncol(V)];

    ## covariance matrix of transcript model estimates
    Sigma <- (Cinv%*%(B%*%Ainv%*%V11%*%Ainv%*%t(B)+
                      V21%*%Ainv%*%t(B)+
                      B%*%Ainv%*%V12+
                      V22)%*%Cinv)/n;

    ## p-values for SNPs
    b <- coef(fit.t)[2:(1+q)];
    b.cov <- Sigma[2:(1+q),2:(1+q)];
    X2 <- t(b)%*%solve(b.cov)%*%b;
    p <- 1-pchisq(X2,q);
    return(list(b=b,b.cov=b.cov,X2=X2,p=p));
}

expit <- function(x) { 1/(1+exp(-x)); } # pi()
## Y = n x 1 outcome
## P = marginal prevalence of Y
## G = n x p genes
## S = n x q SNPs
## X = n x r clinical covariates
## family = type of outcome; same as glm argument
igwas.cc <- function(Y,P,G,S,X=NULL,family)
{
    ## make everything into a matrix
    Y <- as.matrix(Y); G <- as.matrix(G); S <- as.matrix(S);
    if(is.null(X[1]))
    { X <- matrix(NA,nrow=nrow(S),ncol=0); } else
    { X <- as.matrix(X); }

    ## remove missing data
    cc <- complete.cases(Y,G,S,X);
    Y <- as.matrix(Y[cc,]);
    G <- as.matrix(G[cc,]);
    S <- as.matrix(S[cc,]);
    X <- as.matrix(X[cc,]);
    p <- ncol(G); q <- ncol(S); r <- ncol(X);
    n1 <- sum(Y==1); n0 <- sum(Y==0); n <- n1+n0;

    ## fit outcome model and transcript model
    if(r==0){ fit.o <- glm(Y~G,family=family); } else
    { fit.o <- glm(Y~G+X,family=family); }
    ystar <- G%*%coef(fit.o)[2:(1+p)];
    w1 <- P*n/n1; w0 <- (1-P)*n/n0;
    ws <- Y; ws[Y==1] <- w1; ws[Y==0] <- w0;
    if(r==0){ fit.t <- lm(ystar~S,weights=ws); } else
    { fit.t <- lm(ystar~S+X,weights=ws); }
    
    ## calculate covariance matrix of \hat{\beta}_S
    Ainv <- bread(fit.o); ## already scaled by n
    covs <- cbind(1,S,X);
    B <- w1/n*t(covs[Y==1,])%*%as.matrix(G[Y==1,])+
         w0/n*t(covs[Y==0,])%*%as.matrix(G[Y==0,]);
    B <- cbind(rep(0,nrow(B)),
               B,
               matrix(0,nrow=nrow(B),ncol=r));
    Cinv <- solve(w1/n*t(covs[Y==1,])%*%covs[Y==1,]+
                  w0/n*t(covs[Y==0,])%*%covs[Y==0,]);
    ee <- cbind(estfun(fit.o),estfun(fit.t));
    V <- t(ee)%*%ee/n;
    V11 <- V[1:(1+p+r),1:(1+p+r)];
    V12 <- V[1:(1+p+r),(2+p+r):ncol(V)];
    V21 <- V[(2+p+r):ncol(V),1:(1+p+r)];
    V22 <- V[(2+p+r):ncol(V),(2+p+r):ncol(V)];

    ## covariance matrix of transcript model estimates
    Sigma <- (Cinv%*%(B%*%Ainv%*%V11%*%Ainv%*%t(B)+
                      V21%*%Ainv%*%t(B)+
                      B%*%Ainv%*%V12+
                      V22)%*%Cinv)/n;

    ## p-values for SNPs
    b <- coef(fit.t)[2:(1+q)];
    b.cov <- Sigma[2:(1+q),2:(1+q)];
    X2 <- t(b)%*%solve(b.cov)%*%b;
    p <- 1-pchisq(X2,q);
    return(list(b=b,b.cov=b.cov,X2=X2,p=p));
}
