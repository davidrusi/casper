setGeneric("asymmetryCheck", function(x,...) standardGeneric("asymmetryCheck"))
setMethod("asymmetryCheck", signature(x="ExpressionSet"), function(x, ...) { asymmetryCheck(exprs(x)) } )
setMethod("asymmetryCheck", signature(x="data.frame"), function(x, ...) { asymmetryCheck(as.matrix(x)) } )
setMethod("asymmetryCheck", signature(x="matrix"), function(x, ...) {
  m <- rowMeans(x); s <- rowSD(x)
  xsim <- s*matrix(rnorm(nrow(x)*ncol(x)), nrow=nrow(x), ncol=ncol(x)) + m
  sk <- skew(t(x))
  sk.sim <- skew(t(xsim))
  boxplot(data.frame(Observed=sk,Simulated=sk.sim),outline=FALSE,ylab='Asymmetry coefficient',...)
}
)


rowVar <- function(x,...) { return((rowMeans(x^2,...)-rowMeans(x,...)^2)*ncol(x)/(ncol(x)-1)) }
rowSD <- function(x,...) { return(sqrt(rowVar(x,...))) }


#skew function copied from R package skew
skew <- function (x, na.rm = TRUE, type = 3) {
    if (length(dim(x)) == 0) {
        if (na.rm) { x <- x[!is.na(x)] }
        sdx <- sd(x, na.rm = na.rm)
        mx <- mean(x)
        n <- length(x[!is.na(x)])
        switch(type, {
            skewer <- sqrt(n) * (sum((x - mx)^3, na.rm = na.rm)/(sum((x - mx)^2, na.rm = na.rm)^(3/2)))
        }, {
            skewer <- n * sqrt(n - 1) * (sum((x - mx)^3, na.rm = na.rm)/((n - 2) * sum((x - mx)^2, na.rm = na.rm)^(3/2)))
        }, {
            skewer <- sum((x - mx)^3)/(n * sd(x)^3)
        })
    }
    else {
        skewer <- rep(NA, dim(x)[2])
        if (is.matrix(x)) {
            mx <- colMeans(x, na.rm = na.rm)
        }
        else {
            mx <- apply(x, 2, mean, na.rm = na.rm)
        }
        sdx <- apply(x, 2, sd, na.rm = na.rm)
        for (i in 1:dim(x)[2]) {
            n <- length(x[!is.na(x[, i]), i])
            switch(type, {
                skewer[i] <- sqrt(n) * (sum((x[, i] - mx[i])^3, 
                  na.rm = na.rm)/(sum((x[, i] - mx[i])^2, na.rm = na.rm)^(3/2)))
            }, {
                skewer[i] <- n * sqrt(n - 1) * (sum((x[, i] - mx[i])^3, na.rm = na.rm)/((n - 2) * sum((x[, i] - mx[i])^2, na.rm = na.rm)^(3/2)))
            }, {
                skewer[i] <- sum((x[, i] - mx[i])^3, na.rm = na.rm)/(n * sdx[i]^3)
            })
        }
    }
    return(skewer)
}
