#quantile and cdf for truncated normal. Copied from R package msm

qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, log.p = FALSE) {
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    ret <- numeric(length(p))
    ret[p == 1] <- upper
    ret[p == 0] <- lower
    ret[p < 0 | p > 1] <- NaN
    ret[upper < lower] <- NaN
    ind <- (p > 0 & p < 1 & lower <= upper)
    if (any(ind)) {
        hind <- seq(along = p)[ind]
        h <- function(y) {
            (ptnorm(y, mean, sd, lower, upper) - p)[hind[i]]
        }
        ptmp <- numeric(length(p[ind]))
        for (i in 1:length(p[ind])) {
            interval <- c(-1, 1)
            while (h(interval[1]) * h(interval[2]) >= 0) interval <- interval + 
                c(-1, 1) * 0.5 * (interval[2] - interval[1])
            ptmp[i] <- uniroot(h, interval, tol = .Machine$double.eps)$root
        }
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) 
        warning("NaNs produced")
    ret
}

ptnorm <- function (q, mean = 0, sd = 1, lower = -Inf, upper = Inf, lower.tail = TRUE, log.p = FALSE) {
    ret <- numeric(length(q))
    ret[q < lower] <- 0
    ret[q > upper] <- 1
    ind <- q >= lower & q <= upper
    if (any(ind)) {
        denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, 
            sd)
        if (lower.tail) 
            qtmp <- pnorm(q, mean, sd) - pnorm(lower, mean, sd)
        else qtmp <- pnorm(upper, mean, sd) - pnorm(q, mean, 
            sd)
        if (log.p) 
            qtmp <- log(qtmp) - log(denom)
        else qtmp <- qtmp/denom
        ret[q >= lower & q <= upper] <- qtmp[ind]
    }
    ret
}
