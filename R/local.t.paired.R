local.t.paired <-
function (X.mat, y.vec, ...)
{
    if ((sum(sort(abs(y.vec[y.vec < 0])) != sort(y.vec[y.vec >
        0])) > 0) | length(unique(abs(y.vec))) != length(y.vec)/2) {
        stop("y.vec is improperly defined for a paired t-test",
            call. = FALSE)
    }
    return(function(data, vector = y.vec, ...) {
        x <- as.matrix(data[, vector > 0][, order(vector[vector >
            0])])
        y <- as.matrix(data[, vector < 0][, order(-vector[vector <
            0])])
        n <- dim(x)[[2]]
        x.m <- x %*% rep(1/n, n)
        y.m <- y %*% rep(1/n, n)
        x.res <- x - x.m %*% rep(1, n)
        y.res <- y - y.m %*% rep(1, n)
        t <- (x.m - y.m) * (sqrt(n * (n - 1)/((x.res - y.res)^2 %*%
            rep(1, n))))
        return(as.numeric(t))
    })
}
