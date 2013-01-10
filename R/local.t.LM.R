local.t.LM <-
function (X.mat, y.vec, ...)
{
    n <- length(y.vec)
    return(function(data = X.mat, vector = y.vec, n2 = n, ...) {
        m.x <- (data %*% rep(1, n2))/n2
        s.xy <- data %*% vector
        v.x <- ((data^2 %*% rep(1, n2)) - n2 * m.x^2)/(n2 - 1)
        r <- (s.xy - n2 * m.x * mean(vector))/sqrt(v.x * var(vector))/(n2 -
            1)
        return(as.numeric(r * sqrt((n - 2)/(1 - r^2))))
    })
}
