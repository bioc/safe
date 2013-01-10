error.FWER.WY <-
function (P.matrix) {
    order <- order(order(P.matrix[1, ]))
    P.matrix <- P.matrix[, order(P.matrix[1, ])]
    num.cat <- dim(P.matrix)[[2]]
    num.perms <- dim(P.matrix)[[1]]
    WY <- rep(NA, num.cat)
    for (i in 1:(num.cat - 1)) {
        WY[i] <- sum(apply(P.matrix, 1, min) <= P.matrix[1, 1])/num.perms
        P.matrix <- P.matrix[, -1]
    }
    WY[num.cat] <- sum(P.matrix <= P.matrix[1])/num.perms
    for (i in 2:num.cat) WY[i] <- max(WY[i], WY[i - 1])
    return(WY[order])
}
