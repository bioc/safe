error.FWER.Holm <-
function (P.matrix) {
    num.cat <- dim(P.matrix)[[2]]
    order <- order(order(P.matrix[1, ]))
    H <- P.matrix[1, order(P.matrix[1, ])] * (num.cat - 0:(num.cat -
        1))
    for (i in 2:num.cat) H[i] <- max(H[i], H[i - 1])
    return((H^(H <= 1))[order])
}
