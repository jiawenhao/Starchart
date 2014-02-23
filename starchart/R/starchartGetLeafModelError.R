starchartGetLeafModelError <-
function(fit) {
    # NOTE: Change the expression below to use other error metrics.
    sum(fit$residuals^2)
}
