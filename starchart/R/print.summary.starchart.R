print.summary.starchart <-
function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("Number of nodes:", x$nodes, "\n")
    cat("Number of leaf nodes:", x$leaves, "\n")
    if ("residuals" %in% names(x)) {
        cat("Error distribution:\n")
        print(summary(x$residuals))
    }
}
