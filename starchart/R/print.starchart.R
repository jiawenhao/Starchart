print.starchart <-
function(x, ...) {
    cat("Printing a Starchart tree, one node per row:\n")
    print(data.frame(splitParameter = x$nodes, splitThreshold = x$thresholds,
                     leftChild = x$lefts, rightChild = x$rights))
}
