predict.starchart <-
function(object, data = NULL, ...) {
    p = function(sample) {
        # Go down the tree and find the leaf node this sample falls in.
        index = 1
        while(object$nodes[index] != "") {
            if (sample[object$nodes[index]] <= object$thresholds[index])
                index = object$lefts[index]
            else
                index = object$right[index]
        }
        # Return the average response of all samples in the leaf node.
        return(mean(starchartGetNode(object, index)$y))
    }
    if (is.null(data))
        x = object$x
    else {
        if (!is.data.frame(data))
            cat("The input data should be a data frame!")
        x = as.data.frame(data)
    }
    y = apply(x, 1, p)
    return(y)
}
