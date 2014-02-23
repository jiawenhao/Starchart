starchartSplitLeaf <-
function(x, y, splittable) {
    # Theoretically, a leaf node can contain as few as one sample.
    # Pratically, it's probably better to enforce a minimum size requirement to
    # avoid over-growing a tree.
    if (nrow(x) < splittable) return(NULL)
    
    # Only parameters with at least two unique values are eligible for splits.
    parameters = names(x)[apply(x, 2, function(p) length(unique(p)) > 1)]
    if (length(parameters) == 0) return(NULL)
    
    fit = starchartFitModelToLeaf(x, y)
    # For each parameter, pick the split threshold that gives the maximum total
    # error reduction.
    # Each parameter has a number of thresholds to split on.
    # These thresholds lead to a number of (error) reductions.
    # The largest reduction is put in errorReductions and the corresponding
    # threshold is put in splitThresholds.
    splitThresholds = vector()
    errorReductions = vector()
    for (parameter in parameters) {
        if (is.factor(x[[parameter]])) {
            # TODO: Implement the logic for categorical parameters here.
            # For now, just ignore categorical parameters.
            errorReductions = c(errorReductions, -Inf)
            splitThresholds = c(splitThresholds, NA)
        } else {            
            # For a numerical parameter with n unique sampled values, there
            # are (n - 1) ways to divide the samples into two sets.
            # Note the apply() above has ensured n >= 2 for every parameter.
            values = x[[parameter]]
            thresholds = sort(unique(values[values < max(values)]))
            reductions = vector()
            
            # For each threshold value, tentatively split the leaf node into
            # two child nodes (left and right) based on that value and compute
            # the resulting total error reduction.
            for (threshold in thresholds) {
                leftX = x[values <= threshold, , drop = FALSE]
                leftY = y[values <= threshold]
                rightX = x[values > threshold, , drop = FALSE]
                rightY = y[values > threshold]
                stopifnot(nrow(leftX) > 0 && nrow(rightX) > 0)
                
                leftFit = starchartFitModelToLeaf(leftX, leftY)
                rightFit = starchartFitModelToLeaf(rightX, rightY)
                reductions = c(reductions,
                               starchartGetLeafModelError(fit) -
                              (starchartGetLeafModelError(leftFit) +
                               starchartGetLeafModelError(rightFit)))
            }
            errorReductions = c(errorReductions, max(reductions))
            splitThresholds = c(splitThresholds,
                                thresholds[which.max(reductions)])
        }
    }
    
    # The split attempt fails if no splits can reduce the total error at all.
    if (max(errorReductions) <= 0) return(NULL)
    
    splitParameter = parameters[which.max(errorReductions)]
    splitThreshold = splitThresholds[which.max(errorReductions)]
    if (is.factor(x[[splitParameter]])) {
        # TODO: Fill in the logic for categorical parameters.
        # For now, this should never happen.
        stopifnot(FALSE)
    } else {
        leftX = x[x[[splitParameter]] <= splitThreshold, , drop = FALSE]
        leftY = y[x[[splitParameter]] <= splitThreshold]
        rightX = x[x[[splitParameter]] > splitThreshold, , drop = FALSE]
        rightY = y[x[[splitParameter]] > splitThreshold]
    }
    stopifnot(nrow(leftX) > 0 && nrow(rightX) > 0)
    list(splitParameter = splitParameter, splitThreshold = splitThreshold,
         currentError = starchartGetLeafModelError(fit),
         leftX = leftX, leftY = leftY,
         rightX = rightX, rightY = rightY)
}
