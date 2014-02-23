starchartBuildTree <-
function(x, y, splittable) {
    # A Starchart tree (of the "starchart" class) is stored in a list with the
    # following named vectors; these vectors mimic a struct array due to the
    # lack of support for native struct and pointer types in R.
    # nodes: The name of the parameter this node splits on, "" for leaf nodes.
    # thresholds: The split threshold, NA for leaf nodes.
    # lefts: This node's left child's index in nodes[], NA for leaf nodes.
    # rights: This node's right child's index in nodes[], NA for leaf nodes.
    # erros: The sum of squared errors of this node.
    # x and y: The original data sets used to generate the tree.
    nodes = vector()
    thresholds = vector()
    lefts = vector()
    rights = vector()
    errors = vector()
    
    # Start the iterative process of building a regression tree.
    # All samples are treated as a single leaf node at the beginning.
    nodes = c(nodes, "")
    thresholds = c(thresholds, NA)
    lefts = c(lefts, NA)
    rights = c(rights, NA)
    errors = c(errors, NA)
    # This is the list of outstanding work items.
    leavesToSplit = list(list(id = length(nodes), x = x, y = y))
    while (length(leavesToSplit) > 0) {
        # As long as there are work items left, start working on the first one.
        id = leavesToSplit[[1]]$id
        result = starchartSplitLeaf(leavesToSplit[[1]]$x, leavesToSplit[[1]]$y,
                                    splittable)
        # No matter whether the split succeeds, remove the work item first.
        leavesToSplit = leavesToSplit[-1]
        # If the split is unsuccessful, nothing needs to be done.
        if (is.null(result)) next;
        # If the split is succeful, insert child nodes and generate work items.
        nodes[id] = result$splitParameter
        thresholds[id] = result$splitThreshold
        errors[id] = result$currentError
        # Add the left child node.
        nodes = c(nodes, "")
        thresholds = c(thresholds, NA)
        lefts = c(lefts, NA)
        lefts[id] = length(nodes)
        rights = c(rights, NA)
        fit = starchartFitModelToLeaf(result$leftX, result$leftY)
        errors = c(errors, starchartGetLeafModelError(fit))
        leavesToSplit[[length(leavesToSplit) + 1]] = list(id = length(nodes),
                                                          x = result$leftX,
                                                          y = result$leftY)
        # Add the right child node.
        nodes = c(nodes, "")
        thresholds = c(thresholds, NA)
        lefts = c(lefts, NA)
        rights = c(rights, NA)
        rights[id] = length(nodes)
        fit = starchartFitModelToLeaf(result$rightX, result$rightY)
        errors = c(errors, starchartGetLeafModelError(fit))
        leavesToSplit[[length(leavesToSplit) + 1]] = list(id = length(nodes),
                                                          x = result$rightX,
                                                          y = result$rightY)
    }
    list(nodes = nodes, thresholds = thresholds,
         lefts = lefts, rights = rights, errors = errors,
         x = x, y = y)
}
