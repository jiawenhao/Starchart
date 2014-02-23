###############################################################################
# Copyright 2014: Wenhao Jia <wjia@princeton.edu>, Princeton University
# All rights reserved.
# For documentation, see http://www.princeton.edu/~wjia/starchart/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright notice, 
#       this list of conditions, and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions, and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Princeton University nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY PRINCETON UNIVERSITY "AS IS" AND ANY EXPRESS OR 
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
# EVENT SHALL PRINCETON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

###########################################################################
# Customizable functions that help build the tree.
###########################################################################

# By default, Starchart uses the mean measured response of all samples in a
# leaf node as the predicted response for any design drawn from the part of
# the design space corresponding to that leaf node.
# In other words, the default model is a linear regression model with only an
# intercept term.
starchartFitModelToLeaf = function(x, y) {
    data = data.frame(x, response = y)
    # NOTE: Change the lm model below to customize leaf node models.
    lm(response ~ 1, data)
}

# By default, Starchart uses sum of squared errors (SSE) to assess how good a
# particular fit is.
starchartGetLeafModelError = function(fit) {
    # NOTE: Change the expression below to use other error metrics.
    sum(fit$residuals^2)
}

# For samples in a given leaf node, find the best way (i.e. highest reduction
# in total error) to split the leaf node into two child nodes.
# Return a list containing split results when the attempt succeeds.
# Return NULL when the attempt fails, e.g. when no beneficial splits exist.
starchartSplitLeaf = function(x, y, splittable) {
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

###############################################################################
# The main function starts here.
###############################################################################

# This is the main function that builds and returns a Starchart tree.
# x is a data frame containing all samples represented as parameter values.
# y is a vector of numeric measurements (i.e. responses) of the samples.
starchartBuildTree = function(x, y, splittable) {
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

###############################################################################
# Essential functions that complete a tree's functionality.
###############################################################################

predict.starchart = function(object, data = NULL, ...) {
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
    
# Borrow R's built-in dendrogram support to plot a Starchart tree.
plot.starchart = function(x, edge.labels = FALSE, height.sse = FALSE, 
                          top.splits = 0, ...) {
    # The Starchart tree was generated top-down. However, the hierarchical
    # cluster class (hclust) in R, which is the base of dendrograms, requires
    # a bottom-up order in describing hierarchical relations. Thus this
    # recursive function is used to convert the node order.
    # See ?hclust for definitions of its components.
    # The index-th node in tree is the node currently being processed.
    addNodeToHclust = function(hc, tree, index, limit) {
        # This function should never be directly invoked on leaf nodes.
        stopifnot(tree$nodes[index] != "")
        # Process the left child node.
        if (tree$nodes[tree$lefts[index]] == "") {
            node = starchartGetNode(tree, tree$lefts[index])
            label = sprintf("%d [%.2f]", nrow(node$x), mean(node$y))
            hc$labels = c(hc$labels, label)
            hc$leaftexts = c(hc$leaftexts, paste(tree$nodes[index], "<=",
                                                 tree$thresholds[index]))
            leftID = -length(hc$labels)
        } else {
            leftHeight = getUniqueHeight(tree, tree$lefts[index], limit)
            if (leftHeight <= 0) {
                # Treat the left non-leaf node as a leaf node.
                node = starchartGetNode(tree, tree$lefts[index])
                label = sprintf("%d [%.2f]", nrow(node$x), mean(node$y))
                hc$labels = c(hc$labels, label)
                hc$leaftexts = c(hc$leaftexts, paste(tree$nodes[index], "<=",
                                                     tree$thresholds[index]))
                leftID = -length(hc$labels)
            } else {
                hc = addNodeToHclust(hc, tree, tree$lefts[index], limit)
                hc$edgetexts[leftHeight] = 
                    paste(tree$nodes[index], "<=", tree$thresholds[index])
                leftID = nrow(hc$merge)
            }
        }
        # Process the right child node.
        if (tree$nodes[tree$rights[index]] == "") {
            node = starchartGetNode(tree, tree$rights[index])
            label = sprintf("%d [%.2f]", nrow(node$x), mean(node$y))
            hc$labels = c(hc$labels, label)
            hc$leaftexts = c(hc$leaftexts, paste(tree$nodes[index], ">",
                                                 tree$thresholds[index]))
            rightID = -length(hc$labels)
        } else {
            rightHeight = getUniqueHeight(tree, tree$rights[index], limit)
            if (rightHeight <= 0) {
                # Treat the right non-leaf node as a leaf node.
                node = starchartGetNode(tree, tree$rights[index])
                label = sprintf("%d [%.2f]", nrow(node$x), mean(node$y))
                hc$labels = c(hc$labels, label)
                hc$leaftexts = c(hc$leaftexts, paste(tree$nodes[index], ">",
                                                     tree$thresholds[index]))
                rightID = -length(hc$labels)
            } else {
                hc = addNodeToHclust(hc, tree, tree$rights[index], limit)
                hc$edgetexts[rightHeight] = 
                    paste(tree$nodes[index], ">", tree$thresholds[index])
                rightID = nrow(hc$merge)
            }
        }
        hc$merge = rbind(hc$merge, c(leftID, rightID))
        if (height.sse)
            # If the height.sse option is on, splits are plotted at heights
            # proportional to their respective error reductions.
            hc$height = c(hc$height, tree$errors[index])
        else {
            # Otherwise, they are plotted in the same vertical order but
            # evenly spaced.
            hc$height = c(hc$height, getUniqueHeight(tree, index, limit))
        }
        return(hc)
    }
    
    # Recursively apply labels to internal edges, marking split conditions.
    labelEdges = NULL
    local({
        i = 0
        labelEdges <<- function(x, edgetexts, leaftexts) {
            # This option is valid only when height.sse is turned off.
            if (height.sse) return(x)
            if (is.leaf(x)) {
                i <<- i + 1
                attributes(x)$edgetext = leaftexts[i]
                attributes(x)$edgePar = c(p.border = NA, t.cex = 0.75)
                attributes(x)$nodePar = c(pch = NA, lab.cex = 0.75)
                return(x)
            }
            attributes(x)$edgetext = edgetexts[attributes(x)$height]
            attributes(x)$edgePar = c(p.border = NA, t.cex = 0.75)
            return(x)
        }
    })
    
    # This function is used to find a unique height of a node.
    # It is particularly useful to resolve two nodes with equal SSEs.
    getUniqueHeight = function(tree, index, limit) {
        total = sum(tree$nodes != "")
        if (limit > 0)
            total = min(total, limit)
        all = which(tree$errors[index] == tree$errors[tree$nodes != ""])
        base = sum(tree$errors[index] <= tree$errors[tree$nodes != ""])
        if (length(all) > 1) {
            # There are multiple identical SSEs.
            offset = which(all == index) - 1
            height = total - (base + offset) + 1
        } else {
            height = total - base + 1
        }
        return(height)
    }
    
    # Prepare an empty hclust object to receive the conversion result.
    hc = list(merge = matrix(nrow = 0, ncol = 2), labels = vector(),
              height = vector(), edgetexts = vector(), leaftexts = vector())
    # Add nodes to hc in a breadth-first tree traversal order.
    hc = addNodeToHclust(hc, x, 1, top.splits)
    # Due to the way how hc is constructed, order[] is trivial to make.
    hc$order = seq(length(hc$labels))
    # Convert hc to a dengrogram object for further customization.
    class(hc) = "hclust"
    dg = as.dendrogram(hc)
    if (edge.labels)
        dg = dendrapply(dg, labelEdges, hc$edgetexts, hc$leaftexts)
    plot(dg, edge.root = TRUE, yaxt = 'n')
}

###############################################################################
# Helper functions that make manipulating a tree easier.
###############################################################################

# For the index-th node in the given tree, return its data selector as a
# concatinated string and the training/response samples contained in that node.
starchartGetNode = function(tree, index) {
    # Traverse the given tree and return the node matching the given index.
    nodesToVisit = list(list(index = 1, selector = "", x = tree$x, y = tree$y))
    while(length(nodesToVisit) > 0) {
        current = nodesToVisit[[1]]
        nodesToVisit = nodesToVisit[-1]
        if (current$index == index)
            return(current)
        # A leaf node doesn't have any child nodes to traverse.
        if (tree$nodes[current$index] == "")
            next
        # A non-leaf node generates new work from its child nodes.
        left = list(index = tree$lefts[current$index])
        left$selector = paste(current$selector, tree$nodes[current$index],
                              "<=", tree$thresholds[current$index])
        mask = current$x[tree$nodes[current$index]] <=
            tree$thresholds[current$index]
        left$x = current$x[mask, , drop = FALSE]
        left$y = current$y[mask]
        nodesToVisit[[length(nodesToVisit) + 1]] = left
        right = list(index = tree$rights[current$index])
        right$selector = paste(current$selector, tree$nodes[current$index],
                               ">", tree$thresholds[current$index])
        mask = current$x[tree$nodes[current$index]] >
            tree$thresholds[current$index]
        right$x = current$x[mask, , drop = FALSE]
        right$y = current$y[mask]
        nodesToVisit[[length(nodesToVisit) + 1]] = right
    }
    # Raise an exception if the given ID doesn't match any node.
    stopifnot(FALSE)
}

###############################################################################
# Boilerplate code to implement a nice S3 class interface.
###############################################################################

starchart = function(x, ...) {
    UseMethod("starchart")
}


# NOTE: A leaf node needs to have at least 10 samples to be considered for
# splitting. This can be overridden with a runtime argument: splittable.
starchart.default = function(x, y,
                             splittable = 10,
                             fitted = FALSE, ...) {
    x = as.data.frame(x)
    y = as.numeric(y)
    
    tree = starchartBuildTree(x, y, splittable = splittable)
    
    # (Optionally) compute a few additional fields for some generic functions.
    if (fitted) {
        tree$fitted.values = predict.starchart(tree, data = x)
        tree$residuals = y - tree$fitted.values
    }
    tree$call = match.call()
    
    class(tree) = "starchart"
    return(tree)
}

print.starchart = function(x, ...) {
    cat("Printing a Starchart tree, one node per row:\n")
    print(data.frame(splitParameter = x$nodes, splitThreshold = x$thresholds,
                     leftChild = x$lefts, rightChild = x$rights))
}

summary.starchart = function(object, ...) {
    result = list(call = object$call,
                  nodes = length(object$nodes),
                  leaves = sum(object$nodes == ""))
    if ("residuals" %in% names(object))
        result$residuals = object$residuals
    class(result) = "summary.starchart"
    return(result)
}

print.summary.starchart = function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("Number of nodes:", x$nodes, "\n")
    cat("Number of leaf nodes:", x$leaves, "\n")
    if ("residuals" %in% names(x)) {
        cat("Error distribution:\n")
        print(summary(x$residuals))
    }
}
