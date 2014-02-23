plot.starchart <-
function(x, edge.labels = FALSE, height.sse = FALSE, 
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
