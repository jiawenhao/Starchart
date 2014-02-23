starchartGetNode <-
function(tree, index) {
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
