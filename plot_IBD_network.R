#!/usr/bin/env Rscript
library(igraph)

## additional shapes taken from igraph shapes manual
mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(
    x = coords[, 1],
    y = coords[, 2],
    bg = vertex.color,
    stars = cbind(vertex.size, vertex.size, vertex.size),
    add = TRUE,
    inches = FALSE
  )
}
# clips as a circle
add_shape("triangle",
          clip = shapes("circle")$clip,
          plot = mytriangle)

# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(
    coords[, 1],
    coords[, 2],
    vertex.color,
    vertex.size,
    norays,
    FUN = function(x, y, bg, size, nor) {
      symbols(
        x = x,
        y = y,
        bg = bg,
        stars = matrix(c(size, size / 2), nrow = 1, ncol = nor * 2),
        add = TRUE,
        inches = FALSE
      )
    }
  )
}
# no clipping, edges will be below the vertices anyway
add_shape(
  "star",
  clip = shape_noclip,
  plot = mystar,
  parameters = list(vertex.norays = 5)
)


.reorder.metadata <-
  function(IBD,
           metadata,
           pair1.col = "sample1",
           pair2.col = "sample2",
           sample.col = "Sample") {
    samples <- unique(c(IBD[, pair1.col], IBD[, pair2.col]))
    sample.order <- match(metadata[, sample.col], samples)
    actual.length <- length(na.omit(sample.order))
    
    metadata <- metadata[order(sample.order), ]
    metadata[1:actual.length, ]
  }


.get.labels <- function(metadata, label.col, label.palette = NULL) {
  labels <- metadata[, label.col]
  labels[is.na(labels)] <- NaN
  
  labels <- factor(labels)
  
  if (is.null(label.palette)) {
    labels
  } else {
    factor(labels, levels = names(label.palette))
  }
}


.generate.label.palette <-
  function(labels, palette = "Tableau 10") {
    label.names <- levels(labels)
    label.palette <-
      palette.colors(length(label.names), palette = palette)
    names(label.palette) <- label.names
    
    label.palette
  }


.generate.label.colours <-
  function(labels, label.palette) {
    label.palette[labels]
  }


.create.edgelist <-
  function(IBD,
           pair1.col = "sample1",
           pair2.col = "sample2",
           ibd.estimate.col = "fract_sites_IBD") {
    IBD[, c(pair1.col, pair2.col, ibd.estimate.col)]
  }


.create.vertices <-
  function(IBD,
           metadata = NULL,
           pair1.col = "sample1",
           pair2.col = "sample2",
           sample.col = "Sample") {
    if (is.null(metadata)) {
      samples <- unique(c(IBD[, pair1.col], IBD[, pair2.col]))
      vertices <- data.frame(samples)
      
    } else {
      metadata <-
        .reorder.metadata(IBD, metadata, sample.col = sample.col)
      vertices <- data.frame(metadata[, sample.col])
    }
    
    names(vertices) <- sample.col
    
    vertices
  }


.plot.IBD <-
  function(IBD,
           IBD.cutoff,
           metadata,
           label.col,
           pair1.col = "sample1",
           pair2.col = "sample2",
           ibd.estimate.col = "fract_sites_IBD",
           label.palette = NULL,
           coords.cache = NULL,
           unlabelled = TRUE) {
    # extract relevant data to create a graph of IBDs
    edgelist <-
      .create.edgelist(
        IBD,
        pair1.col = pair1.col,
        pair2.col = pair2.col,
        ibd.estimate.col = ibd.estimate.col
      )
    metadata <- .reorder.metadata(IBD, metadata)
    vertices <-
      .create.vertices(
        IBD,
        metadata = metadata,
        pair1.col = pair1.col,
        pair2.col = pair2.col,
        ibd.estimate.col = ibd.estimate.col
      )
    
    
    # subset the data to the specified cutoff
    d <- edgelist[edgelist[, ibd.estimate.col] >= IBD.cutoff, ]
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    
    # store metadata about the IBD graph
    if (!is.null(coords.cache)) {
      coords <- coords.cache[[as.character(IBD.cutoff)]]
    } else {
      coords <- layout_(IBD.graph, nicely())
    }
    
    if (is.null(coords)) {
      coords <- layout_(IBD.graph, nicely())
      coords.cache[[as.character(IBD.cutoff)]] <- coords
    }
    
    
    # add label and colour to the nodes
    if (is.null(label.palette)) {
      labels <- .get.labels(metadata, label.col)
      label.palette <- .generate.label.palette(labels)
      
    } else {
      labels <-
        .get.labels(metadata, label.col, label.palette = label.palette)
    }
    
    label.colours <-
      .generate.label.colours(labels, label.palette)
    
    IBD.graph <-
      set_vertex_attr(IBD.graph, "color", value = label.colours)
    
    # plot the IBD graph
    percent.cutoff <- round(IBD.cutoff * 100, digits = 2)
    
    if (unlabelled) {
      plot(
        IBD.graph,
        layout = coords,
        vertex.size = 4,
        vertex.label = NA,
        main = paste0("IBD >=", percent.cutoff, "%")
      )
    } else {
      plot(
        IBD.graph,
        layout = coords,
        vertex.size = 4,
        vertex.label.cex = 0.3,
        main = paste0("IBD >=", percent.cutoff, "%")
      )
    }
    
    
    coords.cache
  }

plot.IBD <-
  function(file,
           IBD,
           IBD.cutoffs,
           metadata,
           label.col,
           pair1.col = "sample1",
           pair2.col = "sample2",
           ibd.estimate.col = "fract_sites_IBD",
           label.palette = NULL,
           onefile = TRUE,
           coords.cache = NULL,
           unlabelled = TRUE,
           plot.legend = TRUE) {
    # extract relevant data to create a graph of IBDs
    edgelist <-
      .create.edgelist(
        IBD,
        pair1.col = pair1.col,
        pair2.col = pair2.col,
        ibd.estimate.col = ibd.estimate.col
      )
    metadata <- .reorder.metadata(IBD, metadata)
    
    # add label and colour to the nodes
    if (is.null(label.palette)) {
      labels <- .get.labels(metadata, label.col)
      label.palette <- .generate.label.palette(labels)
      
    } else {
      labels <-
        .get.labels(metadata, label.col, label.palette = label.palette)
    }
    
    
    if (onefile) {
      sqrt.n.plots <- ceiling(sqrt(length(IBD.cutoffs)))
      bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
      size <- 6 * sqrt.n.plots
      
      pdf(file = file,
          width = size,
          height = size)
      
      par(
        mar = c(1, 11, 1, 1) + 0.1,
        mfrow = c(sqrt.n.plots, sqrt.n.plots),
        cex.main = 1.8
      )
    }
    
    for (i in seq_along(IBD.cutoffs)) {
      IBD.cutoff <- IBD.cutoffs[i]
      
      if (!onefile) {
        newfile <- paste0(file, "_IBD", IBD.cutoff, ".pdf")
        pdf(file = newfile)
        par(mar = c(1, 6, 1, 1) + 0.1)
      }
      
      coords.cache <- .plot.IBD(
        IBD,
        IBD.cutoff,
        metadata,
        label.col,
        label.palette = label.palette,
        coords.cache = coords.cache,
        unlabelled = unlabelled
      )
      
      if (plot.legend) {
        if (!onefile) {
          # add the legend about the IBD graph
          legend(
            "bottomleft",
            legend = levels(labels),
            fill = label.palette,
            cex = 0.8,
            title = label.col,
            inset = c(-0.2, 0),
            xpd = NA
          )
          dev.off()
        } else {
          if (i == bottom.left.index) {
            # add the legend about the IBD graph
            legend(
              "bottomleft",
              legend = levels(labels),
              fill = label.palette,
              cex = 1.5,
              title = label.col,
              inset = c(-0.3, 0),
              xpd = NA
            )
          }
        }
      }
    }
    
    if (onefile)
      dev.off()
    
    invisible(coords.cache)
  }


## custom.palette usage
# custom.palette is a character vector of colour values
# with categories set as the vector names
#
# Example:
# custom.palette <- palette.colors(n = 3, palette = "Tableau 10")
# names(custom.palette) <- c("Control", "7-day", "14-day")


## coords.cache usage
# coords.cache is a list that will store coordinates of the nodes
# depending on the IBD cutoffs used. Used when plotting
# the same data but with different labels
#
# Example:
# IBD.cutoffs <- c(0.125, 0.25, 0.5, 1)
# coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
# names(coords.cache) <- as.character(IBD.cutoffs)
#
# coords.cache <-
#   plot.IBD(IBD.plot.file,
#            IBD,
#            IBD.cutoffs,
#            metadata,
#            label1,
#            coords.cache = coords.cache)
#
# plot.IBD(IBD.plot.file,
#          IBD,
#          IBD.cutoffs,
#          metadata,
#          label2,
#          coords.cache = coords.cache)
