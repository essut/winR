#!/usr/bin/env Rscript

.add.label.to.IBD <-
  function(IBD, metadata, label.col, by.col = "Sample") {
    metadata <- metadata[, c(by.col, label.col)]
    
    IBD <- merge(
      IBD,
      metadata,
      by.x = paste0(by.col, 1),
      by.y = by.col,
      all.x = TRUE,
      sort = FALSE
    )
    colnames(IBD)[length(IBD)] <- paste0(label.col, 1)
    
    IBD <- merge(
      IBD,
      metadata,
      by.x = paste0(by.col, 2),
      by.y = by.col,
      all.x = TRUE,
      sort = FALSE
    )
    colnames(IBD)[length(IBD)] <- paste0(label.col, 2)
    
    IBD
  }


plot.IBD.distribution <-
  function(y,
           x,
           data,
           filename,
           col = NULL,
           width = 2160,
           height = 2160,
           res = 300) {
    formula <- as.formula(paste(y, "~", x))
    
    if (is.null(col))
      col <- "lightgray"
    
    png(
      filename = boxplot.file,
      width = width,
      height = height,
      res = res,
      type = "cairo"
    )
    
    # initialise plot
    boxplot(
      formula,
      data = data,
      xlab = x,
      ylab = y,
      col = NA
    )
    
    # plot jitters
    stripchart(
      formula,
      data = data,
      method = "jitter",
      jitter = 0.4,
      vertical = TRUE,
      pch = 16,
      col = col,
      add = TRUE,
      axes = FALSE
    )
    
    # add box on top
    boxplot(formula,
            data = data,
            col = NA,
            add = TRUE)
    
    dev.off()
  }


synchronise.data <-
  function(IBD,
           metadata,
           pair1.col = "sample1",
           pair2.col = "sample2",
           sample.col = "Sample") {
    IBD[IBD[, pair1.col] %in% metadata[, sample.col] &
          IBD[, pair2.col] %in% metadata[, sample.col], ]
  }

synchronise.metadata <-
  function(metadata,
           IBD,
           pair1.col = "sample1",
           pair2.col = "sample2",
           sample.col = "Sample") {
    metadata[metadata[, sample.col] %in% c(IBD[, pair1.col], IBD[, pair2.col]), ]
  }
