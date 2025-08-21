chainsPlot1 <- function (samplesList,
                         var = NULL
                         ,
                         ind = NULL,
                         burnin = NULL, 
          scale = FALSE, 
          line = NULL, 
          ncols = NULL,
          width = 7,
          height = NULL, 
          legend = !is.null(names(samplesList)), 
          legend.location = "topright", 
          cex = 1,
          traceplot = TRUE, 
          densityplot = TRUE,
          file = NULL,
          titleName=NULL) 
{
  if (!traceplot && !densityplot) 
    return(invisible(NULL))
  if (!(class(samplesList)[1] %in% c("list", "mcmc.list"))) 
    samplesList <- list(samplesList)
  if (!is.null(var)) 
    samplesList <- lapply(samplesList, function(samples) {
      var <- gsub("\\[", "\\\\\\[", gsub("\\]", "\\\\\\]", 
                                         var))
      theseVar <- unlist(lapply(var, function(n) grep(paste0("^", 
                                                             n, "(\\[.+\\])?$"), colnames(samples), value = TRUE)))
      ret <- samples[, theseVar, drop = FALSE]
      if (dim(ret)[2] == 0) 
        stop("variable names misspelled", call. = FALSE)
      ret
    })
  chainParamNamesList <- lapply(samplesList, function(s) colnames(s))
  nChains <- length(samplesList)
  paramNamesAll <- unique(unlist(lapply(samplesList, function(s) colnames(s))))
  nParamsAll <- length(paramNamesAll)
  if (!is.null(line)) 
    if (!is.numeric(line) | length(line) != nParamsAll) 
      stop(paste0("line argument must be numeric vector of length ", 
                  nParamsAll), call. = FALSE)
  if (!is.null(ind) && !is.null(burnin)) 
    stop("only specify either ind or burnin")
  if (!is.null(ind)) 
    samplesList <- lapply(samplesList, function(samples) samples[ind, 
                                                                 , drop = FALSE])
  if (!is.null(burnin)) 
    samplesList <- lapply(samplesList, function(samples) samples[(burnin + 
                                                                    1):nrow(samples), , drop = FALSE])
  if (traceplot + densityplot == 1) {
    if (is.null(ncols)) 
      ncols <- min(nParamsAll, 3)
    nrows <- ceiling(nParamsAll/ncols)
  }  else {
    ncols <- 2
    nrows <- nParamsAll
  }
  if (is.null(height)) 
    height <- if (nrows == 1) 
      3
  else if (nrows == 2) 
    4
  else if (nrows == 3) 
    5
  else if (nrows == 4) 
    6
  else 6.5
  if (!is.null(file)) 
    pdf(file, width = width, height = height)
  else if (inherits(try(eval(parse(text = "knitr::opts_chunk$get('dev')")[[1]]), 
                        silent = TRUE), "try-error") || is.null(eval(parse(text = "knitr::opts_chunk$get('dev')")[[1]]))) 
    #dev.new(width = width, height = height)
  par.save <- par(no.readonly = TRUE)
  oma1 <- if (nrows == 1) 
    0.4
  else if (nrows == 2) 
    0.4
  else 0.5
  mai1 <- if (traceplot & !densityplot) 
    0.1
  else 0.5
  par(mfrow = c(nrows, ncols), mai = c(mai1, 0.4, 0.4, 0.2), 
      oma = c(oma1, 0.4, 0, 0), mgp = c(2, 0.3, 0))
  for (iParam in 1:nParamsAll) {
    thisParamName <- paramNamesAll[iParam]
    
    if(is.null(thisParamName)){
      titleName1 <- paramNamesAll[iParam]
    }else{
      titleName1 <- titleName[iParam]
    }
    if (traceplot) {
      cols <- rainbow(nChains, alpha = 0.75)
      xlim <- c(1, max(unlist(lapply(samplesList, function(s) if (thisParamName %in% 
                                                                  colnames(s)) dim(s)[1] else NULL))))
      ylim <- range(unlist(lapply(samplesList, function(s) if (thisParamName %in% 
                                                               colnames(s)) s[, thisParamName] else NULL)))
      plot(-1000, -1, xlim = xlim, ylim = ylim, xlab = "", 
           ylab = "", main = titleName1, cex.main = cex, 
           cex.axis = 0.8 * cex, tcl = -0.2, xaxt = "n")
      if (iParam == 1 & legend & !is.null(names(samplesList))) 
        legend(legend.location, legend = names(samplesList), 
               lty = 1, col = cols, cex = cex)
      for (iChain in 1:nChains) {
        if (!(thisParamName %in% colnames(samplesList[[iChain]]))) 
          next
        ys <- samplesList[[iChain]][, thisParamName]
        lines(seq_along(ys), ys, col = cols[iChain])
      }
      if (!is.null(line) && !is.na(line[iParam])) 
        segments(x0 = -length(ys), x1 = length(ys), y0 = line[iParam], 
                 col = "black")
    }
    if (densityplot) {
      xMin <- xMax <- yMax <- NULL
      for (iChain in 1:nChains) {
        if (!(thisParamName %in% colnames(samplesList[[iChain]]))) 
          next
        d <- density(samplesList[[iChain]][, thisParamName])
        xMin <- min(xMin, d$x)
        xMax <- max(xMax, d$x)
        yMax <- max(yMax, d$y)
      }
      plot(-1000, -1, xlim = c(xMin, xMax), ylim = c(0, 
                                                     yMax), type = "n", main = titleName1, xlab = "", 
           ylab = "", tcl = -0.2, yaxt = "n", cex.main = cex, 
           cex.axis = 0.8 * cex)
      if (iParam == 1 & legend & !is.null(names(samplesList))) 
        legend(legend.location, legend = names(samplesList), 
               fill = rainbow(nChains, alpha = 0.5), bty = "n", 
               cex = cex)
      for (iChain in 1:nChains) {
        if (!(thisParamName %in% colnames(samplesList[[iChain]]))) 
          next
        ys <- samplesList[[iChain]][, thisParamName]
        polygon(density(ys), col = rainbow(nChains, alpha = 0.2)[iChain], 
                border = rainbow(nChains, alpha = 0.2)[iChain])
      }
      if (!is.null(line) && !is.na(line[iParam])) 
        segments(x0 = line[iParam], y0 = 0, y1 = yMax, 
                 col = "black")
    }
  }
  invisible(par(par.save))
  if (!is.null(file)) 
    dev.off()
}