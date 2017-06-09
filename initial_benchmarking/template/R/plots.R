
# panel.cor puts correlation in upper panels, size proportional to correlation


#' panel.cor
#' 
#' Put correlation in panels of \code{pairs} plot, size proportional to correlation.
#'
#' @param x       \code{num} vector
#' @param y       \code{num} vector
#' @param digits  \code{int} scalar of correlation coefficient digits to show
#' @param prefix  \code{chr} string for correlation coefficient
#' @param cex.cor \code{num} for controling correlation coefficient font size
#'
#' @return 
#' @export
#'
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}



#' ggnice
#' 
#' Addon for ggplots: make theme_classic, optional continuous x or y log2 transformation, and a title.
#'
#' @param title  \code{chr} scalar for plot title
#' @param ld     \code{chr} scalar for log2 transformations ("", "x", "y", "xy")
#'
#' @return 
#' @export
#'
ggnice <- function( title = NULL, ld = c("", "x", "y", "xy") )
{
  ld <- match.arg(ld)
  trafo <- switch(
    ld, 
    x = scale_x_continuous(trans = "log2"),
    y = scale_y_continuous(trans = "log2"),
    xy = c(scale_x_continuous(trans = "log2"), scale_y_continuous(trans = "log2"))
  )
  out <- list(ggtitle(title), theme_classic())
  out <- if(ld == "") out else list(trafo, out) 
  return(out)
}



#' Hist
#' 
#' Convenience histogram plot with option for vlines and log2 trafo.
#'
#' @param vec    \code{num} data vector
#' @param xlab   \code{chr} scalar for x axis title
#' @param vlines \code{num} vector for vertical line positionings, or \code{NULL} for none
#' @param ld     \code{bool} whether to do a log2 transformation
#' @param main   \code{chr} scalar for plot title
#' @param breaks \code{int} scalar for the number of bins
#'
#' @return 
#' @export
#'
Hist <- function( vec, xlab = "x", vlines = NULL, ld = FALSE, main = NULL, breaks = 100 )
{
  df <- data.frame(x = vec)
  trafo <- if(ld) "log2" else "identity"
  gg <- ggplot(df, aes(x = x)) + geom_histogram(bins = breaks) +   
    scale_x_continuous(name = xlab, trans = trafo) +
    theme_classic() + ggtitle(main)
  if(is.null(vlines)) gg else gg + geom_vline(xintercept = vlines, color = "red", linetype = 2)
}



#' Boxplot
#' 
#' Convenience boxplot with options for hlines and log2 trafo
#'
#' @param xx     \code{factor} vector for x axis mapping
#' @param yy     \code{num} vector for y axis
#' @param col    \code{factor} vector for color mapping
#' @param main   \code{chr} scalar for plot title
#' @param xlab   \code{chr} scalar for x axis title
#' @param ylab   \code{chr} scalar for y axis title
#' @param hlines \code{num} vector for horizontal line positionings, or \code{NULL} for none
#' @param ld     \code{bool} whether to do a log2 transformation
#'
#' @return 
#' @export
#'
Boxplot <- function(xx, yy, col, main = NULL, xlab = "Groups", ylab = "metric", hlines = NULL, ld = FALSE)
{
  trafo <- if(ld) "log2" else "identity"
  df <- data.frame(xx =xx, yy = yy, col = col)
  gg <- ggplot(df, aes(x = xx, y = yy, fill = col)) +
    geom_boxplot() + 
    scale_y_continuous(name = ylab, trans = trafo) +
    scale_x_discrete(name = xlab) +
    ggtitle(main) +
    theme_classic()
  if(is.null(hlines)) gg else gg + geom_hline(yintercept = hlines, color = "red", linetype = 2)
}




#' ScatterDensity
#' 
#' Scatter plot with point density mapped to heat color palette.
#'
#' @param x      \code{num} vector with x axis data
#' @param y      \code{num} vector with y axis data
#' @param xl     \code{chr} scalar for x axis title
#' @param yl     \code{chr} scalar for y axis title
#' @param title  \code{chr} title of plot (main)
#' @param ld     \code{chr} optional log2 transformation ("x", "y", "xy")
#'
#' @return 
#' @export
#'
ScatterDensity <- function( x, y, xl = "x", yl = "y", title = "", ld = c("", "x", "y", "xy") )
{
  data <- data.frame(x = x, y = y)
  ld <- match.arg(ld)
  if( ld == "xy" ){
    data$x <- log2(data$x)
    data$y <- log2(data$y)
    xl <- paste("log2", xl)
    yl <- paste("log2", yl)
  } else if( ld == "x" ){
    data$x <- log2(data$x)
    xl <- paste("log2", xl)
  } else if( ld == "y" ){
    data$y <- log2(data$y)
    yl <- paste("log2", yl)
  }
  arr <- densCols(data$y, data$x, colramp= colorRampPalette(c("black", "white")))
  data$dens <- col2rgb(arr)[1, ] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  data$col <- cols[data$dens]
  plot(y ~ x, data = data[order(data$dens), ], 
       pch = 20, col = col, cex = 0.3, bty = "l", xlab = xl, ylab = yl, main = title)
  legend("bottomright", c("high", "low"), 
         fill = c("#FF3100", "#000099"), border = "white", bty = "n", title = "Density")
}


#' ScatterControls
#' 
#' Scatter plot with spike-ins and mitochondiral genes highlighted in red and green.
#'
#' @param x      \code{num} vector with x axis data
#' @param y      \code{num} vector with y axis data
#' @param xl     \code{chr} scalar for x axis title
#' @param yl     \code{chr} scalar for y axis title
#' @param title  \code{chr} title of plot (main)
#' @param spike  \code{bool or int} vector defining points which are spike-ins (\code{NULL} for none)
#' @param mt     \code{bool or int} vector defining points which are mitochondrial (\code{NULL} for none)
#' @param ld     \code{chr} optional log2 transformation ("x", "y", "xy")
#'
#' @return 
#' @export
#'
ScatterControls <- function( x, y, xl = "x", yl = "y", title = "", spike = NULL, mt = NULL, ld = c("", "x", "y", "xy") )
{
  data <- data.frame(x = x, y = y)
  ld <- match.arg(ld)
  if( ld == "xy" ){
    data$x <- log2(data$x)
    data$y <- log2(data$y)
    xl <- paste("log2", xl)
    yl <- paste("log2", yl)
  } else if( ld == "x" ){
    data$x <- log2(data$x)
    xl <- paste("log2", xl)
  } else if( ld == "y" ){
    data$y <- log2(data$y)
    yl <- paste("log2", yl)
  }
  f <- y ~ x
  plot(f, data = data, pch = 20, cex = 1, col = alpha("darkgray", 0.3), bty = "l", xlab = xl, ylab = yl, main = title)
  if(!is.null(spike)) points(f, data = data[spike, ], pch = 20, cex = 0.3, col = "#db4437")
  if(!is.null(mt)) points(f, data = data[mt, ], pch = 20, cex = 1, col = "#0f9d58")
  legend("bottomright", c("Spike-In", "MT"), fill = c("#db4437", "#0f9d58"), border = "white", bty = "n", title = "Control")
}




#' VennDiag
#' 
#' Plots a Venn diagram from co-occuring elements in vectors given in a list.
#'
#' @param l      \code{list} of vectors, each vector is a set
#' @param labs   \code{chr vec} giving set names (default is \code{names(l)})
#' @param cols   \code{chr vec} of colors for each set (will be recycled)
#'
#' @return 
#' @export
#'
VennDiag <- function( l, labs = names(l), cols = c("#4285f4", "#db4437", "#0f9d58", "#f4b400") ){
  cols <- rep(cols, length.out = length(l))
  names(l) <- labs
  grid::grid.newpage()
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grid::grid.draw(VennDiagram::venn.diagram(
    l, filename = NULL, na = "remove", cex = 1.5, lty = 1, lwd = 1, cat.cex = 1,
    fill = cols, 
    main.cex = 1.5, main.fontfamily = 1, cat.fontfamily = 1
  ))  
}



#' PlotNLS
#' 
#' Convenience plotting function which summarizes a nls fit.
#' 4 plots: resid vs fit, real vs fit, quantile real vs quantile fit, ecdf.
#'
#' @param fit    \code{nls} fit object
#' @param real   \code{num} vector giving target values of the fit
#'
#' @return 
#' @export
#'
PlotNLS <- function( fit, real )
{
  opar <- par(mfrow = c(2, 2))
  plot(fitted(fit), resid(fit, type = "p"), pch = 20, ylab = "standardised residuals", xlab = "fitted values")
  abline(h = 0, col = "red")
  plot(fitted(fit), real, pch = 20, xlab = "fitted values", ylab = "real values")
  abline(0, 1, col = "red")
  s <- seq(0, 1, length.out = 100)
  plot(quantile(fitted(fit), s), quantile(real, s), pch = 20, xlab = "quantile fitted values", ylab = "quantile real values")
  abline(0, 1, col = "red")
  plot(ecdf(fitted(fit)), main = "", ylab = "cummulative probability", xlab = "quantile")
  curve(ecdf(real)(x), add = TRUE, col = "red")
  par(opar)
}



#' PlotPredict
#' 
#' Convenience plotting function for plotting data points in grey
#' with the fitted model as red line over them.
#'
#' @param real_x      \code{num} x vector for data points
#' @param real_y      \code{num} y vector for data points
#' @param newdata     \code{num} vector for new data
#' @param prediction  \code{num} vector for predicted data
#' @param xlab        \code{chr} string for x axis label
#' @param ylab        \code{chr} string for y axis label
#' @param main        \code{chr} string for plot title
#'
#' @return 
#' @export
#'
PlotPredict <- function( real_x, real_y, newdata, prediction, xlab = "x", ylab = "y", main = "" )
{
  new <- data.frame(x = newdata, y = prediction)
  real <- data.frame(x = real_x, y = real_y)
  ggplot(new) + geom_point(aes(x, y), data = real, color = "gray") +
    geom_line(aes(x, y), color = "red") + 
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab) +
    ggnice(main)
}



#' Dens
#' 
#' Convenience plotting functions plots density estimations with a color mapping.
#'
#' @param vec   \code{num} vector of data points
#' @param col   \code{factor} vector for mapping of colors
#' @param xlab  \code{chr} string for x axis title
#' @param leg   \code{chr} string for legend title
#' @param main  \code{chr} string for plot title
#' @param ld    \code{bool} whether on log2 scale
#'
#' @return 
#' @export
#'
Dens <- function(vec, col, xlab = "x", leg = "color", main = "", ld = FALSE)
{
  df <- data.frame(
    vec = vec,
    col = col
  )
  trafo <- if(ld) "log2" else "identity"
  ggplot(df, aes(x = vec, fill = col, color = col)) + geom_density(alpha = 0.2) + 
    scale_color_discrete(leg) +
    scale_fill_discrete(leg) +
    scale_x_continuous(xlab, trans = trafo) +
    ggnice(main) 
}



#' PlotDistro
#' 
#' Convenience plotting function plots histogram in grey with lines describing density curves on top.
#' If \emph{prediction} is \code{NULL} then density from \emph{newdata} is estimated.
#'
#' @param real_x      \code{num} vector of x data points to create histogram from
#' @param newdata     \code{num} vector of new data with x axis values
#' @param prediction  \code{num} vector of density prediction for newdata or \code{NULL}
#' @param col         \code{factor} vector for color mapping
#' @param xl          \code{chr} string for x axis title
#' @param leg         \code{chr} string for legend title
#' @param main        \code{chr} string for plot title
#' @param ld          \code{bool} whether on log2 scale
#'
#' @return 
#' @export
#'
PlotDistro <- function(real_x, newdata, prediction = NULL, col = "1", xl = "x", leg = "color", main = "", ld = FALSE)
{
  df <- data.frame(
    vec = newdata,
    col = col
  )
  if( !is.null(prediction) ) df$pred <- prediction
  model <- if( is.null(prediction) ){
    geom_density(aes(x = vec, col = col, fill = col), data = df, alpha = 0.2)
  } else { 
    geom_line(aes(x = vec, y = pred, color = col), data = df)
  }
  trafo <- if(ld) "log2" else "identity"
  data <- data.frame(x = real_x)
  ggplot(data) + geom_histogram(aes(x = x, y = ..density..), bins = 100, fill = "gray") +
    model +
    scale_color_discrete(leg) +
    scale_fill_discrete(leg) +
    scale_x_continuous(xl, trans = trafo) +
    ggnice(main)
}



#' PlotPredictMult
#' 
#' Convenience plotting function for plotting data points in grey
#' with multiple fitted models as colored lines over them.
#'
#' @param real_x      \code{num} x vector for data points
#' @param real_y      \code{num} y vector for data points
#' @param newdata     \code{num} vector for new data
#' @param prediction  \code{num} vector for predicted data
#' @param col         \code{factor} vector for mapping color
#' @param xlab        \code{chr} string for x axis label
#' @param ylab        \code{chr} string for y axis label
#' @param main        \code{chr} string for plot title
#'
#' @return 
#' @export
#'
PlotPredictMult <- function( real_x, real_y, newdata, prediction, col, xlab = "x", ylab = "y", leg = "", main = "" )
{
  new <- data.frame(x = newdata, y = prediction, col = col)
  ylim <- range(real_y)
  real <- data.frame(x = real_x, y = real_y)
  ggplot(new) + geom_point(aes(x, y), data = real, color = "gray") +
    geom_line(aes(x, y, color = col)) + 
    scale_color_discrete(leg) +
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab, limits = ylim) +
    ggnice(main)
}




#' Barplot
#' 
#' Convenience plotting function for bar plot with color mapping.
#'
#' @param fac         \code{factor} vector for data points
#' @param col         \code{factor} vector for mapping color
#' @param xlab        \code{chr} string for x axis label
#' @param lg          \code{chr} string for legend title
#' @param main        \code{chr} string for plot title
#' @param ld          \code{bool} whether on log2 scale
#'
#' @return 
#' @export
#'
Barplot <- function(fac, col, xl = "x", lg = "col", main = "", ld = FALSE)
{
  df <- data.frame(fac = fac, col = col)
  trafo <- if(ld) "log2" else "identity"
  ggplot(df, aes(x = fac, fill = col)) + geom_bar(position = "dodge") +
    scale_fill_discrete(lg) +
    scale_x_discrete(xl) +
    scale_y_continuous(trans = trafo) + 
    ggnice(main)
}
