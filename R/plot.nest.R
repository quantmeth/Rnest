#' Plot results of Next Eigenvalues Sufficiency Test (NEST)
#'
#' @description Scree plot of the eigenvalues and the \code{(1-alpha)*100\%} confidence intervals derived from the resampled eigenvalues supplied to \code{nest}.
#' 
#' @param x an object of class "nest".
#' @param pa show results of Parallel Analysis.
#' @param ... further arguments for other methods, ignored for "nest".
#' 
#' @return A ggplot output.
#' 
#' @note This function is more interesting with many \code{alpha} values.
#'
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom grDevices grey 
#' @importFrom grDevices rgb
#' @export
#'
#' @examples
#' results <- nest(ex_2factors, n = 100, alpha = c(.01, .05, .01))
#' plot(results)
#' # Return the data used to produce the plot
#' df <- plot(results)$data
plot.nest <- function(x, pa = FALSE, ...){
  
  df <- data2plot(x, pa)
  # Position <- df$Position
  # Eigenvalues <- df$Eigenvalues
  # Confidence <- df$Confidence 
  
  if(pa == FALSE){
    col.pal <- c(grey(seq(0.75, .2, length.out = length(x$alpha))), "blue")
  } else {
    col.pal <- c(grey(seq(0.75, .2, length.out = length(x$alpha))), "blue", c(rgb(1, 
                                                                                  seq(0, .5, length.out = length(x$alpha)), 
                                                                                  seq(0, .5, length.out = length(x$alpha))))[length(x$alpha):1])
  }

  out <- ggplot2::ggplot(data = df,
                       mapping = aes(x = .data$Position,
                                     y = .data$Eigenvalues,
                                     color = .data$Confidence)) +
    geom_line(linetype = "dashed") +
    geom_point() +
    scale_color_manual(values = col.pal) +
    scale_x_continuous(breaks = scales::pretty_breaks())+
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    labs(title = paste(x$stopping.rule)) +
    theme(legend.position = c(.8, .8))
  
  out$data <- df
  
  out
}



data2plot <- function(x, pa = FALSE){
  
  nv <- length(x$values) 
  
  df <- data.frame(Position = 1:nv,
                   Eigenvalues = c(x$values),
                   Confidence = "Original")
  
  if(length(x$Eig) == 1){
    
    pa <- TRUE
    
  } else {
    
    for(i in 1:length(x$Eig)){
      df <- rbind(df, data.frame(Position = rep(i, length(x$alpha)),
                                 Eigenvalues = x$Eig[[i]][,i],
                                 Confidence = paste0((1 - sort(x$alpha)) * 100,"%")))
    }
  }
  
  if(pa == TRUE){
    df <- rbind(df, data.frame(Position = rep(1:nv, length(x$alpha)),
                               Eigenvalues = c(t(x$Eig[[1]])),
                               Confidence = rep(paste0("PA ", (1 - sort(x$alpha)) * 100,"%"), each = nv)))
  }
  
  rownames(df) <- NULL

  return(df)
}

