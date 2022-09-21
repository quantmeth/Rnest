#' Print results of NEST
#'
#' @description Scree plot of the eigenvalues and the \code{(1-alpha)*100\%} confidence intervals derived from the resampled eigenvalues supplied to \code{nest}.
#' 
#' @param x An object of class "nest".
#' @param y Further arguments for other methods, ignored for "nest".
#' @param ... Further arguments for other methods, ignored for "nest".
#' 
#' @note This function is more interesting with many \code{alpha} values.
#'
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom grDevices grey
#' @export
#'
#' @examples
#' results <- nest(ex_2factors, n = 100, alpha = c(.01, .05, .01))
#' plot(results)
#' # Return the data used to produce the plot
#' df <- plot(results)$data
plot.nest <- function(x, y, ...){
  
  df <- data2plot(x)
  
  ggplot2::ggplot(df,
                  mapping = aes_string(x = "Position",
                                       y = "Eigenvalues",
                                       color = "Confidence")) +
    geom_line(linetype = "dashed") +
    geom_point() +
    scale_color_manual(values = c(grey(seq(0.75, .2, length.out = length(x$alpha))), "blue")) +
    scale_x_continuous(breaks = scales::pretty_breaks())+
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    labs(title = "Next Eigenvalues Sequential Test") +
    theme(legend.position = c(.8, .8))
  
  
}

data2plot <- function(x){
  df <- data.frame(Position = 1:length(c(x$values)),
                   Eigenvalues = c(x$values),
                   Confidence = c("Original")
  )
  
  for(i in 1:length(x$Eig)){
    df <- rbind(df, data.frame(Position = rep(i, length(x$alpha)),
                               Eigenvalues = x$Eig[[i]][,i],
                               Confidence = paste0((1 - sort(x$alpha)) * 100,"%")))
  }
  
  rownames(df) <- NULL
  
  return(df)
}

