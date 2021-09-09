ddplot <- function(drp, ...){
      drp <- na.omit(drp)
      temp <- sample(c(1:length(drp)), size = length(drp), replace = F)
      plot(drp ~ temp, col = rgb(123, 168, 194, maxColorValue = 255), bty = "L", pch = 20 , xlab = "observation", ylab = "FU", ylim = c(min(drp), max(drp)), ...)}
