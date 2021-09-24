ddplot <- function(drp, well = 1, channel = 1,...){
      drp1 <- if(channel == 1){na.omit(drp$Ch1[, well])}else{na.omit(drp$Ch2[, well])}
      temp <- sample(c(1:length(drp1)), size = length(drp1), replace = F)
      plot(drp1 ~ temp, col = rgb(123, 168, 194, maxColorValue = 255), bty = "L", pch = '.' , cex = 5, xlab = "observation", ylab = "FU", ...)}

tdplot <- function(drp, well = 1,...){
      drp1 <- na.omit(drp$Ch1[, well])
      drp2 <- na.omit(drp$Ch2[, well])
      plot(drp1 ~ drp2, col = rgb(123, 168, 194, maxColorValue = 255), bty = "L", pch = '.' , cex = 5, xlab = "FU_Ch2", ylab = "FU_Ch1", ...)}
