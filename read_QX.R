read.QX <- function(directory, nr.r = NA){

##‘directory’ is a character vector of length 1 indicating the location of the CSV files
##`nr.r‘ is an optional number indicating the number of rows in the matrix (missing values will be padded with NA)
 #      (if nr.r = NA, the largest number of droplets found will be used)
## The QX .csv files should be the ONLY .csv files in the directory. Otherwise errors will occur!

##Prepping to read
  bup <- getwd()   #save current work directory
  setwd(directory) #go to data directory

##step one: read the files required from the directory
  d.max <- round(20/0.00085, digits = -1) #maximum number of droplet reads possible
  if(isTRUE(nr.r > d.max)){ #in case someone wants more padding than droplets are possible
                 message("this is not a good idea")
                 message(paste("nr.r reduced to ",d.max))
                 nr.r <- d.max
                }
  r.max <- 0 #highest number of droplets found in this run
  files <- Sys.glob("*.csv") #save all filenames
  output.Ch1 <- matrix(, nrow = d.max, ncol = length(files)) #generate output file for channel 1
  output.Ch2 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 2
  for(i in seq_along(files)){        #read each file, add variable zero padding to max size and store in matrix
        temp <- read.csv(files[i], skip = 1, header = F) #read csv file
        n.dr <- dim(temp)[1] #check how many droplets were generated in this run
        #save data into output matrices
        output.Ch1[1:n.dr, i] <- temp[,1]
        output.Ch2[1:n.dr, i] <- temp[,2]
        #check if current number of droplets is higher the run max
        if(n.dr > r.max){
                r.max <- n.dr #we have a new winner
                }
        }

##step two: trimming
  if(isTRUE(nr.r > r.max)){ #in this case we leave padding
   output.Ch1 <- output.Ch1[1:nr.r, ]
   output.Ch2 <- output.Ch2[1:nr.r, ]
  }else{
   output.Ch1 <- output.Ch1[1:r.max, ]
   output.Ch2 <- output.Ch2[1:r.max, ]
   if(isTRUE(nr.r < r.max)){
           message("more droplets than expected")
           message(paste("nr.r increased to ",r.max))}
  }

##step three: output
  setwd(bup) #return to previous wd
  return(list(Ch1 = output.Ch1, Ch2 = output.Ch2))
}
