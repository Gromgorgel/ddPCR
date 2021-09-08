## Function to read the raw fluorescence output files from the Biorad QX200
## digital PCR platform. Files can be exported by choosing options>export amplitude and cluster data
## The program will create one .csv file per reaction (well). All files
## should be placed into one folder. The latter's path (eg. "D:/Documents/dPCR Run")
## is submitted as a string to the function below.
## The QX .csv files should be the ONLY .csv files in the directory.
## Otherwise errors will occur!

 # Updates will appear at: github.com/Gromgorgel

## Input:
 #‘directory’ is a character vector of length 1 indicating the location of the CSV files
 #`nr.r‘      is an optional number indicating the number of rows in the matrix (missing values will be padded with NA)
 #            (if nr.r = NA, the largest number of droplets found will be used)
 # autoname   logical. If TRUE the well names are collected from the csv files
 #            and used as column names. Only works if there are 96 or less csv files
 #            If there are more, a warning message will be reported. Problems
 #            may occur if additional underscores ("_") are present in the name of the run/experiment.

## Output:
 # The function creates a list containing two matrices with fluorescence
 # readings (one for each channel: $Ch1 and $Ch2). The matrix has as many columns
 # as there are samples. The number of rows can be either set (max: approx 23500)
 # otherwise, the number of rows equals the largest number of readings found.
 # Missing values are reported as NA.

read.QX <- function(directory, nr.r = NA, autoname = TRUE){
##Prepping to read
  bup <- getwd()   #save current work directory
  setwd(directory) #go to data directory

##step one: read the files required from the directory
  d.max <- round(20/0.00085, digits = -1) *1.275 #maximum number of droplet reads possible, add 27.5% extra as some reactions have been yielding more droplest than possible
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

##step two (optional): extract colnames from files
  if(isTRUE(autoname)){
    wells <- vector()
    for(i in seq_along(files)){  #extract colnames from each file
        wells <- append(wells, sub("?.*_(.*?)_.*","\\1",files[i]))
    }
   #check for errors and report:
    if(length(wells) > 96 | length(wells) > length(files)){
      if(length(wells) > 96){
          message("more than 96 names found, switching off autoname")
      }else{
          message("more names than columns found, switching off autoname")
      }
    }else{  # add names to matrices
      colnames(output.Ch1) <- wells
      colnames(output.Ch2) <- wells
    }
  }
##step three: trimming
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

##step four: output
  setwd(bup) #return to previous wd
  return(list(Ch1 = output.Ch1, Ch2 = output.Ch2))
}
