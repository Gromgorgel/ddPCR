## Function to read the raw fluorescence output files from the Stilla Naica
## digital PCR platform. Files can be exported by choosing options>export amplitude and cluster data
## The program will create a folder (droplet-level-data/RawData) with several .csv files per reaction (well).
## To read data from multiple runs, all files should be placed into one folder. The latter's path (eg. "D:/Documents/dPCR Run")
## is submitted as a string to the function below.
## The ddPCR .csv files should be the ONLY .csv files in the directory.
## Otherwise errors may occur!

 # Updates will appear at: github.com/Gromgorgel

## Input:
 #‘directory' is a character vector of length 1 indicating the location of the CSV files, the algorithm will attempt to analyze
 #            any csv files present in the directory (unless 'file' is specified) and expects all csv files to contain droplet amplitudes.
 #            defaults to "C:\\Quantalife\\Data"
 #‘files' is a character vector. if specified, the algorithm will analyse only the filenames in the vector and will look for them
 #       in the current working directory. The file extension should be present in the file name.
 #‘nr.r'     is an optional number indicating the number of rows in the matrix (missing values will be padded with NA)
 #            (if nr.r = NA, the largest number of droplets found will be used)
 # autoname   logical. If TRUE the well names are collected from the csv files
 #            and used as column names. Only works if there are 96 or less csv files
 #            If there are more, a warning message will be reported. Problems
 #            may occur if additional underscores ("_") are present in the name of the run/experiment.

## Output:
 # The function creates a list containing matrices with fluorescence readings
 # one for each channel: $Ch1, $Ch2, etc.. The matrix has as many columns
 # as there are samples. The number of rows can be either set (max: approx 23500)
 # otherwise, the number of rows equals the largest number of readings found.
 # Missing values are reported as NA.

 ## TO DO
 ## select chip type
 ## check data structure with less channels

read.Naica <- function(directory = "C:\\droplet-level-data\\RawData", files = NA, nr.r = NA, autoname = TRUE){
##Prepping to read
  bup <- getwd()   #save current work directory
  if(is.na(files)){  # files not specified
    setwd(directory) # go to data directory
    files <- Sys.glob("*.csv") #save all filenames
    files <- files[!seq_along(files) %in% grep('Channel', files)] # select files with all channel data (these do not contain 'Channel' in the filename)
  }else{             # files is specified, we look in the current working directory
    directory <- bup
  }

 # check rownumber required  #note 0.25 nl is an estimate I made, it is not the actual droplet volume, I will add this once I find it.
  d.max <- round(5/0.00025, digits = -1)*1.275 #maximum number of droplet reads possible, add 27.5% extra as some reactions have been yielding more droplest than possible
  if(isTRUE(nr.r > d.max)){ #in case someone wants more padding than droplets are possible
                 message("this is not a good idea")
                 message(paste("nr.r reduced to ",d.max))
                 nr.r <- d.max
                }
  r.max <- 0 #highest number of droplets found in this run

  ##step one: prep the list to store all data
  list.out <- list()
  list.out$Ch1 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 1
  list.out$Ch2 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 2
  list.out$Ch3 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 3
  list.out$Ch4 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 4
  list.out$Ch5 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 5
  list.out$Ch6 <- matrix(nrow = d.max, ncol = length(files)) #generate output file for channel 6

  for(i in seq_along(files)){        #read each file, add variable zero padding to max size and store in matrix
        temp <- read.csv(files[i], skip = 1, header = F) #read csv file
        n.dr <- dim(temp)[1]  # check how many droplets were generated in this sample
        n.ch <- dim(temp)[2] - 3 # check number of channels there are in this sample (there are three non-channel columns: x & y coordinates, & index)
        #save data into output matrices
    for(j in 1:n.ch){
        list.out[[j]][1:n.dr, i] <- temp[,j + 2]
        #check if current number of droplets is higher the run max
        if(n.dr > r.max){
                r.max <- n.dr #we have a new winner
        }}} # close all statments

##step two (optional): extract colnames from files
  if(isTRUE(autoname)){
    if(length(files) > 96){
          message("more than 96 files found, switching off autoname")
    }else{
      wells <-vector(mode="character", length(files))
      for(i in seq_along(files)){  #extract colnames from each file
          wells[i] <- sub("?.*_(.*?)_RawData.csv","\\1",files[i])
      }
     #check for errors and report:
      if(length(wells) > length(files)){
            message("more names than columns found, switching off autoname")
       }else{
        for(j in 1:n.ch){  # I looked for an alternative to a loop here, but I could not get it to work with matrices.
         colnames(list.out[[j]]) <- wells
    }}}}# close all statements

##step three: trimming
  r.trim <- if(isTRUE(nr.r > r.max)){nr.r}else{r.max} # if an exact number of rows was specificied (and it is larger than r.max, we keep it)
  for(j in 1:n.ch){
   list.out[[j]] <- list.out[[j]][1:r.trim, ]
  }
##step four: output
  setwd(bup) #return to previous wd
  return(list.out)
}