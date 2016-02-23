# Digital Droplet PCR analysis
This branch contains two R scripts related to digital PCR analysis. Both scripts were written to handle output from the Biorad QX200 platform, but the main analysis algortighm (cloudy.R) should be able to handle any type of digital PCR output provided the data have the correct format (more on that later).
The two scripts are:
- cloudy.R
- read_QX.R

The following paragraphs contain a brief overview of the functionality of both algorithms. Details on their use is also contained as commented text in the .R files themselves.

## the 'read.QX' algorithm
This function reads the flourescence values, as exported from the Biorad Quantasoft software, into R. To export the fluorescence readings, choose `options > export amplitude and cluster data` in the `setup` tab of the Quantasoft software (V1.6.6.0320). The software will prompt you to choose folder, for use with read_QX.R it is recommended that you choose an empty folder. For each well in the analysis, a .csv file will be created that contains the fluorescence reads for each channel as well as the assigned cluster.

The function output is a list containing two matrices (one for each channel). Each column in the matrix is a reaction and each row contains fluorscence readings. As the number of readings per reaction may differ, reactions with less accepted droplets are padded with `NA` to fill the matrices. The order of the reactions in the matrix is the order in which they were read from the folder (usually alphabetical) it is therefore recommended to set `autoname = TRUE` in order to be able to identify each column correctly. The cluster data is currently not imported into R as the cloudy algortihm has no use for it.

The basic function call to read.QX takes the following arguments:
```
read.QX(directory, nr.r = NA, autoname = TRUE)
```
- `directory` is a character vector of length 1 indicating the location of the CSV files (either path or name of the folder in the working directory)
- `nr.r` is an optional numerical of length 1 indicating the number of rows to be used in the matrix (missing values will be padded with NA). If `nr.r = NA`, the largest number of droplets in the analysis will be used. This option is usefull if you want to read several folders of data and then merge the resulting matrices using, for instance, `cbind`.
- `autoname` is a logical, if `TRUE` the well names are collected from the .csv files and used as column names. Essentially, the script looks for whatever is contained between two underscores ("_") and uses this as the column name (The biorad software adds a suffix to the plate name contianing the the well).  Therefore, it is strongly recommended **never** to use underscores in the plate name, this may crash the algortihm. If you did and still want to use the algorithm, you can either change to code to fit your plate names, or look into a batch-rename program (e.g. (rename master)[http://rename-master.en.softonic.com/]). For similar reasons the algorithm currently only works if there are 96 or less csv files  (to prevent duplicate column names). If there are more, a warning message will be reported. 

## The 'cloudy' algorithm
bleh bleh

Ni Ni
