# Digital Droplet PCR analysis
This branch contains two R scripts related to digital PCR analysis. Both scripts were written to handle output from the Biorad QX200 platform, but the main analysis algortighm (cloudy.R) should be able to handle any type of digital PCR output provided the data have the correct format (more on that later).
The two scripts are:
- cloudy.R
- read_QX.R

Both algorithm were originally provided as supplemental material to the publication "*Measuring digital PCR quality: Performance Parameters and their Optimization*" by Antoon Lievens, Sara Jacchia, Dafni Kagkli, Cristian Savini, and Maddalena Querci (*link to Plos One article goes here once publication is finalized*). In addition, this branch also contains the full dataset on which the publication was based (`Dataset.zip`) which may be used to further explore the functions of the scripts provided.  

The following paragraphs contain a brief overview of the functionality of the algorithms. Details on their use is also contained as commented text in the .R files themselves. The last paragraph contains information on the structure of the dataset provided (i.e. `Dataset.zip`).

## The 'read.QX' algorithm
This function reads the fluorescence values, as exported from the Biorad Quantasoft software, into R. To export the fluorescence readings, choose `options > export amplitude and cluster data` in the `setup` tab of the Quantasoft software (V1.6.6.0320). The software will prompt you to choose folder, for use with read_QX.R it is recommended that you choose an empty folder. For each well in the analysis, a .csv file will be created that contains the fluorescence reads for each channel as well as the assigned cluster.

The function output is a list containing two matrices (one for each channel). Each column in the matrix is a reaction and each row contains fluorscence readings. As the number of readings per reaction may differ, reactions with less accepted droplets are padded with `NA` to fill the matrices. The order of the reactions in the matrix is the order in which they were read from the folder (usually alphabetical) it is therefore recommended to set `autoname = TRUE` in order to be able to identify each column correctly. The cluster data is currently not imported into R as the cloudy algortihm has no use for it. Note that the algorithm will attempt to read **any** .csv file present in a folder, so make sure that only .csv files containing raw ddPCR data are present.

The basic function call to read.QX takes the following arguments:
```
read.QX(directory, nr.r = NA, autoname = TRUE)
```
- `directory` is a character vector of length 1 indicating the location of the CSV files (either path or name of the folder in the working directory)
- `nr.r` is an optional numerical of length 1 indicating the number of rows to be used in the matrix (missing values will be padded with NA). If `nr.r = NA`, the largest number of droplets in the analysis will be used. This option is usefull if you want to read several folders of data and then merge the resulting matrices, for instance, using `cbind`.
- `autoname` is a logical, if `TRUE` the well names are collected from the .csv files and used as column names. Essentially, the script looks for whatever is contained between two underscores ("_") and uses this as the column name (The biorad software adds a suffix to the plate name containing the well name).  Therefore, it is strongly recommended **never** to use underscores in the plate name, this may crash the algortihm. If you did and still want to use the algorithm, you can either change to code to fit your plate names, or look into a batch-rename program (e.g. [rename master](http://rename-master.en.softonic.com/)). For similar reasons the algorithm currently only works if there are 96 or less csv files  (to prevent duplicate column names). If there are more, a warning message will be reported. 

### Example
Unzip the folder contained in `example.zip` into you working directory and have all .R and files from this branch placed there as well before running the example. 

```
source("read_QX.R")

# read files from folder
d.data <- read.QX(directory = "example QXfiles")

# plot a single reaction
plot(sample(d.data$Ch1[, 1]), ylab = "Amplitude", bty = "L")
```

## The 'cloudy' algorithm

This is the algorithm for digital PCR analysis from raw compartment fluorescence measurements

The original algorithm is supplemental material to the publication "Measuring digital PCR quality: Performance Parameters and their Optimization" by Antoon Lievens, Sara Jacchia, Dafni Kagkli and Cristian Savini, and Maddalena Querci (*link to Plos One article goes here once publication is finalized*). For details on the interpretation of results and background on the algorithm we refer the user to the original text of the publication.

The standard function output is a list with the following components:
- `targets.in.sVol`, a numerical vector of length three containing the estimated number of targets in the sample volume (`targ.in.sVol`) and its Poisson confidence bounds (`upper`, `lower`).
- `lambda`, a numerical vector of length three containing the estimated average number of targets per compartment (`lambda`) and its Poisson confidence bounds (`upper`, `lower`).
- `performance`, a numerical vector of length three containing the performance paramters calculated from the input data: `resolution`, ratio of rain to total compartments (`%rain`), and the degree of compartmentalization (`%comparted`)
- `droplets`, a numerical vector of length three containing the number of compartments counted in each category (`positive`, `negative`, & `rain`).
- `populations`, a numerical of length 1, the number of fluorescence populations as detected by the algorithm 

Note that if `vec` is set to `TRUE`, all of the above parameters are returned in a single vector rather than as a list
 
The basic function call to cloudy takes the following arguments:
```
cloudy(drp, dVol = 0.85, sVol = 20, plots = FALSE, silent = TRUE, vec = FALSE)
```
- `drp` is a numeric vector of all (endpoint) fluorescence measurments in a digital reaction. Readings do **not** have to be ordored in any particular way although Quantasoft export is usually sorted from small to large. `NA` values are allowed (will be removed). Negative values are allowed as well (baseline subtraction may cause these in the Quantasoft export).
- `dVol` is a numerical of length 1, the compartment (droplet) volume in nanoliter (standard = 0.85)
- `sVol` is a  numerical of length 1, the sample volume in microliter (standard = 20) 
- `plots` is a logical, if set to `TRUE` plots will be generated (see below)
- `silent` is a logical, if is set to `FALSE` warning messages will be generated detailing the choices the algorithm makes when deviating from the standard analysis routine. 
- `vec` is a logical, if set to `TRUE` the results will be returned in a vector instead of a list. This is useful when you batch analyse dPCR reactions using `apply` and want the results to be returned as a matrix.

When `plots = TRUE` the algortihm will produce two graphical windows with three visual representations of the analysis and its results (two in the first window, one in the second):
- The first window's **plot A** contains a kernel density plot of the fluorescence values, population boundaries and peaks are indicated with vertical dashed lines. 
- The first window's **plot B** contains a dot plot of the fluorescence readings (randomized order), each popultion is coloured differently, a horizontal line represents the threshold that was set by the algorithm
- The second window contains a barplot with the relative values of the performance paramters. Four parameter values are shown (and scaled to be comparable). 
  * `Resolution` (devided by 10), the horizontal line represents the 2.5 (0.25) acceptance limit
  * `Percentage rain` (multiplied by 10), the horizontal line represents the 2.5% (0.25) rejection limit
  * `percentage of sample compartmentalized` (horizontal lines represent the 0.3 and 0.5 limits)
  * `droplet count`: positive and negative droplets (stacked) scaled to the total number of droplets in the analysis

### Example
Please have all .R and .RData files from this branch in your working directory before running the example. Installation of the package `SuppDists` may be required before the algorithm functions.

```
load("dpcr-Example.RData")
source("Cloudy-V2-04.R")

# Analyse a single reaction
cloudy(dd.pcr$Ch1[, 1])

# Analyse a single reaction with plots
cloudy(dd.pcr$Ch1[, 1], plots = T)

# Example of a reaction with multiple populations
cloudy(dd.pcr$Ch1[, 25], plots = T)

# Example of batch analysis
results <- apply(dd.pcr$Ch1, 2, cloudy, vec = T)
head(results[, 1:2])
```

### Using cloudy for non-Biorad data
Even though the algorithm was design for use with data produced by the Biorad QX platform, there is no immediate reason the procedure cannot be applied to dPCR data from another source. The algorithm expects a vector as input in which each value is the endpoint fluorescence reading from a single chamber/droplet. As long as the input meets these expectations, the output should make sense. 

## 'Dataset.zip': contents and strucuture
The dataset is provided as a zipped `.csv` file. The is essentially a matrix of 16996 rows by 473 columns. The first column provides the rownames, each subsequent column is a digital PCR reaction. The rows contain the following information:
  * `react-ID` : 1 to 472. Numerical identifier of each reaction.
  * `plate-ID`: 1 to 9. Numerical identifier of each plate. A total of 9 plates was run, representing the different experiments in the dataset. **Plate 1** tests all PCR targets in 4 repeats using qPCR validated conditions. **Plate 2** is the primer/probe concentration gradient. **Plate 3** is the "rain" dilution series. **Plate 4** tests the PCR enhancers. **Plate 5** is the cycle gradient. **Plate 6** is the sonication gradient. **Plate 7** is the annealing temperature gradient. **Plate 8** tests digital touchdown PCR. **Plate 9** is the final run of reactions under optimized conditions.
  * `target` : name of the PCR target
  * `Primers` : Primer concentration in nM
  * `Probe` : Probe concentration in nM
  * `Anneal` : Annealing temperature used during the run.
  * `Touchdown` : Was digital touchdown PCR used? Either `Y` or `N`.
  * `Enhancer` : PCR enhancer that was added to the reaction. `NA` if no enhancer was used.
  * `Cycles` : Number of cycles run.
  * `Sonication` : Amount of sonication (in seconds). `NA` if the sample was not sonicated. 
  * `FU0001` to `FU16985` : The droplet fluorescence values. In case there are less than 16985 droplets in a reaction, `NA` is used to padd the column.
  
 If these descriptions do not seem meaningful, you might want to read our publication (*link goes here one more time*) for more background on the experiments.
