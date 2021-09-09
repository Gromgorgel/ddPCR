#########################################################################################################
## This is the algorithm for digital PCR analysis from raw compartment fluorescence measurements
## This algorithm is supplemental material to the publication "Measuring digital PCR quality: Performance Parameters and their Optimization"
## by Antoon Lievens, Sara Jacchia, Dafni Kagkli and Cristian Savini, and Maddalena Querci

 # Updates will appear at: github.com/Gromgorgel

## ARGUMENTS ###############
## 'drp' = named list, the input should contain 2 values ('Ch1' and 'Ch2') each of which should be a matrix
 #         containing the (endpoint) fluorescence measurments (one column per reaction). The readings do not
 #         have to be ordored in any particular way, NA are allowed (will be removed)
## 'well' = vector (numerical or character), selects which reactions are analyzed. Either by column number (eg c(1,2,4,5)) or
 #          by column name (eg c("A01", "C02", "H12")) in the 'drp' input. defaults to 'all' (analyzes all reactions in the dataset)
## 'method' = string, type of analysis to be made. One of the following: "simplex", "simplex2", "eva", or "duplex"
 #       'simplex' = standard 'cloudy' analysis of the flourescence in channel 1
 #       'simplex2' = standard 'cloudy' analysis of the flourescence in channel 2
 #       'eva' = simplex analysis of channel 1 with channel crosstalk correction using channel 2 (designed for Evagreen)
 #               this generally improves the resolution of the reaction but is noticeably slower.

## 'dVol' = numerical, the compartment volume in nanoliter
## 'sVol' =  numerical, the sample volume in microliter
## 'plots' = logical, if set to 'TRUE' plots will be generated
## 'silent' = logical, if is set to 'FALSE' messages will be generated
## 'vec' = logical, if set to 'TRUE' the results will be returned in a vector instead of a list
## 'threshold' = numerical, OPTIONAL: manually set threshold for counting positive and negative droplets (in Fluorescence Units).
 #               The algorithm will still try to define the populations in order to calculate the statistics.
 #               Note that even when wrong or invalid the manual threshold will not be overruled
## 'neg.ref' = is a numerical of length 1. OPTIONAL: this provides a reference where the algorithm
 #              expects the negative popultation to be. This will only be used when the algorithm can only find a single
 #              population: it wil serve as a reference to decide it the droplets are all negative or all positive.

## OUTPUT ###############
## standard output is a list with the following components
 # targets.in.sVol = the estimated number of targets in the sample volume (targ.in.sVol) and its Poisson confidence bounds (upper, lower)
 # lambda  = the estimated average number of targets per compartment (lambda) and its Poisson confidence bounds (upper, lower)
 # performance = Performance paramters calculated from the input data: resolution, ratio fo rain to total compartments (%rain), degree of compartmentalization (%comparted)
 # droplets = number of compartments counted in each category (positive, negative, & rain)
 # populations = number of fluorescence populations as detected by the algorithm
 # threshold = Fluorescence value used as threshold
## if 'vec' is set to true, all the above parameters are returned in a single vector rather than as a list

## PLOTS ###############
## PLOT 1: Kernel density plot, population boundaries and peak are indicated with vertical dashed lines and are shaded
## PLOT 2: dot plot of the flourescence readings (random order), each popultion is coloured differently, horizontal line represents the threshold set

## Version history: this is stand alone version 3.04
                  # V02
                    #01 added a step to determine the number of populations (eg: only one => abort, two => OK, three => warn and continue)
                    #02 further fine-tuned to 'only one population' procedure for the cases where there are very few pos or neg droplets
                    #03 restructuring of the algorithm for stability and efficiency reasons
                    #04 added control step for overlapping bounderies & use of minima as a general backup for non-sd based boundaries
                    #05 added further robustness when dealing with overlapping populations, added manual threshold option, added threshold to output
                    #06 improvements for stability, removal of some typos, fixed a few bugs related to manual threshold
                  # V03 intended to be able to handle & rotate 2D data
                    #00 re-organisation of the function for the V03 purpose. We organize the 'threshold placement' into a separate function so we can
                        # call it before and after rotation. rotation will be done using a robust linear modeling of the largest population
                    #01 is a complete overhaul of the function structure in an attempt to de-spaghettify large portions of the code
                        # the input mode has been changed & recursion has been added
                        # several minor bugs have been removed
                    #02 changed the rotation to pivot around the median droplet FUs
                        # added neg.ref as an effort to make a better decision on 'single populations'
                        # updated graphics with ggpubr
                    #03 further graphical updates (all plots now use ggplot, which makes transfer to shiny more straightforward
                        # perfplots has been removed, it was not getting the point accross anyway
                    #04 addition of the 'well' input variable allowing to run the algorithm on only a few reactions in a set
                        # corrected some of the input variable descriptions
                        # adjustments to 'X11()' call to get rid of warning messages
                        # fixed an issue where a certain scenario would lead to plotting errors.
## TO DO ##
 # add a way to check the p-value of the slope, a lot of 'good' analysis don't benefit from rotation, so it should be auto-skipped
 # else the double call to treshfinder will only slow the analysis
 # add rotation to output (i.e the slope)
 # quick check if rotation improved Rs, if not, undo?
 # ADD NEGREF MESSAGING

 # currently functional: simplex, eva
 # next up: simplex2 (simplex analysis of channel2)
 # and ultimately duplex, probably needs another rotation. we may postpone that to version 4 as we have no reactions to test it on

 # maybe add a check to see if both matrices have the same size in ch1 and ch2

## Main function ###############
#########################################################################################################
cloudy <- function(drp, well = "all", method = c("simplex", "simplex2", "eva", "duplex")[1], neg.ref = 13500, threshold = NA,
                   dVol = 0.85, sVol = 20, plots = FALSE, silent = TRUE, vec = FALSE, ...){
   require(SuppDists)# required for the 'moments()' function
   require(MASS)     # required for the 'lqs()' function
   require(ggplot2)  # required for graphics
   require(ggpubr)   # required for graphics
   require(gridExtra)# required for graphics
##################################################################################
#Well selection
##################################################################################
if(all(well == "all")){
     break()
   }else{
     if(is.character(well)){
        if(any(well != "all")){
          if(all(well %in% colnames(drp$Ch1))){
            well <- sapply(well, function(x) which(x == colnames(drops$Ch1)))
          }else{
            message("WARNING: well names not found. Analyzing ALL data")
     }}
     }else{ #assuming well is numeric
          if(!all(well %in% seq(from = 1, to = ncol(drops$Ch1), by = 1))){
            message("WARNING: well numbers out of bounds. Analyzing ALL data")
            well <- seq(from = 1, to = ncol(drops$Ch1), by = 1)
     }}
    drp <- list(Ch1 = drp$Ch1[, well], Ch2 = drp$Ch2[, well])
   }

##################################################################################
#Recursion
##################################################################################
if(is.matrix(drp$Ch1)){
  multiput <- matrix(NA, ncol = ncol(drp$Ch1), nrow = 14)# allocate
  rownames(multiput) <- c("targ.in.sVol", "targ.upper", "targ.lower",
                          "lambda", "lambda.upper", "lambda.lower",
                          "resolution", "%rain", "%comparted",
                          "positive", "negative", "rain", "populations", "threshlold")
  for(s in 1:ncol(drp$Ch1)){
    subdrp <- list(Ch1 = drp$Ch1[, s], Ch2 = drp$Ch2[, s])
    multiput[, s] <- cloudy(subdrp, method, neg.ref, threshold, dVol, sVol, plots, silent, vec = TRUE)
  }
  return(multiput)
}else{
##################################################################################
#input consistency check
##################################################################################

checks <- c(length(method == 1), method %in% c("simplex", "simplex2", "eva", "duplex"), length(c(dVol, sVol)) == 2, is.numeric(dVol) & is.numeric(sVol),
             all(is.logical(c(plots, silent, vec))), length(c(plots, silent, vec)) == 3, is.na(threshold) | is.numeric(threshold))

## only continue into function unless all are true
if(all(checks)){

##################################################################################
#colorschemes for the colorblind##################################################
     two.color <- c(rgb(68,  140,  59, maxColorValue = 255), rgb(117, 170,  93, maxColorValue = 255), rgb(160, 204, 124, maxColorValue = 255),
                    rgb(219, 255, 182, maxColorValue = 255), rgb(216, 196, 213, maxColorValue = 255), rgb(192, 160, 202, maxColorValue = 255),
                    rgb(171, 124, 189, maxColorValue = 255), rgb(133,  84, 168, maxColorValue = 255))
     spectral  <- c(rgb(241,  17,  39, maxColorValue = 255), rgb(255, 119,  40, maxColorValue = 255), rgb(255, 163,  41, maxColorValue = 255),
                    rgb(253, 208,  42, maxColorValue = 255), rgb(245, 255,  46, maxColorValue = 255), rgb(139, 205, 255, maxColorValue = 255),
                    rgb(132, 132, 223, maxColorValue = 255), rgb(123, 168, 194, maxColorValue = 255), rgb(110,   9, 123, maxColorValue = 255))
##################################################################################
#small internal functions:
##################################################################################
     s.funk <- function(k){return(3.8 + 0.35 * log(k) + 0.045 * log(k)^2 +0.75)} # kurtosis to spread
     nonas <- function(d){return(d[!is.na(d)])} # remove NAs from vector
##################################################################################
#Internal function: PEAKFINDER (finds all local minima/maxima in a range of data )
##################################################################################
     # changes have been made to findpeaks to work in this specific context, it is no longer a general purpose function
     # i.e. it makes certain assumptions: the start and end of the 'spectrum' are baseline values
     #      it returns valleys for the purpose of delimiting the peaks it finds, so downstream peak calculations can happen
     #      between these valley limits. There may technically be more valleys located throughout the spectrum, but we don't care.
     findpeaks <- function(vec, bw = 1, trim = 0.01, ...){  #where bw = is box width, setting the sensitivity of the search
                    ### initiate vectors
                   	pos.x.max <- NULL ;	pos.x.min <- 1
                    ###Start of for loop:    we walk down the vector with a window of size "bw"
                    for(i in 1:(length(vec)-1)){
                      #check if we have reached the end of the vector
                      if((i+1+bw)>length(vec)){sup.stop <- length(vec)}else{sup.stop <- i+1+bw}
                      #check if we are at beginning of the vector
                  		if((i-bw) < 1){inf.stop <- 1}else{inf.stop <- i-bw}
                      #select window in two parts: values beyond i (superior), and values before i (inferior)
                  		subset.sup <- vec[(i+1):sup.stop]
                  		subset.inf <- vec[inf.stop:(i-1)]
                    ##############################################################
                  		is.max   <- sum(subset.inf > vec[i]) == 0    #are ALL trailing data smaller than i?
                   		is.nomin <- sum(subset.sup > vec[i]) == 0    #are ALL leading data smaller than i?
                   		no.max   <- sum(subset.inf > vec[i]) == length(subset.inf) #are ALL trailing data larger than i?
                  		no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup) #are ALL leading data larger than i?
                    ##############################################################
                      #a maximum is found if  all data before and after i are smaller than i
                  		if(is.max & is.nomin){ pos.x.max <- c(pos.x.max, i) }#
                      #a minimum is found if  all data before and after i are smaller than i
                  		if(no.max & no.nomin){ pos.x.min <- c(pos.x.min, i) }#
                     }#end of for loop
                    ### trim off tiny peaks
                      pos.x.max <- pos.x.max[vec[pos.x.max] > max(vec[pos.x.max])*trim]
                    ###Consolidate peak-valley matrix  # 2 indicates max, 1 indicates valley
                      pvmat <- matrix(c(pos.x.max, pos.x.min, rep(2, length(pos.x.max)), rep(1, length(pos.x.min))),
                                      nrow = 2, ncol = length(c(pos.x.max, pos.x.min)), byrow = T)
                      pvmat <- pvmat[, order(pvmat[1,])] # order
                     # we want n + 1 valeys for n peaks (start & then one valley after avery peak)
                     # if there are multiple possibilities we want the valley centered between 2 peaks
                     check <- 1 ; j <- 1 ; busy <- TRUE
                     while(busy){ # loop over pv matrix, since we'll be adding columns to the matrix we use a while loop
                       if(pvmat[2, j] != check){ # we do not find what we expect!
                         if(pvmat[2, j] == 1){ # previous was also valley, merge both valleys
                            # move  position to halfway between the two valeys, except if j = 2
                            if(j != 2){ # the first valley should always be the vector start
                             pvmat[1, j - 1] <- pvmat[1, j-1] + floor((pvmat[1, j] - pvmat[1, j-1])/2) }#
                            # remove current position (we've updated the previous valley)
                             pvmat <- pvmat[, -j]
                            # since we removed the col we should not advance j
                         }else if(pvmat[2, j] == 2){  # previous was also peak, add a valley
                             # use minimum value between this peak and the previous one
                              pvmat <- cbind(pvmat[, 1:(j-1)],
                                    matrix(c(pvmat[1, j-1] + which.min(vec[(pvmat[1, j-1] + 1):(pvmat[1, j] - 1)]), 1),
                                    ncol = 1),pvmat[, -c(1:(j-1))])
                              j <- j + 2
                         }# END check
                       }else{ # order is respected, switch check & advance
                         check <- c(2 ,1)[check]
                         j <- j + 1
                       }# END if != check
                      if(j > ncol(pvmat)){ busy <- F }# check if we reached the end of the matrix
                     }# END while
                    ###Final check: last element shoud be a valley, if so set to vector end, if not add it
                    if(tail(pvmat[2,], 1) == 1){
                      pvmat[1, ncol(pvmat)] <- length(vec)
                    }else{
                      pvmat <- cbind(pvmat, matrix(c(length(vec), 1), ncol = 1))
                    }#END if tail
                 ###Output
                  pea <- pvmat[1, pvmat[2, ] == 2]
                  val <- pvmat[1, pvmat[2, ] == 1]
                  return(list("peaks" = rbind(pea, vec[pea]), "valleys" = rbind(val, vec[val])))
     }#End of peakfinder function########################################################

##################################################################################
#Internal function: THRESHFINDER (finds population boundary limits)
##################################################################################
    # it is important that findthresh is defined here (i.e. internal to cloudy) otherwise
    # the '<<-' operator saves to the wrong environment as functions retain a link to the
    # environment in which they were created.
    findthresh <- function(drops, silent, threshold, ...){  #
      skip <- 0
      ##################################################################################
      # KERNEL DENSITY SPLIT
      ##################################################################################
      #be sure bandwidth for kernel density is at least 50
      temp <- bw.nrd0(drops)
       if(temp < 50){
          temp <- 50
       }
        #kernel density estimation
         krn <- density(drops, bw = temp)
         krn <- rbind(krn$y, krn$x)
         rownames(krn) <- c("y", "x")
        #we use findpeaks with a high bandwith to single out the populations
         piek <- findpeaks(krn[1, ], bw=20)
      ##################################################################################
      ##First major control point: we only proceed if there are no NA
      if(!any(is.na(unlist(piek)))){ ##Control LVL 1
      ##################################################################################
      if(length(piek$peaks[1,]) == 1){ # only one majore peak was found, we'll check if we can deduce if it's pos or neg
         #Since there still may be few negatives / positives that got smoothed out.
         #We should be able to detect them as 'outliers' from the single population
         #therefore we calculate z-scores for all droplets
         zs <- (drops - median(drops))/mad(drops)
         s1 <- s.funk(abs(moments(drops)[4]))
         # check for outliers and check if they are higher or lower to decide if the pop is pos or neg
         # if we have BOTH positive and negative outliers we cannot make a meaningful decision without user input
         # in the latter case, we'll just give up. Also, we must make sure there is more than one outlier,
         # otherwise none of the below makes sense. In fact, we'll demand at least 3.
         outzs <- abs(zs) > (1.5 * s1)
         # pos and neg outzs may exist, we'll choose one. Either based on the neg.ref, or if it is NA we use
         # decide based on the number of observations. The latter is not an ideal approach, but it should
         # work for 99% of the cases where this question pops up
         if(!is.na(neg.ref)){
           if(median(drops) > neg.ref){ # ALL positive, we focus on the lower outliers (possibly negative)
             outzs[outzs & zs > 0] <- F # set all positive outzs to F
           }else{
             outzs[outzs & zs < 0] <- F # set all negative outzs to F
           }#
         }else{ # no neg.ref, we decide based on number of outliers
           if(sum(zs[outzs] > 0) > sum(zs[outzs] < 0)){
              #there are more positives
             outzs[outzs & zs < 0] <- F # set all negative outzs to F
           }else{
              #there are more negatives
             outzs[outzs & zs > 0] <- F # set all positive outzs to F
           }# END if-else
         }# END if-else

         if(length(zs[outzs]) > 2){
          #we have outliers! we'll need to update piek so we can proceed as normal
          if(all(zs[outzs] > 0)){ #outliers have higher fluorescence and are thus positive
             popmean <- mean(drops[outzs])
             #retrieve position in krn
             popkrn <- which.min((krn[2, ] - popmean)^2)
             #find lowest spot between the 2 peaks
             minkrn <-  piek$peaks[1,1] + 1 + which.min(krn[1, (piek$peaks[1,1] + 1):(popkrn - 1)])
             #update pieks & valleys
             piek$peaks <- cbind(piek$peaks, c(popkrn, krn[1, popkrn]))
             piek$valleys <- cbind(piek$valleys[,1], c(minkrn, krn[1, minkrn]), piek$valleys[,2])
          }else{ #outliers have lower fluorescence and are thus negative
             popmean <- mean(drops[outzs])
             #retrieve position in krn
             popkrn <- which.min((krn[2, ] - popmean)^2)
             #find lowest spot between the 2 peaks
             minkrn <-  popkrn + 1 + which.min(krn[1, (popkrn + 1):(piek$peaks[1,1] - 1)])
             #update pieks & valleys
             piek$peaks <- cbind(c(popkrn, krn[1, popkrn]), piek$peaks)
             piek$valleys <- cbind(piek$valleys[,1], c(minkrn, krn[1, minkrn]), piek$valleys[,2])
          }
        ### we have no outliers, but at least the threshold allows us to see if the drops are all pos or neg
        #############################
         }else if(!is.na(threshold)){
           pops <<- 1
           tresh <<- threshold
           # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
           piik <- matrix(c(threshold - 100, threshold + 100, min(drops) - 100, max(drops) + 100, threshold -0.1, threshold + 0.1),
                          nrow = 3, ncol = 2, byrow = T)
           rownames(piik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
           # we should, however, not attempt further analysis of the data or we screw up our ersatz piik
           skip <- 1
        ### we have no outliers, but at least the neg.ref may allows us to see if the drops are all pos or neg
        #############################
         }else if(!is.na(neg.ref)){
           pops <<- 1
           if(median(drops) > neg.ref){ # ALL positive, we need at least 1 negative
             tresh <<- mean(tail(sort(drops, decreasing = T), 2)) # put threshold between two lowest droplets
             # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
             piik <- matrix(c(min(drops), median(drops), min(drops) - 100, max(drops) + 100, tresh -0.1, tresh + 0.1),
                            nrow = 3, ncol = 2, byrow = T)
             rownames(piik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
             }else{ # all negative, zero positives is perfectly ok, we set the threshold above the max droplet
             tresh <<- max(drops) + 2 * mad(drops)
             # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
             piik <- matrix(c(median(drops), tresh + 50, min(drops) - 100, tresh + 100, tresh - 0.1, tresh + 0.1),
                            nrow = 3, ncol = 2, byrow = T)
             rownames(piik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
             } # END if-else
           # we should, however, not attempt further analysis of the data or we screw up our ersatz piik
           skip <- 1
        ### no outliers, no threshold: we give up
        #############################
         }else{
           piik <- matrix(NA, ncol = 2, nrow = 3)  # make empty piik matrix
           rownames(piik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
           pops <<- 1
           fail <<- 1 # notify failure
           skip <- 1 # make sure we skip
           badluck[2] <<- TRUE
         } #end of no outlier case
       }# end of piek == 1 check
     ## #check if we should continue or give up#########################################
     if(fail == 0 & skip == 0){#####################################################################
       if(length(piek$peaks[1,]) >= 3){ #give a warning if we find more than 2 populations
          #mark warning
          badluck[3] <<- TRUE
        }# END if peaks >= 3
       #*universal* population delimiter
       ################################
       #set the number of populations
        pops <<- length(piek$peaks[1,])
       # set the peak limits using valleys defined by peakfinder , row 1: limits low FU peak, row 2: limits high FU peak
        pea.lim <- matrix(krn[2, c(head(piek$valleys[1,], 2), tail(piek$valleys[1,], 2))], nrow = 2, byrow = T) # in FU

       #get started on the main population matrix
        piik <- matrix(c(krn[2, piek$peaks[1, range(seq_along(piek$peaks[1,]))]], # FU of first and last peak
                         pea.lim[1,1], pea.lim[2,2], pea.lim[1,2], pea.lim[2,1]),ncol = 2, byrow = T)
        # for clarity we add row names:
        rownames(piik) <- c("x.FU",        # first row is population center position in Fluorescence units
                            "outer.bounds",# second row is the outer limits (lowest and hightst FU over ALL)
                            "inner.bounds")# third row is the inner population bounds (between the populations)
         # I know the order of the values in piik is not completely logical & I have no idea why i did it this way
         # I suppose I had a valid reason when I first wrote this & I'll keep it this way for compatibility's sake

       #allright, using this maximal delimination of the clouds we can calculate more accurate values
       #by iterating the update process three times (i.e. use current estimate to calculate pos and neg,
       # then recaculte sd and re-delimite bands). Three iterations is generaly enough to reach
       # stability in the first decimal place.
         for ( i  in 1 : 3 ) {
          #split in positive and negatives
           pdrp <- drops[drops >= piik[3, 2] & drops <= piik[2, 2]]
           ndrp <- drops[drops >= piik[2, 1] & drops <= piik[3, 1]]
           #update s multiplier
           s1 <- s.funk(abs(moments(ndrp)[4]))
           s2 <- s.funk(abs(moments(pdrp)[4]))
           #back-up
           piik.bup <- piik
           #update
           piik[2, 1] <-  piik[1, 1] - mad(ndrp) * s1 #leftmost extreme (start of neg band)
           piik[2, 2] <-  piik[1, 2] + mad(pdrp) * s2 #rigthmost extreme (end of pos band)
           piik[3, 1] <-  piik[1, 1] + mad(ndrp) * s1 #first inner extreme (end of the neg band)
           piik[3, 2] <-  piik[1, 2] - mad(pdrp) * s2 #second inner extreme (start of the pos band)
           # check if we didn't screw up (which may happen with very few pos/neg)
           if(any(is.na(piik))){
             #maybe a warning?
             piik <- piik.bup # error error, restore
             break # stop optimization
           }# END NA check
          }#END of the iterative update

         #If there is a threshold specified we'll use that, otherwise we calculate our own
         if(is.na(threshold)){
           #We now calculate the threshold for counting positives
            tresh <<- piik[1, 1] + (s1 + 1/2 * s1) * mad(ndrp)
           #should the threshold fail, we provide a backup
            tresh.bup1 <<- pea.lim[1, 2]
         }else{
            tresh <<- threshold
         }#END if is.na(threshold)
      # a last check before we return the output:
      ##############################################################################
      ###Control LVL 2 : a problem related to standard deviations is overlapping bounderies
      #check : upper neg boundary is situated above lower pos boundary?
      if(piik[3, 1] > piik[3, 2]){##Control LVL 5
          badluck[4] <<- TRUE
          #we can use the 'valleys' returned by findpeaks to draw new population bounderies
          piik[3, 1] <- pea.lim[1, 2]
          piik[3, 2] <- pea.lim[2, 1]
      }##Control LVL 2 END
       } #END of fail== 0 & skip == 0
      ################################
      # END of the *universal* population delimiter
      ##############################################################################
      }else{ ##Control LVL 1
        # to arrive here, piek must contain NA values. All is lost, throw an error and set failure to 1.
        badluck[1] <<- TRUE
        piik <- matrix(NA, ncol = 2, nrow = 3)
        rownames(piik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
        fail <<- 1 #notify failure
      }##Control LVL 1 END
      ##############################################################################
     return(piik)
     }# FINDTRESH END
################################################################################
#messaging system
     public.address <- c("WARNING: something is not right, please check data \n",    #1
                         "WARNING: only one population detected, returning NA",      #2
                         "WARNING: multiple population detected \n",                 #3
                         "WARNING!: population boundary overlap \n   boundaries rewritten, not based on population parameters \n", #4
                         "Rotation seems to have screwed up! ...restoring droplets... \n", #5
                         "WARNING!: threshold placement failed \n   Check data, arbitrary threshold was set \n", #6
                         "WARNING!: threshold inside positive cloud \n   Forced lower threshold, not based on population parameters \n" #7
                       )##
##################################################################################
## ACTUAL FUNCTION START #########################################################
##################################################################################
   # Choose your adventure:
    if(method == "simplex"){ # simplex analysis of channel 1, no channel 2 needed
      drp1 <- na.omit(drp$Ch1)
    }else if(method == "simplex2"){ # simplex analysis of channel 2
      drp1 <- na.omit(drp$Ch2)
    }else if(method == "eva"){ # analysis of channel 1 with crosstalk correction using channel 2
      drp1 <- na.omit(drp$Ch1)
      drp2 <- na.omit(drp$Ch2)
    }#
   # create some values for the internal function to update (<<-)
     fail <- 0 #create failure tracker
     pops <- 0 #create population tracker
     tresh      <- NA
     tresh.bup1 <- NA
     badluck <- rep(FALSE, times = 7)

   # call to treshfinder function
     piik <- findthresh(drp1, silent, threshold)   # TEMP !! change piik <<- back to <- for final version!!

   # use piik to define the populations
     pospop <- drp1 >= piik[3, 2] & drp1 <= piik[2, 2]
     negpop <- drp1 >= piik[2, 1] & drp1 <= piik[3, 1]

##################################################################################
# Rotation!
##################################################################################
# note that at this point the initial call to threshfinder may have failed and maybe fail == 1
# However, if rotate == T we may still have a chance at making an analysis by improving contrast.
# We'll simply use all droplets in the regression (since they are so close together that they appear undescernable anyway)
    if(method %in% c("eva", "duplex")){
      # check if fail == 1
      if(fail == 1){
        targpop <- rep(T, length(drp1)) # use all droplets!
        # we'll also reset fail
        fail <- 0
      }else{ # no fail, all is good, we use the populations defined by piik
        # we take the largest population in channel1 to make our correlation analysis
        targpop <- if(sum(pospop) > sum(negpop)){ pospop }else{ negpop }#
      }# END if fail == 1
      # robust linear regression using LTS
      lq.drp <- lqs(drp1[targpop] ~ drp2[targpop], method = "S")
        # backup in case the rotation doesn't work out and we need to restore
        drp1.bup <- drp1
        piik.bup <- piik
        pops.bup <- pops
        neg.ref.bup <- neg.ref
      # actual correction
      drp1 <- drp1 - drp2 * coef(lq.drp)[2] + median(drp2) * coef(lq.drp)[2] # to prevent changing the FU too much we add the median value
      neg.ref <- neg.ref * median(drp1[targpop])/median(drp1.bup[targpop]) # we try to keep the distance from the neg.ref to the targpop the same (relatively)
      # another call to threshfinder
      bad.bup <- badluck # backup
      badluck[1:4] <- rep(FALSE, times = 4) #reset badluck
      piik <- findthresh(drp1, silent, threshold)

       # on the off chance that the rotation totally screwed up the data we need to undo all the changes!
        if(!any(is.na(piik.bup)) & fail == 1){ ## if that is the case, fail will be 1, in all other cases a valid piik will be produced
                                                # actually, that's not true, threshfinder may have failed before, so we only restore if piik.bup
                                                # does not contain any NA. Otherwise we keep fail == 1 & give up
           badluck[5] <- TRUE
           # restore variables
           badluck <- bad.bup
           drp1 <- drp1.bup
           piik <- piik.bup
           pops <- pops.bup
           neg.reg <- neg.ref.bup
           rm(piik.bup, drp1.bup, pops.bup, bad.bup, neg.ref.bup) # cleanup
           # reset fail
           fail <- 0
        }# END fail check

      # use the new and improved piik to re-define the populations
      pospop <- drp1 >= piik[3, 2] & drp1 <= piik[2, 2]
      negpop <- drp1 >= piik[2, 1] & drp1 <= piik[3, 1]
        # show a quick plot of before and after rotation (IF rotation was succesful):
           if(isTRUE(plots) & exists("drp1.bup")){
            #x11(14,7); par(mfrow = c(1,2))
             #plot(drp1.bup ~ drp2, col = spectral[6], bty = "L", pch = 20 , xlab = "Channel2", ylab = "Channel1", ylim = c(min(drp1.bup), max(drp1.bup)), main = "before")
             #abline(lq.drp)
             #abline(h = neg.ref.bup, lty = 2)
            #plot(drp1 ~ drp2, col = spectral[6], bty = "L", pch = 20 , xlab = "Channel2", ylab = "Channel1", ylim = c(min(drp1), max(drp1)), main  = "after")
             #abline(h = neg.ref, lty = 2)
        # Grouped Scatter plot with marginal density plots
        ggscatterhist(data.frame(ch1 = drp1, ch2 = drp2), x = "ch2", y = "ch1")
        }# END plots
    }#END if isTRUE rotator

  #  print("check1")

##################################################################################
## we can only continue if the threshfinder was able to define the populations
## otherwise we skip to the end.
if(fail == 0){ ##Control LVL 2
   # split droplets in positives and negatives
   pdrp <- drp1[pospop]
   ndrp <- drp1[negpop]
# that's it, we should have everything to produce an output
# from here on it's a short run to the end of the function
# however, we should do some additional checks to make sure:
##################################################################################
   ##third to fifth control points: checks for stability reasons
   ##################################################################################
   #in case of failure (tresh is NA) we fall back on the backup threshold
   if(is.na(tresh)){##Control LVL 3
        badluck[6] <- TRUE
        tresh <-  tresh.bup1
        ndrp <- drp1[drp1 < tresh]
        pdrp <- drp1[drp1 > tresh]
        ## if tresh is NA the piik may be NA as well, check & redo all the calculations:
        if(any(is.na(piik))){
        piik <- c(median(ndrp),median(pdrp))
        #Expand piik to take all values
        piik <- rbind (piik, matrix(c(min(ndrp),max(pdrp),max(ndrp),min(pdrp)),nrow=2,ncol=2, byrow = T))
        }# END piik rescue
    }##Control LVL 3 END
    ##############################
    # the next control is skipped if a threshold is specified (we don't want it to be overruled, even if it's wrong)
    if(is.na(threshold)){
    ##############################
    #excessive rain can make the standard deviations go haywire
    #we build in this simple check to make sure the treshold is not placed into the positive band
    if(tresh > piik[3, 2]){ ##Control LVL 4
        badluck[7] <- TRUE
        tresh <-  tresh.bup1
    }##Control LVL 4 END
    }# end if(is.na(threshold))

##################################################################################
#we can now start the trivial calculations:
##################################################################################
    #count the raindrops (i.e. the population of the 'clearance band' between the positive and negative band)
     #again to avoid including negatives in the rain we use the treshold in stead of the peak base
     #rain <- length(which(drp > piik[3, 1] & drp < piik[3, 2]))
     # we do not calculate rain in case of 3 or more populations
     # the middle populations end up as 'rain' and give a biased metric
     # Since we allow for manual threshold placement it is now possible that
     # reactions with 1 population arrive at this point. However,
     # we also cannot calculate the rain if there is only 1 population
     if(pops >= 3 | pops == 1){
      rain <- NA
     }else{
      rain <- sum(drp1 > tresh & drp1 < piik[3, 2])
     }#
    # make a dataframe to label the positive and negatives so we can apply colours on the plot
     dfdrops <- data.frame(FU = drp1, positive = drp1 > tresh, category = 0)
     tdrp <- length(drp1) # total number of droplets

     # the following should only be calculated if there is more than 1 population
     if(pops > 1){
      #also define 'cerntainly pos', 'rain' and 'drift'
       dfdrops[!dfdrops$positive, 'category'] <- 1  # negative droplets
       dfdrops[dfdrops$FU > tresh & dfdrops$FU < piik[3, 2], 'category'] <- 2    #rain is the population of the clearance band (counted positive)
       dfdrops[dfdrops$FU > piik[3, 2] & dfdrops$FU < piik[2, 2], 'category']  <- 3 #c-pos is the population of the pos band only
       dfdrops[dfdrops$FU > piik[2, 2], 'category']  <- 4  #drift is population above the pos band (high-Fluorescence outliers)
 ###Calculation of Performance parameters
      #1#  Resolution
       Resol <- 2 * (piik[1, 2] - piik[1, 1]) / ((piik[3, 1] - piik[2, 1]) + (piik[2, 2] - piik[3, 2]))
      }else{ # in case of only 1 population Resolution & so on should still be defined
       Resol <- NA
       dfdrops[!dfdrops$positive, 'category'] <- 1  # negative droplets
       dfdrops[dfdrops$positive, 'category'] <- 3  # positive droplets
      }# END pops > 1

    #2# percentage rain
     pRain <- rain / tdrp
    #3# Compartmentalization
     Comp <- tdrp * dVol * 0.001 / sVol
 ###Calculation of Lambda and confidence bounds
     lambda <- -log(sum(!dfdrops$positive)/tdrp)
      lamlo <- lambda - 1.96 * sqrt((tdrp - sum(!dfdrops$positive))/(tdrp * sum(!dfdrops$positive)))   #lower 95% confidence bound for lambda 1 (cpd)
      lamhi <- lambda + 1.96 * sqrt((tdrp - sum(!dfdrops$positive))/(tdrp * sum(!dfdrops$positive)))   #upper 95% confidence bound for lambda 1 (cpd)
      lambs <- c(lambda, lamhi, lamlo)
  #total number of targets present in droplets analysed
    targets <- tdrp * lambs
  #total number of targets in sample volume
    targall <- lambs * (sVol/(dVol * 0.001))
 ###gather and name all output
      names(targets) <- c("targets", "upper", "lower")
      names(targall) <- c("targ.in.sVol", "upper", "lower")
      names(lambs) <- c("lambda", "upper", "lower")
     drops <- c(sum(dfdrops$positive), sum(!dfdrops$positive), rain)
      names(drops) <- c("positive", "negative", "rain")
     perform <- c(Resol, pRain, Comp)
     names(perform) <- c("resolution", "%rain", "%comparted")
     names(pops) <- "populations"
     names(tresh) <- "threshold"

##################################################################################
### Make plots (optional)
##################################################################################
    if (isTRUE(plots)) {
     # since we're using ggplot, we'll have to turn piik into a dataframe for convenient plotting
     dfpiik <- as.data.frame(t(piik))
     dfpiik$category <- c(1,3)
     # we also have to add a random order to dfdrops to plot them randomly
     dfdrops$ordr <- sample(c(1:nrow(dfdrops)), size = nrow(dfdrops), replace = F)
     # unfortunately our plan for using geom_density does not work out well, both populations get scaled differently which gives a very misleading plot
     # so we must make a separate plotting data frame:
     dens <- density(drp1)
     plotdf <- data.frame(x=dens$x, y=dens$y)
     plotdf$dropsplit <- factor(findInterval(plotdf$x, tresh))
     #1# kernel density
      #plot1 <-  ggplot() + geom_density(data = dfdrops, aes(x = FU, fill = positive)) +  geom_density(alpha=0.4) +
      #         scale_fill_manual(labels = c("negative", "positive"), limits = c('FALSE', 'TRUE'), values = c('#F8766D', '#00BFC4'), name = "Population") +
      plot1 <- ggplot(plotdf, aes(x,y)) + geom_ribbon(aes(ymin=0, ymax=y, fill=dropsplit)) +  geom_line() +
               scale_fill_manual(labels = c("negative", "positive"), limits = c(0, 1), values = c('#F8766D', '#00BFC4'), name = "Population") +
               labs(tag = "A", x = "FU", y = "density")
      if(fail == 0){ #piik lines can only be added if piik is not NA
               plot1 <- plot1 + geom_vline(data = dfpiik, aes(xintercept = x.FU), linetype = "dashed") +
                     geom_vline(data = dfpiik, aes(xintercept = outer.bounds), linetype = "dashed") +
                     geom_vline(data = dfpiik, aes(xintercept = inner.bounds), linetype = "dashed")
      }#

     #2# droplet view
      plot2 <- ggplot() + geom_point(data = dfdrops, aes(y = FU, x = ordr, colour = category)) +
               scale_colour_gradientn(breaks = c(1, 2, 3, 4), colors = c("#F8766D", "#FFC425", "#00BFC4", "#0076C4"), labels = c("negative", "rain", "positive", "outlier"), limits = c(0.9,4.1), name = "Population") +
               geom_hline(yintercept = tresh, linetype = "dashed") +
               labs(tag = "B", x = "Observation")

     ### put on screen
      x11(width = 14, height = 6)
      grid.arrange(plot1, plot2, nrow = 1)

  }# END plots

##################################################################################
#Standard INFO messages
##################################################################################
    if(!isTRUE(silent)){
      if(isTRUE(Resol < 2.5)){  message("Resolution: FAILED")}else{message("Resolution: PASSED")} #isTRUE wrapper as Resol can be NA
      if(is.na(pRain)){message("Rain: multiple populations, rain not calculated")
                      }else if(pRain < 0.025){message("Rain: PASSED")
                      }else{message("Rain: FAILED")}
      if(Comp < 0.3){message("Compartmentalization: FAILED")}else{
      if(Comp < 0.5){message("Comparmentalization: MEDIOCRE")}else{
                     message("Comparmentalization: PASSED")}}
      if(lambs[1] < 0.01){message("Confidence bounds: NOT optimal")}else{
      if(lambs[1] < 5){message("Confidence bounds: OPTIMAL")}else{
                     message("Comparmentalization: NOT optimal")}}
      }#end of messages

##################################################################################
#ERROR messages
##################################################################################
   }else{
    message("WARNING: treshold finding failed, returning NA")
   }##Control LVL 2 END
}else{
    message("WARNING: input inconsistent, returning NA")
    error.log <- c("You may only choose 1 method, argument should be of length 1 \n",
                   "Available methods: simplex, simplex2, eva, and duplex. Other choices are invalid \n",
                   "dVol & sVol should both be of length 1 \n",
                   "dVol & sVol should both be a numerical \n",
                   "plots, silent, and vec all must be a logical \n",
                   "plots, silent, and vec all should be of length 1 \n",
                   "threshold should be either NA or a numerical of length 1")
    message(error.log[checks])
}## consistency check
## report messages
if(!isTRUE(silent)){
  if(any(badluck)){message("\n --------- \n")}
  message(public.address[badluck])
  }##
##################################################################################
   #Failure (NA) output
   if(fail == 1){
   #write all NA vectors:
    targall <- c(NA, NA, NA)
    names(targall) <- c("targ.in.sVol", "upper", "lower")
    lambs  <- c(NA, NA, NA)
    names(lambs) <- c("lambda", "upper", "lower")
    perform  <- c(NA, NA, NA)
    names(perform) <- c("resolution", "%rain", "%comparted")
    drops <- c(NA, NA, NA)
    names(drops) <- c("positive", "negative", "rain")
    pops <- NA
    names(pops) <- c("populations")
    tresh <- NA
    names(tresh) <- "threshlold"
   } #End of failure (NA) procedure
##################################################################################
# actual output ##################################
##################################################################################
    if(isTRUE(vec)){ #output if
        return(c(targall, lambs, perform, drops, pops, tresh))
    }else{
        return(list("targets.in.sVol" = targall, "lambda" = lambs, "performance" = perform, "droplets" = drops, "populations" = pops, "threshold" = tresh))
    }#output if END
##################################################################################
 }## END of recursion
}##END of 'cloudy'
##################################################################################
##################################################################################
##################################################################################
