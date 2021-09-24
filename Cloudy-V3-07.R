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
                    #05 added a second threhold calculation method based on the Generalized Extreme Value (GEV) approach from Trypsteen et al.
                        # I have implemented what is essentially a less statistically sound application of the GEV approach
                        # it had to be less statistically sound for it to be applied on a per-reaction-basis as I do here
                        # the proper thing to do is to establish the threshold on a set of NTCs, while I use a rough guess of the negatives
                        # as a consequence I had to make some adjustments to deal with 'false negative' droplets (i.e. rain) that mess with the gev model fitting
                        # The main change 'under the hood' of the algorithm is that the threshfinder function is now split into 2 parts:
                          # a rough estimator (also used for gev)
                          # the remainder of the cloudty threhold approach
                        # currently the gev-threshold is not returned as output, it is merely drawn on the plot.
                        # this entire feature may be removed in future versions if it turns out to be unstable
                    #06 implementation of 'trimming' in cloudy threshold function
                        # this should improve population delimitation in asymetric droplet populations
                        # in the previous versions, the population limits were put symetrical around the peak
                        # but the tail of asymetric populations would then cause a large 'empty' region on the other side of the peak
                        # the trimming looks for these empty spots and removes them
                    # 07 proper implementation of duplex analysis
                        # until now analysis of two channels had to be carried out manually as the code would crash when calling duplex
                        # this has been solved and multi-channel analysis is now a thing. Once initiated, the algorithm will loop over all
                        # channels present. Since read QX only reads two channels that is the current limit, support for more channels will
                        # follow but requires re-doing the recursion bit where currently only 2 channels are processed (see line 148).
                        # as a result, the output has also changed sligthly when 2 channels are processed. For a single reaction a row has
                        # been added to the output to accomodate the second channel results. For multiple reactions, instead of returning a
                        # single matrix, it now returns a list with two items (Ch1 and Ch2) each containing a matrix with results for that channel.
                        # also added a progress bar when analysing multiple reactions
## TO DO ##
 # add a way to check the p-value of the slope, a lot of 'good' analysis don't benefit from rotation, so it should be auto-skipped
 # else the double call to treshfinder will only slow the analysis
 # quick check if rotation improved Rs, if not, undo?
 # ADD NEGREF MESSAGING
 # bootstrap approach for GEV threshold in case of few negatives??
 # maybe add a check to see if both matrices have the same size in ch1 and ch2

## Main function ###############
#########################################################################################################
cloudy <- function(drp, well = "all", method = c("simplex", "simplex2", "eva", "duplex")[1], neg.ref = 7500,
                   threshold = c("cloudy", "gev", "mixed")[1], dVol = 0.85, sVol = 20, plots = FALSE, silent = TRUE, vec = FALSE, ...){
##################################################################################
#dependencies
##################################################################################
   require(SuppDists)# required for the 'moments()' function
   require(MASS)     # required for the 'lqs()' function
   require(ggplot2)  # required for graphics
   require(ggpubr)   # required for graphics
   require(gridExtra)# required for graphics
   require(evd)      # required for 'fgev()' and  'qgev()' functions
   require(tcltk)    # required for progressbar
##################################################################################
#Well selection
##################################################################################
if(all(well == "all")){
     # do nothing
   }else{
     if(is.character(well)){
        if(any(well != "all")){
          if(all(well %in% colnames(drp$Ch1))){
            well <- sapply(well, function(x) which(x == colnames(drp$Ch1)))
          }else{
            message("WARNING: well names not found. Analyzing ALL data")
     }}
     }else{ #assuming well is numeric
          if(!all(well %in% seq(from = 1, to = ncol(drp$Ch1), by = 1))){
            message("WARNING: well numbers out of bounds. Analyzing ALL data")
            well <- seq(from = 1, to = ncol(drp$Ch1), by = 1)
     }}
    drp <- list(Ch1 = drp$Ch1[, well], Ch2 = drp$Ch2[, well])
   }
   # note: these two messages are kept out of the public address system and cannot be silenced

##################################################################################
#Recursion
##################################################################################
if(is.matrix(drp$Ch1)){
   # the approach is different depenging on the requested analysis
   # but we can already define the progress bar as it is common to all
   #pb <- txtProgressBar(min = 0, max = ncol(drp$Ch1), style = 2)
   pb <- tkProgressBar("info for you", "Cloudy is busy", min = 0, max = ncol(drp$Ch1))
   if(method != "duplex"){ # analysis of a single channel, keep legacy output
    multiput <- matrix(NA, ncol = ncol(drp$Ch1), nrow = 14, dimnames = list(c("targ.in.sVol", "targ.upper", "targ.lower",
                       "lambda", "lambda.upper", "lambda.lower", "resolution", "%rain", "%comparted", "positive",
                       "negative", "rain", "populations", "threshold"), colnames(drp$Ch1)))# allocate space to save output
    for(s in 1:ncol(drp$Ch1)){
      subdrp <- list(Ch1 = drp$Ch1[, s], Ch2 = drp$Ch2[, s])
      multiput[, s] <- cloudy(subdrp, well = "all", method, neg.ref, threshold, dVol, sVol, plots, silent, vec = TRUE)   # set well to "all" to prevent recursive subsampling
      #setTxtProgressBar(pb, s)
      info <- paste(s, " reactions analysed", sep = "")
      setTkProgressBar(pb, s, info)
    }# for END
   }else{
    multiput     <- list()
    multiput$Ch1 <- matrix(NA, ncol = ncol(drp$Ch1), nrow = 14, dimnames = list(c("targ.in.sVol", "targ.upper", "targ.lower",
                          "lambda", "lambda.upper", "lambda.lower", "resolution", "%rain", "%comparted", "positive",
                          "negative", "rain", "populations", "threshold"), colnames(drp$Ch1)))# allocate space to save output
    multiput$Ch2 <- matrix(NA, ncol = ncol(drp$Ch1), nrow = 14, dimnames = list(c("targ.in.sVol", "targ.upper", "targ.lower",
                          "lambda", "lambda.upper", "lambda.lower", "resolution", "%rain", "%comparted", "positive",
                          "negative", "rain", "populations", "threshold"), colnames(drp$Ch1)))# allocate space to save output
    for(s in 1:ncol(drp$Ch1)){
      subdrp <- list(Ch1 = drp$Ch1[, s], Ch2 = drp$Ch2[, s])
      subres <- cloudy(subdrp, well = "all", method, neg.ref, threshold, dVol, sVol, plots, silent, vec = TRUE)   # set well to "all" to prevent recursive subsampling
      multiput$Ch1[, s] <- subres[1, ]
      multiput$Ch2[, s] <- subres[2, ]
      #setTxtProgressBar(pb, s)
      info <- paste(s, " reactions analysed", sep = "")
      setTkProgressBar(pb, s, info)
    }# for END
   }# if else END
    close(pb) # close progress bar
    return(multiput)
}else{
##################################################################################
#input consistency check
##################################################################################

checks <- c(length(method == 1), method %in% c("simplex", "simplex2", "eva", "duplex"), length(c(dVol, sVol)) == 2, is.numeric(dVol) & is.numeric(sVol),
             all(is.logical(c(plots, silent, vec))), length(c(plots, silent, vec)) == 3, threshold %in% c("cloudy", "gev", "mixed") | is.numeric(threshold), length(threshold) == 1)

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
     s.funk <- function(k){return(3.8 + 0.35 * log(k) + 0.045 * log(k)^2 + 0.75)} # kurtosis to spread
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
#Internal function: ROUGH.CUT (finds approximate population boundaries)
##################################################################################
    # it is important that the threshold functions are defined here (i.e. internal to cloudy) otherwise
    # the '<<-' operator saves to the wrong environment as functions retain a link to the
    # environment in which they were created.
    rough.cut <- function(drops, silent, threshold, ...){  #
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
         # pos and neg outzs may exist, we'll choose one. Either based on the neg.ref, or if it is NA we
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
           tresh.cld <<- threshold
           # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
           pjik <- matrix(c(threshold - 100, threshold + 100, min(drops) - 100, max(drops) + 100, threshold -0.1, threshold + 0.1),
                          nrow = 3, ncol = 2, byrow = T)
           rownames(pjik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
           # we should, however, not attempt further analysis of the data or we screw up our ersatz piik
           skip <<- 1
        ### we have no outliers, but at least the neg.ref may allows us to see if the drops are all pos or neg
        #############################
         }else if(!is.na(neg.ref)){
           pops <<- 1
           if(median(drops) > neg.ref){ # ALL positive, we need at least 1 negative
             tresh.cld <<- mean(tail(sort(drops, decreasing = T), 2)) # put threshold between two lowest droplets
             # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
             pjik <- matrix(c(min(drops), median(drops), min(drops) - 100, max(drops) + 100, tresh.cld -0.1, tresh.cld + 0.1),
                            nrow = 3, ncol = 2, byrow = T)
             rownames(pjik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
             }else{ # all negative, zero positives is perfectly ok, we set the threshold above the max droplet
             tresh.cld <<- max(drops) + 2 * mad(drops)
             # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
             pjik <- matrix(c(median(drops), tresh.cld + 50, min(drops) - 100, tresh.cld + 100, tresh.cld - 0.1, tresh.cld + 0.1),
                            nrow = 3, ncol = 2, byrow = T)
             rownames(pjik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
             } # END if-else
           # we should, however, not attempt further analysis of the data or we screw up our ersatz piik
           skip <<- 1
        ### no outliers, no threshold: we give up
        #############################
         }else{
           pjik <- matrix(NA, ncol = 2, nrow = 3)  # make empty piik matrix
           rownames(pjik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
           pops <<- 1
           fail <<- 1 # notify failure
           skip <<- 1 # make sure we skip
           badluck[2] <<- TRUE
         } #end of no outlier case
       }# end of piek == 1 check
      ##################################################################################
      if(length(piek$peaks[1,]) >= 3){ #give a warning if we find more than 2 populations
         #mark warning
         badluck[3] <<- TRUE
       }# END if peaks >= 3
      ##################################################################################
      # we've handled length(piek$peaks[1,]) of 1 and warned if 3-or-more
      # now make the initial piik matrix for the two main populations.
      ##################################################################################
      if(skip == 0){
       #set the number of populations
        pops <<- length(piek$peaks[1,])
       # set the peak limits using valleys defined by peakfinder , row 1: limits low FU peak, row 2: limits high FU peak
        pea.lim <<- matrix(krn[2, c(head(piek$valleys[1,], 2), tail(piek$valleys[1,], 2))], nrow = 2, byrow = T) # in FU
       #make the 'rough guess' population matrix
        pjik <- matrix(c(krn[2, piek$peaks[1, range(seq_along(piek$peaks[1,]))]], # FU of first and last peak
                         pea.lim[1,1], pea.lim[2,2], pea.lim[1,2], pea.lim[2,1]),ncol = 2, byrow = T)
        # for clarity we add row names:
        rownames(pjik) <- c("x.FU",        # first row is population center position in Fluorescence units
                            "outer.bounds",# second row is the outer limits (lowest and hightst FU over ALL)
                            "inner.bounds")# third row is the inner population bounds (between the populations)
         # I know the order of the values in piik is not completely logical & I have no idea why i did it this way
         # I suppose I had a valid reason when I first wrote this & I'll keep it this way for compatibility's sake
         }
     ##############################################################################
     }else{ # Control LVL 1
          # to arrive here, rough.cut must have failed. All is lost, throw an error and set failure to 1.
          badluck[1] <<- TRUE
          pjik <- matrix(NA, ncol = 2, nrow = 3)
          rownames(pjik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!
          fail <- 1 #make sure to notify failure
     }##Control LVL 1 END
     return(pjik)
    }# ROUGH.CUT END


##################################################################################
#Internal function: CLOUDY.THRESH(finds exact population boundary limits)
##################################################################################
    cloudy.thresh <- function(drops, pjjk, threshold, pii.lim, trim = T){
      # continous the work rough.cut has started
      ##################################################################################
       #*iterative* population delimiter
       ################################
       # allright, using this initial delimination of the clouds we can calculate more accurate values
       # by iterating the update process three times (i.e. use current estimate to calculate pos and neg,
       # then recaculte sd and re-delimite bands). Three iterations is generaly enough to reach
       # stability in the first decimal place.
         for ( i  in 1 : 3 ) {
          #split in positive and negatives
           pdrp <- drops[drops >= pjjk[3, 2] & drops <= pjjk[2, 2]]
           ndrp <- drops[drops >= pjjk[2, 1] & drops <= pjjk[3, 1]]
           #update s multiplier
           s1 <- s.funk(abs(moments(ndrp)[4]))
           s2 <- s.funk(abs(moments(pdrp)[4]))
           #back-up
           pjjk.bup <- pjjk
           #update
           pjjk[2, 1] <-  pjjk[1, 1] - mad(ndrp) * s1 #leftmost extreme (start of neg band)
           pjjk[2, 2] <-  pjjk[1, 2] + mad(pdrp) * s2 #rigthmost extreme (end of pos band)
           pjjk[3, 1] <-  pjjk[1, 1] + mad(ndrp) * s1 #first inner extreme (end of the neg band)
           pjjk[3, 2] <-  pjjk[1, 2] - mad(pdrp) * s2 #second inner extreme (start of the pos band)
           # check if we didn't screw up (which may happen with very few pos/neg)
           if(any(is.na(pjjk))){
             #maybe a warning?
             pjjk <- pjjk.bup # error error, restore
             break # stop optimization
           }# END NA check
          }#END of the iterative update
       ################################
       #trimming: adjust for asymetric populations
       #############################################
       # added a switch to turn of trimming, this is necessary if we intend to run GEV after cloudy, removing 'bona fide' outliers will affect gev
       # anyway, the trimming is just to get nicer plots, it does not affect the threshold placement of cloudy
       if(isTRUE(trim)){
        ## update populations
         pdrp <- drops[drops >= pjjk[3, 2] & drops <= pjjk[2, 2]]
         ndrp <- drops[drops >= pjjk[2, 1] & drops <= pjjk[3, 1]]
        ## we will only attempt trimming if both populations have more than 10 droplets, otherwise hist may crash
        if(length(pdrp) > 10 & length(ndrp) > 10){
        ## use hist to break the population into breaks so we can find trailing empty breaks
         # use breaks of approx 1/10 th of the median average deviation
         nhist <- hist(ndrp, breaks = floor(diff(range(ndrp))/(mad(ndrp)*0.1)), plot = F)
         phist <- hist(pdrp, breaks = floor(diff(range(pdrp))/(mad(pdrp)*0.1)), plot = F)
         # retrieve counts
         ncount <- nhist$counts # retrieve counts
         pcount <- phist$counts
         ncount[ncount == 1] <- 0 #also set single droplet breaks to 'empty'
         pcount[pcount == 1] <- 0
         ncount <- rle(ncount) # replace ncount by it's run length encoding
         pcount <- rle(pcount) # replace pcount by it's run length encoding
        ## we can now inspect the first and last element of the rle and see if it has value zero, if so we can subtrace its length * break width from the relevant piik entry
         # check for empty bins on the lower population bound
          if(ncount$values[1] == 0){ pjjk[2, 1] <- pjjk[2, 1] + diff(nhist$breaks)[1] * ncount$length[1] } #
          if(pcount$values[1] == 0){ pjjk[3, 2] <- pjjk[3, 2] + diff(phist$breaks)[1] * pcount$length[1]} #
         # check for empty bins on the higher population bound
          if(tail(ncount$values, n=1) == 0){ pjjk[3, 1] <- pjjk[3, 1] - diff(nhist$breaks)[1] * tail(ncount$length, n=1) } #
          if(tail(pcount$values, n=1) == 0){ pjjk[2, 2] <- pjjk[2, 2] - diff(phist$breaks)[1] * tail(pcount$length, n=1) } #
         # correct threshold for trimming
          tresh.corr <- diff(nhist$breaks)[1] * tail(ncount$length, n=1)
        ## clean up superfluous variables
        rm(ncount, nhist, pcount, phist)
        }#
       }# END of Trimming section
       ################################
       #calculate threshold
       #############################################
         #If there is a threshold specified we'll use that, otherwise we calculate our own
         if(is.na(threshold)){
           #We now calculate the threshold for counting positives
            tresh.cld <<- pjjk[1, 1] + (1.5 * s1) * mad(ndrp) # - tresh.corr
           #should the threshold fail, we provide a backup
            tresh.bup1 <<- pii.lim[1, 2]
         }else{
            tresh.cld <<- threshold
         }#END if is.na(threshold)
      # a last check before we return the output:
      ##############################################################################
      ###Control LVL 2 : a problem related to standard deviations is overlapping bounderies
      #check : upper neg boundary is situated above lower pos boundary?
      if(pjjk[3, 1] > pjjk[3, 2]){##Control LVL 5
          badluck[4] <<- TRUE
          #we can use the 'valleys' returned by findpeaks to draw new population bounderies
          pjjk[3, 1] <- pii.lim[1, 2]
          pjjk[3, 2] <- pii.lim[2, 1]
      }##Control LVL 2 END
      ################################
      # END of the *cloudy* population delimiter
       return(pjjk)
    }# CLOUDY.THRESH END

##################################################################################
#Internal function: gev.THRESH(finds threshold using general extreme value theory)
##################################################################################

gev.thresh<- function(sdrp, cutoff.quantile = 0.9995, reps = 10, blocks = 150){
     # as in the original ddpcrquant we do the whole block sampling thing 10 times, we'll then average the threshold we got
  ##############################################################################
  # check if we have enough negative droplets to actually do the analysis (if not we return NA)
  ##########################################################
  if(length(sdrp) >= blocks*10){ #demand at least 10 values per block
    #initialize vectors
    quantgev <- rep(NA, times = reps)
    #start repeated calculations
    #############################
      for(k in 1:reps){
        ### block creation
        sdrp <- sample(sdrp)   # randomize droplet order
        size <- floor(length(sdrp)/blocks) #block size
        the.blocks <- matrix(sdrp[1:(blocks * size)], ncol =  blocks, nrow = size)
        signal.maxima <- apply(the.blocks, 2, max, na.rm = T) # function(x) sort(x, decreasing = T)[2] #function allows to select nth max value
        ## fit the GEV model using ML
        droplet.fit <- try(fgev(signal.maxima, warn.inf = FALSE), silent = TRUE)
        thresh.k <- try(qgev(cutoff.quantile, droplet.fit$estimate[1],
                             droplet.fit$estimate[2], droplet.fit$estimate[3]), silent =TRUE)
        #catch some errors
        quantgev[k] <- if(inherits(thresh.k, "try-error")|(thresh.k > max(sdrp)+3000)){ NA }else{ thresh.k }
      } # END for k
    g.threshold <- mean(quantgev, na.rm = T)
  }else{ # not enough droplets to perform analysis
    g.threshold <- NA
    badluck[8] <- TRUE
  } # END number of droplets check
  return(g.threshold)
  }# GEV.THRESH END

################################################################################
#messaging system
     public.address <- c("WARNING: something is not right, please check data \n",    #1
                         "WARNING: only one population detected, returning NA",      #2
                         "WARNING: multiple population detected \n",                 #3
                         "WARNING!: population boundary overlap \n   boundaries rewritten, not based on population parameters \n", #4
                         "Rotation seems to have screwed up! ...restoring droplets... \n", #5
                         "WARNING!: threshold placement failed \n   Check data, arbitrary threshold was set \n", #6
                         "WARNING!: threshold inside positive cloud \n   Forced lower threshold, not based on population parameters \n", #7
                         "Not enough negative droplets to calculate GEV threshold \n"  #8
                       )##
##################################################################################
## ACTUAL FUNCTION START #########################################################
##################################################################################
   # Choose your adventure:
    if(method == "simplex"){ # simplex analysis of channel 1, no channel 2 needed
      drp[[2]] <- NULL # na.omit(drp$Ch1)
    }else if(method == "simplex2"){ # simplex analysis of channel 2
      drp[[1]] <- NULL #na.omit(drp$Ch2)
    }#else if(method %in% c("eva", "duplex")){ # analysis of channel 1 with crosstalk correction using channel 2
   # Choose your threshold:
    if(is.numeric(threshold)){ #because threshold used to only two possibilities (NA or numeric) we'll have to split it in
                               # 2 variables so that the code below keeps making sence (multiple calls to is.na(threshold) )
     #skip
    }else{
     the.tresh <- threshold
     threshold <- NA
    }# if else END
   ###############################################
   # allocate space to store loop outcome
   rown <- if(method != 'eva'){length(drp)}else{1} # eva uses both channels but has only one output
   the.rows <- c('Ch1', 'Ch2', 'Ch3', 'Ch4', 'Ch5', 'Ch6') # so yeah, we only support two colors for now, but it let's try and be ambitious
   if(method == "simplex2"){the.rows <- the.rows[-1]} # make sure the output gets labeled "Ch2" when only inspecting Ch2
   #write all NA vectors:
    targall <- matrix(NA, nrow = rown, ncol = 3)
      colnames(targall) <- c("targ.in.sVol", "upper", "lower")
      rownames(targall) <- the.rows[1:rown]
    lambs  <- matrix(NA, nrow = rown, ncol = 3)
      colnames(lambs) <- c("lambda", "upper", "lower")
      rownames(lambs) <- the.rows[1:rown]
    perform <- matrix(NA, nrow = rown, ncol = 3)
      colnames(perform) <- c("resolution", "%rain", "%comparted")
      rownames(perform) <- the.rows[1:rown]
    drops <- matrix(NA, nrow = rown, ncol = 3)
      colnames(drops) <- c("positive", "negative", "rain")
      rownames(drops) <- the.rows[1:rown]
    popsi <- matrix(NA, nrow = rown, ncol = 1)
      colnames(popsi) <- c("populations")
      rownames(popsi) <- the.rows[1:rown]
    treshes <- matrix(NA, nrow = rown, ncol = 1)
      colnames(treshes) <- c("threshold")
      rownames(treshes) <- the.rows[1:rown]
    # cleanup
    rm(rown)

  ## LOOP OVER CHANNELS #########################
  ###############################################
  for(d in 1:length(drp)){
   # print the channel name so we know what channel the messages relate to.
    if(!isTRUE(silent)){message(paste("Analysis channel", the.rows[d], sep = " "))}
   # channels select
   drp1 <- na.omit(drp[[d]])
   # in case rotation will be tried, we'll need another channel, so let's get that data
   if(method %in% c("eva", "duplex")){
          drp2 <- if(d == 1){na.omit(drp[[d + 1]])}else{na.omit(drp[[d - 1]])}
          # this should ensure there is a drp2 data available for 'rotation' analysis of Ch2
          # note that I have currenlty disabled rotation for duplex as I cannot guarantee error-free results for now
          }#

   # create some values for the internal function to update (<<-)
     fail <- 0 #create failure tracker
     pops <- 0 #create population tracker
     skip <- 0 #create skip tracker
     pea.lim <- 0 #needed as backup in case things go wrong
     tresh.cld      <- NA # cloudy threshold
     tresh.bup1 <- NA
     badluck <- rep(FALSE, times = 8)

   ### call to rough.cut function
    piik <- rough.cut(drp1, silent, threshold)

   #### #check if we should continue or give up#####################################
     if(fail == 0 & skip == 0){
       ## call to treshfinder functions
        # before we call cloudy and over-write the rough estimate, use it to get the GEV-threshold
        # unless we use the 'mixed' threshold setting, then we use cloudy's population delimitation to feed into gev
         if(the.tresh != "mixed"){
                sdrp <- subset(drp1, drp1 > piik[2, 1] & drp1 < piik[3, 1])       # rough guess of which droplets are negative
                piik <-   cloudy.thresh(drp1, piik, threshold, pea.lim, trim = T) # trim is on!
                }#
         if(the.tresh == "mixed"){
                piik <-   cloudy.thresh(drp1, piik, threshold, pea.lim, trim = F) # turn off trim!
                sdrp <- subset(drp1, drp1 > piik[2, 1] & drp1 < piik[3, 1])       # better guess of which droplets are negative
                } #
         tresh.gev <- gev.thresh(sdrp)
      ##############################################################################
      }else if(fail == 0){# if fail is zero there still is a valid guess for the populations so we can call gev.thresh
        sdrp <- subset(drp1, drp1 > piik[2, 1] & drp1 < piik[3, 1] )  # rough guess of which droplets are negative
        tresh.gev <- gev.thresh(sdrp)
      }#############################################################################

   # use piik to define the populations
     pospop <- drp1 >= piik[3, 2] & drp1 <= piik[2, 2]
     negpop <- drp1 >= piik[2, 1] & drp1 <= piik[3, 1]

##################################################################################
# Rotation!
##################################################################################
# note that at this point the initial call to threshfinder may have failed and maybe fail == 1
# However, if rotate == T we may still have a chance at making an analysis by improving contrast.
# We'll simply use all droplets in the regression (since they are so close together that they appear undescernable anyway)
    if(method %in% c("eva")){    # I have removed the rotation for duplex analysis for now. I do not have the right data to properly troubleshoot it
                                 #  & calling duplex on eva data results in undesired outcomes, so for now I limit the use to 'eva' data where it will always work properly
      # set skip back to zero
      skip <- 0 # if it was one it would mess with the function call to rough.cut
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
      ##############################################################################
      ### another call to the threshold functions!
      bad.bup <- badluck # backup
      badluck[1:4] <- rep(FALSE, times = 4) #reset badluck
      ### again: first call rough.cut (FU values will have shifted, the original piik makes no sense anymore)
       piik <- rough.cut(drp1, silent, threshold)
       #### #check if we should continue or give up#####################################
        if(fail == 0 & skip == 0){
          ## call to treshfinder functions
             # before we call cloudy and over-write the rough estimate, use it to get the GEV-threshold
             # unless we use the 'mixed' threshold setting, then we use cloudy's population delimitation to feed into gev
               if(the.tresh != "mixed"){
                      sdrp <- subset(drp1, drp1 > piik[2, 1] & drp1 < piik[3, 1])       # rough guess of which droplets are negative
                      piik <-   cloudy.thresh(drp1, piik, threshold, pea.lim, trim = T) # trim is on!
                      }#
               if(the.tresh == "mixed"){
                      piik <-   cloudy.thresh(drp1, piik, threshold, pea.lim, trim = F) # turn off trim!
                      sdrp <- subset(drp1, drp1 > piik[2, 1] & drp1 < piik[3, 1])       # better guess of which droplets are negative
                      } #
               tresh.gev <- gev.thresh(sdrp)
         ##############################################################################
         }else if(fail == 0){# if fail is zero there still is a valid guess for the populations so we can call gev.thresh
           sdrp <- subset(drp1, drp1 > piik[2, 1] & drp1 < piik[3, 1] )  # rough guess of which droplets are negative
           tresh.gev <- gev.thresh(sdrp)
         }#############################################################################
      #################################################################################
      # on the off chance that the rotation totally screwed up the data we need to undo all the changes!
        if(!any(is.na(piik.bup)) & fail == 1){ ## if that is the case, fail will be 1, in all other cases a valid piik will be produced
                                                # actually, that's not true, threshfinder may have failed before, so we only restore if piik.bup
                                                # does not contain any NA. Otherwise we keep fail == 1 & give up
           # restore variables
           badluck <- bad.bup
           badluck[5] <- TRUE # let user know we screwed up
           drp1 <- drp1.bup
           piik <- piik.bup
           pops <- pops.bup
           neg.ref <- neg.ref.bup
           rm(piik.bup, drp1.bup, pops.bup, bad.bup, neg.ref.bup) # cleanup
           # reset fail
           fail <- 0
        }# END fail check

      # use the new and improved piik to re-define the populations
      pospop <- drp1 >= piik[3, 2] & drp1 <= piik[2, 2]
      negpop <- drp1 >= piik[2, 1] & drp1 <= piik[3, 1]
      # show a quick plot of before and after rotation (IF rotation was succesful):
         if(isTRUE(plots) & exists("drp1.bup")){
            plot1 <- ggplot(data.frame(ch1 = drp1.bup, ch2 = drp2), aes(x = ch2, y = ch1)) + geom_point(color = "#7BA8C2") + labs(title = "data")
            plot2 <- ggplot(data.frame(ch1 = drp1, ch2 = drp2), aes(x = ch2, y = ch1)) + geom_point(color = "#7BA8C2") + labs(title = "corrected")
            x11(width = 14, height = 6)
            grid.arrange(plot1, plot2, nrow = 1)
         }# END plots
    }#END if isTRUE rotator

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
   if(is.na(tresh.cld)){##Control LVL 3
        badluck[6] <- TRUE
        tresh.cld <-  tresh.bup1
        ndrp <- drp1[drp1 < tresh.cld]
        pdrp <- drp1[drp1 > tresh.cld]
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
    #we build in this simple check to make sure the threshold is not placed into the positive band
    if(tresh.cld > piik[3, 2]){ ##Control LVL 4
       tresh.cld <-  min(tresh.bup1, piik[3, 2])  ### in fact there is no guarantuee that tresh.bup1 is outside the positive cloud.
                                 # therefore, we'll take either tresh.bup1 or the lower bound on the pos cloud, whichever is lower
    }##Control LVL 4 END
    }# end if(is.na(threshold))

##################################################################################
#we can now start the trivial calculations:
##################################################################################
    #Pick threshold, if mixed or gev but GEV failed than we will use cloudy anaway
    if(is.na(tresh.gev) & the.tresh %in% c("gev", "mixed")){
       if(!isTRUE(silent)){message("gev threshold failed, using cloudy threshold")}
       tresh <- tresh.cld
    }else{
       tresh <- c(tresh.cld, tresh.gev, mean(c(tresh.cld, tresh.gev)))[the.tresh == c("cloudy", "gev", "mixed")]
    }# if else END
    #count the raindrops (i.e. the population of the 'clearance band' between the positive and negative band)
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
       dfdrops[dfdrops$FU > piik[2, 2], 'category']  <- 4  #drift is above the pos band (high-Fluorescence outliers)
       dfdrops[dfdrops$FU < piik[2, 1], 'category']  <- 4  #drift is also below the neg band (low-Fluorescence outliers)
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


  ############################
  ### Save to ouput matrices
  ############################

    lambs[d, ]   <- c(lambda, lamhi, lamlo) # lambda estimates
   # targets[d, ] <- tdrp * lambs  #total number of targets present in droplets analysed
    targall[d, ] <- lambs[d, ] * (sVol/(dVol * 0.001))  #total number of targets in sample volume
    drops[d, ]   <- c(sum(dfdrops$positive), sum(!dfdrops$positive), rain)
    perform[d, ] <- c(Resol, pRain, Comp)
    treshes[d, ] <- tresh
    popsi[d, ] <- pops

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
      # for graphical purposes we limit the upper bound of the bandiwth otherwise things can get over-smoothed
       temp <- bw.nrd0(drp1)
       if(temp < 50){temp <- 50}else if(temp > diff(range(drp1))/40){temp <- floor(diff(range(drp1))/40)}
     dens <- density(drp1, bw = temp)
     plotdf <- data.frame(x=dens$x, y=dens$y)
     plotdf$dropsplit <- factor(findInterval(plotdf$x, tresh))

     #1# kernel density
      plot1 <- ggplot(plotdf, aes(x,y)) + geom_ribbon(aes(ymin=0, ymax=y, fill=dropsplit)) +  geom_line() +
               scale_fill_manual(labels = c("negative", "positive"), limits = factor(c(0, 1)), values = c('#F8766D', '#00BFC4'), name = "Population") +
               labs(tag = "A", x = "FU", y = "density")
      if(fail == 0){ #piik lines can only be added if piik is not NA
               plot1 <- plot1 + geom_vline(data = dfpiik, aes(xintercept = x.FU), linetype = "dashed") +
                     geom_vline(data = dfpiik, aes(xintercept = outer.bounds), linetype = "dashed") +
                     geom_vline(data = dfpiik, aes(xintercept = inner.bounds), linetype = "dashed") +
                     geom_vline(xintercept = tresh, size = 0.5)
               if(!is.na(tresh.gev)){plot1 <- plot1 + geom_vline(xintercept = tresh.gev, linetype = "dotted")} # only add GEV threshold if it's not NA
      }#
     #2# droplet view
      plot2 <- ggplot() + geom_point(data = dfdrops, aes(y = FU, x = ordr, colour = category)) +
               scale_colour_gradientn(breaks = c(1, 2, 3, 4), colors = c("#F8766D", "#FFC425", "#00BFC4", "#0076C4"), labels = c("negative", "rain", "positive", "outlier"), limits = c(0.9,4.1), name = "Population") +
               geom_hline(yintercept = tresh.cld, linetype = "dotted") +
               geom_hline(yintercept = tresh, size = 0.5) +
               labs(tag = "B", x = "Observation")
               if(!is.na(tresh.gev)){plot2 <- plot2 + geom_hline(yintercept = tresh.gev, linetype = "dotted")} # only add GEV threshold if it's not NA
     ### put on screen
      x11(width = 14, height = 6)
      grid.arrange(plot1, plot2, nrow = 1)
  }# END plots

##################################################################################
#Standard INFO messages, comment out for now they are a distraction when analysing multiple data
##################################################################################
  #  if(!isTRUE(silent)){
  #    if(isTRUE(Resol < 2.5)){  message("Resolution: FAILED")}else{message("Resolution: PASSED")} #isTRUE wrapper as Resol can be NA
  #    if(is.na(pRain)){message("Rain: multiple populations, rain not calculated")
  #                    }else if(pRain < 0.025){message("Rain: PASSED")
  #                    }else{message("Rain: FAILED")}
  #    if(Comp < 0.3){message("Compartmentalization: FAILED")}else{
  #    if(Comp < 0.5){message("Comparmentalization: MEDIOCRE")}else{
  #                   message("Comparmentalization: PASSED")}}
  #    if(lambs[d, 1] < 0.01){message("Confidence bounds: NOT optimal")}else{
  #    if(lambs[d, 1] < 5){message("Confidence bounds: OPTIMAL")}else{
  #                   message("Comparmentalization: NOT optimal")}}
  #   message("\n --------- \n")
  #   }#end of messages
##################################################################################
#ERROR messages
##################################################################################
   }else{
    message(paste("WARNING: threshold finding failed for Ch", d, ", returning NA", sep = ""))
   }##Control LVL 2 END

 # in case of method 'eva' there are 2 channels but we should only loop once, so we break here
 if(method == 'eva'){break}
 }## END OF LOOP OVER CHANNELS ###################################################

}else{ #consistency check failed: throw some errors
    message("WARNING: input inconsistent, returning NA")
    error.log <- c("You may only choose 1 method, argument should be of length 1 \n Available methods: simplex, simplex2, eva, and duplex. Other choices are invalid \n",
                   "dVol & sVol should both be of length 1 \n",
                   "dVol & sVol should both be a numerical \n",
                   "plots, silent, and vec all must be a logical \n",
                   "plots, silent, and vec all should be of length 1 \n",
                   "threshold should be either numerical or one of the available methods ('cloudy', 'gev', 'mixed' \n",
                   "threshold should be of length 1, you cannot set multiple thresholds in the same function call")
    message(error.log[checks])
}## consistency check
## report messages
if(!isTRUE(silent)){
  if(any(badluck)){message("\n --------- \n")}
  message(public.address[badluck])
  }##
##################################################################################
# actual output ##################################
##################################################################################
    if(isTRUE(vec)){ #output if
        return(cbind(targall, lambs, perform, drops, popsi, treshes))
    }else{
        return(list("targets.in.sVol" = targall, "lambda" = lambs, "performance" = perform, "droplets" = drops, "populations" = popsi, "threshold" = treshes))
    }#output if END
##################################################################################
 }## END of recursion
}##END of 'cloudy'
##################################################################################
##################################################################################
##################################################################################
