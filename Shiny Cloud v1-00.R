library(shiny)
require(SuppDists)#required for the 'moments()' function
require(MASS)     #required for the 'lqs()' function
library(gridExtra)
require(ggplot2)

########################################################
###### Function repository##############################
########################################################

###### plotting for Shiny ##############################
########################################################

ddplot <- function(drop, ...){
  drp <- na.omit(drop)
  temp <- sample(c(1:length(drop)), size = length(drop), replace = F)
  plot(drop ~ temp, col = rgb(123, 168, 194, maxColorValue = 255), bty = "L", pch = 20 , xlab = "observation", ylab = "FU", ylim = c(min(drop), max(drop)), ...)
  }# END of ddplot ########################################################

poplotter <- function(drop, piek, info, manthre, ...){
  fail <- info[1]
  pops <- info[2]
  tresh <- info[3]
  # check for threshold consistency
  if(is.na(tresh) & !is.na(info[4])){
    tresh <- info[4]
  }else if(is.na(tresh) & is.na(info[4])){
    tresh <- 0 # for shiny to produce output, tresh cannot be NA
  }
  # check for manual threshold
  if(!is.na(manthre)){
    tresh <- manthre
  }
  
  drp1 <- na.omit(drop)
  dfdrops <- data.frame(FU = drp1, positive = drp1 > tresh, category = 0)
  tdrp <- length(drp1) # total number of droplets
  
  # the following should only be calculated if there is more than 1 population
  if(pops > 1){
    #also define 'cerntainly pos', 'rain' and 'drift'
    dfdrops[!dfdrops$positive, 'category'] <- 1  # negative droplets
    dfdrops[dfdrops$FU > tresh & dfdrops$FU < piek[3, 2], 'category'] <- 2    #rain is the population of the clearance band (counted positive)
    dfdrops[dfdrops$FU > piek[3, 2] & dfdrops$FU < piek[2, 2], 'category']  <- 3 #c-pos is the population of the pos band only
    dfdrops[dfdrops$FU > piek[2, 2], 'category']  <- 4  #drift is population above the pos band (high-Fluorescence outliers)
  }else{ # in case of only 1 population Resolution & so on should still be defined
    dfdrops[!dfdrops$positive, 'category'] <- 1  # negative droplets
    dfdrops[dfdrops$positive, 'category'] <- 3  # positive droplets
  }# END pops > 1
  
  #piik thing
  dfpiik <- as.data.frame(t(piek))
  dfpiik$category <- c(1,3)
  
  # we also have to add a ranomd order to dfdrops to plot them randomly
  dfdrops$ordr <- sample(c(1:nrow(dfdrops)), size = nrow(dfdrops), replace = F)

  # unfortunately our plan for using geom_density does not work out well, both populations get scaled differently which gives a very misleading plot
  # so we must make a separate plotting data frame:
  dens <- density(drp1)
  plotdf <- data.frame(x=dens$x, y=dens$y)
  plotdf$dropsplit <- factor(findInterval(plotdf$x, tresh))
    
  #1# kernel density
  plot1 <- ggplot(plotdf, aes(x,y)) + geom_ribbon(aes(ymin=0, ymax=y, fill=dropsplit)) +  geom_line() +
    scale_fill_manual(labels = c("negative", "positive"), limits = c(0, 1), values = c('#F8766D', '#00BFC4'), name = "Population") +
    labs(tag = "A", x = "FU", y = "density")
  if(fail == 0){ #piik lines can only be added if piik is not NA
    plot1 <- plot1 + geom_vline(data = dfpiik, aes(xintercept = x.FU), linetype = "dashed") +
      geom_vline(data = dfpiik, aes(xintercept = outer.bounds), linetype = "dashed") +
      geom_vline(data = dfpiik, aes(xintercept = inner.bounds), linetype = "dashed")
  }
  
  #2# droplet view
  plot2 <- ggplot() + geom_point(data = dfdrops, aes(y = FU, x = ordr, colour = category)) +
    scale_colour_gradientn (breaks = c(1, 2, 3, 4), colors = c("#F8766D", "#FFC425", "#00BFC4", "#0076C4"), labels = c("negative", "rain", "positive", "outlier"), limits = c(0.9,4.1), name = "Population") +
    geom_hline(yintercept = tresh, linetype = "dashed") +
    labs(tag = "B", x = "Observation")
  
  ### put on screen
  grid.arrange(plot1, plot2, nrow = 1)
}# END of poplotter ########################################################

###### Adapted internal functions ######################
########################################################

s.funk <- function(k){return(3.8 + 0.35 * log(k) + 0.045 * log(k)^2 +0.75)} # kurtosis to spread

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

findthresh <- function(drops, silent, threshold, neg.ref,...){  #
  skip <- 0
  # create some values for the internal function to update (have to be internally defined for shiny!)
  fail <- 0 #create failure tracker
  pops <- 0 #create population tracker
  tresh      <- NA
  tresh.bup1 <- NA
  badluck <- rep(FALSE, times = 7)
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
        pops <- 1
        tresh <- threshold
        # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
        piik <- matrix(c(threshold - 100, threshold + 100, min(drops) - 100, max(drops) + 100, threshold -0.1, threshold + 0.1),
                       nrow = 3, ncol = 2, byrow = T)
        rownames(piik) <- c("x.FU", "outer.bounds", "inner.bounds") #  the plotting dataframe needs names!        
        # we should, however, not attempt further analysis of the data or we screw up our ersatz piik
        skip <- 1
        ### we have no outliers, but at least the neg.ref may allows us to see if the drops are all pos or neg
        #############################
      }else if(!is.na(neg.ref)){
        pops <- 1
        if(median(drops) > neg.ref){ # ALL positive, we need at least 1 negative
          tresh <- mean(tail(sort(drops, decreasing = T), 2)) # put threshold between two lowest droplets
          # we should construct a makeshift piik matrix so that we can still return an output (and produce some plots)
          piik <- matrix(c(min(drops), median(drops), min(drops) - 100, max(drops) + 100, tresh -0.1, tresh + 0.1),
                         nrow = 3, ncol = 2, byrow = T)
        }else{ # all negative, zero positives is perfectly ok, we set the threshold above the max droplet
          tresh <- max(drops) + 2 * mad(drops)
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
        pops <- 1
        fail <- 1 # notify failure
        skip <- 1 # make sure we skip
        badluck[2] <- TRUE
      } #end of no outlier case
    }# end of piek == 1 check
    ## #check if we should continue or give up#########################################
    if(fail == 0 & skip == 0){#####################################################################
      if(length(piek$peaks[1,]) >= 3){ #give a warning if we find more than 2 populations
        #mark warning
        badluck[3] <- TRUE
      }# END if peaks >= 3
      #*universal* population delimiter
      ################################
      #set the number of populations
      pops <- length(piek$peaks[1,])
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
        tresh <- piik[1, 1] + (s1 + 1/2 * s1) * mad(ndrp)
        #should the threshold fail, we provide a backup
        tresh.bup1 <- pea.lim[1, 2]
      }else{
        tresh <- threshold
      }#END if is.na(threshold)
      # a last check before we return the output:
      ##############################################################################
      ###Control LVL 2 : a problem related to standard deviations is overlapping bounderies
      #check : upper neg boundary is situated above lower pos boundary?
      if(piik[3, 1] > piik[3, 2]){##Control LVL 5
        badluck[4] <- TRUE
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
    badluck[1] <- TRUE
    piik <- matrix(NA, ncol = 2, nrow = 3)
    fail <- 1 #notify failure
  }##Control LVL 1 END
  
  ### extra bit for shiny, as this code cannot fit in valkyr
  if(fail == 0){ ##Control LVL 2
    ##################################################################################
    ##third to fifth control points: checks for stability reasons
    ##################################################################################
    #in case of failure (tresh is NA) we fall back on the backup threshold
    if(is.na(tresh)){##Control LVL 3
      badluck[6] <- TRUE
      tresh <-  tresh.bup1
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
  }##Control LVL 2 END
  ##############################################################################
  sumar <- c(fail, pops, tresh, tresh.bup1)
  names(sumar) <- c("fail", "pops", "tresh", "tbup")
  return(list(piik = piik, infos = sumar, bad = badluck))
}#End of findtresh function########################################################

valkyr <- function(zealot, drops, dVol, sVol, manthre,... ){# function that calculates the values for output from the treshold set by treshfinder
  # split zealot into its components so we can re-use the same code as in the original cloudy algortihm
  # (it's not lazy, it will make future updates more straightforward, we can just copy paste everything!)
  drp1 <- na.omit(drops)
  piik <- zealot$piik
  fail <- zealot$infos[1]
  pops <- zealot$infos[2]
  # we should only continue if fail == 0, otherwise the entire script will fail
  if(fail == 0){
    tresh      <- zealot$infos[3]
    tresh.bup1 <- zealot$infos[4]
    # check for threshold consistency
      if(is.na(tresh) & !is.na(zealot$infos[4])){
        tresh <- zealot$infos[4]
      }else if(is.na(tresh) & is.na(zealot$infos[4])){
        tresh <- 0 # for shiny to produce output, tresh cannot be NA
      }#
      # check for manual threshold
      if(!is.na(manthre)){
        tresh <- manthre
      }#
    
    badluck    <- zealot$bad
    # use piik to define the populations
    pospop <- drp1 >= piik[3, 2] & drp1 <= piik[2, 2]
    negpop <- drp1 >= piik[2, 1] & drp1 <= piik[3, 1]
    pdrp <- drp1[pospop]
    ndrp <- drp1[negpop]
    ## beyond this point is the 'copy-paste' of cloudy code ######################
    
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
    # actual output ##################################
    ##################################################################################
    return(c(targall, lambs, perform, drops, pops, tresh))
  }else{ # fail == 1, we return a lot of NA & zeroes
    return(c(rep(NA, times = 6), 0, NA, length(drp1) * dVol * 0.001 / sVol, length(drp1), 0, 0, pops, 0))
  }##
  #return(list("targets.in.sVol" = targall, "lambda" = lambs, "performance" = perform, "droplets" = drops, "populations" = pops, "threshold" = tresh))
}#End of valkyr function########################################################

# pucblic adress messages
public.address <- c("WARNING: something is not right, please check data <br/>",    #1
                    "WARNING: only one population detected, returning NA <br/>",      #2
                    "WARNING: multiple population detected <br/>",                 #3
                    "WARNING!: population boundary overlap. The boundaries were rewritten and are not based on population parameters <br/>", #4
                    "Rotation seems to have screwed up! ...restoring droplets... <br/>", #5
                    "WARNING!: threshold placement failed. Check the data, an arbitrary threshold was set <br/>", #6
                    "WARNING!: threshold inside positive cloud. Algortihm forced a lower threshol which is not based on population parameters <br/>" #7
)##


########################################################
###### App Code ########################################
########################################################

# Define UI for data upload app ----
ui <- fluidPage(
  titlePanel("welcome to Cloudy"),
#### tab 1
  #tabPanel("Uploading Files",
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        ### Sidebar panel for inputs ----
        sidebarPanel(
          # Input: Select a file ----
          fileInput("file1", "Choose CSV File",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          ## Horizontal line ----
          tags$hr(),
          # Input: Select Channel ----
          numericInput(inputId = "Ch", label = "Channel:",
                       value = 1, min = 1, max = 2),
          ## Horizontal line ----
          tags$hr(),
          # Input: Select variables ----
          sliderInput(inputId = "dVol", label = "droplet volume (nl):",
                       value = 0.85, min = 0.1, max = 1),
          sliderInput(inputId = "sVol", label = "sample volume (ul):",
                      value = 20, min = 5, max = 25),
          ## Horizontal line ----
          tags$hr(),
          numericInput(inputId = "thresh", label = "manual threshold:",
                       value = NA, min = 1000, max = 25000),
          #numericInput(inputId = "Nref", label = "neg.ref:",
          #             value = NA, min = 5000, max = 25000),
        width = 3),
        ### Main panel for displaying outputs ----
        mainPanel(
          # Output: Data file ----
          tabsetPanel(
            #### tab 1
            tabPanel("help", 
                     h3("getting started"),
                     p("Use the 'browse' button to select a *.csv file containing digital droplet data. The algorithm is designed for use with data coming from the Biorad QX200 platform. It will probably choke on any other format."),
                     p("The 'inspect data' tab allows you to preview the data. You can use it to select the correct 'Channel' (either 1 or 2) if you are unsure which one contains the correct values"),
                     p("The other tabs contain the algorithm output. Please be patient as the plots may take a while to load. The 'plot' tab will also display some warning messages (if appropriate, most of the time there shouldn't be any)"),
                     p("I have made it so you can mess about with the further inputs. The first two provide the partition volume (in nanoliter) and the sample volume (in microliter). Changing these will affect the concentrations calculated, but will not affect the threshold setting or the calculation of lambda."),
                     p("Lastly, there is 'manual threshold'. If provided, a manual threshold will override the one calculated by the algorithm, this will affect a lot of the output values, It's quite possible to break the plots if you push it to silly values. Leave blank if you want everything to remain fully automatic.")
            ), #, tabpanel END
            #### tab 1
            tabPanel("Inspect data", 
                     fluidRow(
                      column(5, tableOutput("contents")),
                      column(7, plotOutput(outputId = "ddplot"))
                     )), #fluidrow, tabpanel END
            #### tab 2
            tabPanel("plot", #""),
                     fluidRow(
                       column(12,plotOutput(outputId = "popplot"))
                     ), #fluidrow1 END
                     fluidRow(
                       column(12, htmlOutput(outputId = "badluck"))
                     )), #fluidrow2, tabpanel END
            #### tab 3
            tabPanel("Table", 
                     fluidRow(
                       column(12, tableOutput(outputId = "restable1"))
                     ), #fluidrow1 END
                     fluidRow(
                       column(6, tableOutput(outputId = "restable2")),
                       column(6, tableOutput(outputId = "restable3"))
                     ),#fluidrow2 END
                     fluidRow(
                       h3("A Note on Performance"),
                       p("For the performance parameters to make sense, they should be measured at a lambda of 0.7 (1:1 positive:negative partitions). Especially the resolution, which assumes peaks of equal height and approx. gaussian in shape, may be biased otherwise. Heavily tailed distributions will affect the parameters as well."),
                       p("This does not mean the measurement or analysis result is incorrect. That depends solely on the ability to distinguish positives from negatives, but the performance parameters may not be meaningful")
                     )) #fluidrow3, tabpanel END

          )# tabsetpanel END
        )# mainPanel END
      )# sidebarLayout END
  )# fluidPage END


server <- function(input, output){# server functions determines how inputs are converted to outputs
  # server function rules:
  # 1 # output objects are stored in the output list: output$something <- something
  # output elements in the ui will look into this list to find their target. eg. plotOutput("hist") will look for output$hist
  # 2 # what is saved into output should be built with a render function
  # e.g. output$hist <- renderPlot()
  # 3 # use input values with the inputId
  # e.g. input$num  after   numericInput(inputId = "num")

  # read user selected file
  df <- reactive({req(input$file1)
    read.csv(input$file1$datapath, header = T) })#
  
  # extract droplet FU readings from file according to Ch chosen by user
  drp <- reactive({req(input$file1)
    df()[, input$Ch]  })#
  
  # call to treshfinder function
  piki <- reactive({req(input$file1)
    findthresh(drp(), silent = T, threshold = input$thresh, neg.ref = NA)  })#
    
  # call valkyr to make the calculations
  res <- reactive({req(input$file1)
    valkyr(zealot = piki(), drops = drp(), dVol = input$dVol, sVol = input$sVol, manthre = input$thresh)  })#
  
  # make some preview output for the selected channel (droplet table header, droplet plot)
  output$contents <- renderTable({
    req(input$file1)
    # display file head for inspection
    head(df(), n = 10)})
  
  output$ddplot <- renderPlot({
    req(input$file1)
    ddplot(drp())
  })# END renderplot
 
  # make output plot (population plot)
  output$popplot <- renderPlot({
    req(input$file1)
    poplotter(drp(), piki()$piik, piki()$infos, input$thresh)
  })# END renderplot
  
  # make output tables
  output$restable1 <- renderTable({
    data.frame('parameter' = c("lambda", "targ.in.sVol"), 'value'= res()[c(4,1)], 'upper.limit' = res()[c(5,2)], 'lower.limit' = res()[c(6,3)])
  }, caption = "Sample measurements")# END rendertable1
  output$restable2 <- renderTable({
    data.frame('parameter' = c("threshold", "Resolution", "%rain", "%comparted", "populations"), 'value'= res()[c(14, 7:9, 13)])
  }, caption = "Performance parameters")# END rendertable2
  output$restable3 <- renderTable({
    data.frame('count' = c("total", "positive", "negative", "rain"), 'value'= c((res()[10] + res()[11]),res()[10:12]))
  }, caption = 'droplet counts')# END rendertable3
  
  # make output warning messages
  output$badluck <- renderUI({
    HTML(public.address[piki()$bad])
  })# END rendertext
  
}#END server

shinyApp(ui = ui, server = server)#