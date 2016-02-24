#########################################################################################################
## This is the algorithm for digital PCR analysis from raw compartment fluorescence measurements
## This algorithm is supplemental material to the publication "Measuring digital PCR quality: Performance Parameters and their Optimization"
## by Antoon Lievens, Sara Jacchia, Dafni Kagkli and Cristian Savini, and Maddalena Querci 

 # Updates will appear at: github.com/Gromgorgel

## ARGUMENTS ###############
## 'drp' a numeric vector of all (endpoint) fluorescence measurments in a digital reaction
 #  Readings do not have to be ordored in any particular way, NA are allowed (will be removed)
## 'dVol' = numerical, the compartment volume in nanoliter
## 'sVol' =  numerical, the sample volume in microliter
## 'plots' = logical, if set to 'TRUE' plots will be generated 
## 'silent' = logical, if is set to 'FALSE' messages will be generated
## 'vec' = logical, if set to 'TRUE' the results will be returned in a vector instead of a list

## OUTPUT ###############
## standard output is a list with the following components
 # targets.in.sVol = the estimated number of targets in the sample volume (targ.in.sVol) and its Poisson confidence bounds (upper, lower)
 # lambda  = the estimated average number of targets per compartment (lambda) and its Poisson confidence bounds (upper, lower)
 # performance = Performance paramters calculated from the input data: resolution, ratio fo rain to total compartments (%rain), degree of compartmentalization (%comparted)
 # droplets = number of compartments counted in each category (positive, negative, & rain)
 # populations = number of fluorescence populations as detected by the algorithm 
## if 'vec' is set to true, all the above parameters are returned in a single vector rather than as a list
 
## PLOTS ###############
## PLOT 1: Kernel density plot, population boundaries and peak are indicated with vertical dashed lines and are shaded 
## PLOT 2: dot plot of the flourescence readings (random order), each popultion is coloured differently, horizontal line represents the threshold set
## PLOT 3: Relative performance plot. Four parameter values are represented as barplot (and scaled to be comparable). 
           # Resolution (devided by 10), the horizontal line represents the 2.5 (0.25) acceptance limit
           # Percentage rain (multiplied by 10), the horizontal line represents the 2.5% (0.25) rejection limit
           # percentage of sample compartmentalized (horizontal lines represent the 0.3 and 0.5 limits)
           # droplet count: positive and negative droplets (stacked) scaled to the total number of droplets in the analysis

## Version history: this is stand alone version 2.04
                    #01 added a step to determine the number of populations (eg: only one => abort, two => OK, three => warn and continue)
                    #02 further finetuned to 'only one population' procedure for the cases where there are very few pos or neg droplets
                    #03 restructuring of the algorithm for stability and efficiency reasons
                    #04 added control step for overlapping bounderies & use of minima as a general backup for non-sd based boundaries

## Main function ###############                    
#########################################################################################################
cloudy <- function (drp, dVol = 0.85, sVol = 20, plots = FALSE, silent = TRUE, vec = FALSE){
   ##################################################################################
   #colorschemes for the colorblind##################################################
     two.color <- c(rgb(68,  140,  59, maxColorValue = 255), rgb(117, 170,  93, maxColorValue = 255), rgb(160, 204, 124, maxColorValue = 255),
                    rgb(219, 255, 182, maxColorValue = 255), rgb(216, 196, 213, maxColorValue = 255), rgb(192, 160, 202, maxColorValue = 255),
                    rgb(171, 124, 189, maxColorValue = 255), rgb(133,  84, 168, maxColorValue = 255))
     spectral  <- c(rgb(241,  17,  39, maxColorValue = 255), rgb(255, 119,  40, maxColorValue = 255), rgb(255, 163,  41, maxColorValue = 255),
                    rgb(253, 208,  42, maxColorValue = 255), rgb(245, 255,  46, maxColorValue = 255), rgb(139, 205, 255, maxColorValue = 255),
                    rgb(132, 132, 223, maxColorValue = 255), rgb(123, 168, 194, maxColorValue = 255), rgb(110,   9, 123, maxColorValue = 255))
   ##################################################################################
   #Internal function: PEAKFINDER (finds all local minima/maxima in a range of data )
   ##################################################################################
     findpeaks <- function(vec, bw = 1, x.coo = c(1:length(vec))){  #where bw = is box width, setting the sensitivity of the search
                    ###set all vectors to null
                   	pos.x.max <- NULL ;	pos.y.max <- NULL ;	pos.x.min <- NULL ;	pos.y.min <- NULL
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
                      #are ALL trailing data smaller than i?
                  		is.max   <- sum(subset.inf > vec[i]) == 0
                  		#are ALL leading data smaller than i?
                   		is.nomin <- sum(subset.sup > vec[i]) == 0
                      #are ALL trailing data larger than i?
                   		no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
                      #are ALL leading data larger than i?
                  		no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
                    ##############################################################
                      #a maximum is found if  all data before and after i are smaller than i
                  		if(is.max & is.nomin){
                  			pos.x.max <- c(pos.x.max, x.coo[i])
                   			pos.y.max <- c(pos.y.max, vec[i])
                         }
                      #a maximum is found if  all data before and after i are larger than i
                  		if(no.max & no.nomin){
                  			pos.x.min <- c(pos.x.min, x.coo[i])
                  			pos.y.min <- c(pos.y.min, vec[i])
                        }
                     }#end of for loop
                      ###Output
                       return(list("max.X" = pos.x.max, "max.Y" = pos.y.max, "min.X" = pos.x.min, "min.Y" = pos.y.min))}
   #End if peakfinder function########################################################                                        
   ##################################################################################
   require(SuppDists)     #required for the 'moments()' function
   ##################################################################################
   ## ACTUAL FUNCTION START #########################################################
   ##################################################################################   
   #pre-perp: remove NA from data, inspect bandwidth   
     drp <- na.omit(drp)
     fail <- 0 #create failure tracker
     pops <- 0 #create population tracker
    #be sure bandwidth for kernel density is at least 50
     temp <- bw.nrd0(drp)
     if(temp < 50){
        temp <- 50
     } 
    #kernel density estimation
     krn <- density(drp, bw = temp)
     krn <- rbind(krn$y, krn$x)
    #we use findpeaks with a high bandwith to single out the populations
     piek <- findpeaks(krn[1, ], bw=20)$max.X
      dal <- findpeaks(krn[1, ], bw=20)$min.X
     piek <- rbind(piek, krn[1, piek])    #add peak heights
      dal <- rbind(dal,  krn[1,  dal])    
     piek <- rbind(piek, krn[2, piek[1,]])#add peak x-locations
      dal <- rbind(dal,  krn[2,  dal[1,]])     
    #we also remove all peaks that are smaller than 1% of the max peak height. 
    #Outlying fluorescence values will otherwise cause insignificantly small peaks that screw with algorithm robustness
    if(any(piek[2,] < max(piek[2,])/100)){
      piek <- piek[, -which(piek[2,] < max(piek[2,])/100)] 
      #make sure shit didn't vectorize
      if(is.null(dim(piek))){
        piek <- as.matrix(piek)
       }   
    }
   ##################################################################################
   ##First major control point: we only proceed if there are no NA
   if(!any(is.na(piek))){ ##Control LVL 1
   ##################################################################################  
   #First, we check for three or more populations
     if(length(piek[1,]) >= 3){
      #post warning
       message("WARNING: three population detected") 
        message("---------")
        #set the number of populations
        pops <- length(piek[1,])      
        #split total bandwidth in three parts, using lower and upper for neg and pos populations 
        temp <- c(max(piek[3, ]) - diff(tail(sort(piek[3, ]), n = 2))/2,      #limit between pos and middle
                  min(piek[3, ]) + diff(head(sort(piek[3, ]), n = 2))/2)      #limit between neg and middle                                                                          
        #The presence of the 'middle' population prevents us from using the iterative approach to set the standard deviations   
        #instead, we directly calculate the final values from the limits we've set
        ndrp <- drp[which(drp < temp[2])]
        pdrp <- drp[which(drp > temp[1])]
        #calculate kurtosis of droplets
        k1 <- abs(moments(ndrp)[4])
        k2 <- abs(moments(pdrp)[4])
        #calculate s multiplier
        s1 <- 3.8+0.35*log(k1)+0.045*log(k1)^2   +0.75
        s2 <- 3.8+0.35*log(k2)+0.045*log(k2)^2   +0.75
        #construct piik matrix      
        piik <- c(median(ndrp),median(pdrp))
       #Expand piik to take all values
       piik <- rbind (piik, matrix(c(NA,NA,NA,NA),nrow=2,ncol=2))   
       #fill in the values
       piik[2, 1] <-  piik[1, 1] -  mad(ndrp) * s1 #leftmost extreme (start of neg band)
       piik[2, 2] <-  piik[1, 2] + mad(pdrp) * s2 #rigthmost extreme (end of pos band)
       piik[3, 1] <-  piik[1, 1] + mad(ndrp) * s1 #first inner extreme (end of the neg band)
       piik[3, 2] <-  piik[1, 2] - mad(pdrp) *s2 #second inner extreme (start of the pos band)
        #Extra check: the start of the pos band should not fall in the 'middle' peak, 
        #             if so we replace it with the minimum between both peaks
        if(piik[3, 2] <  temp[1]){    
           piik[3, 2] <- max(dal[3 ,which(dal[3, ] < piik[1, 2])])  #maximum of all valleys smaller than pos pesk
           }                                     
       #We now calculate the threshold for counting positives      
        tresh <- piik[1, 1] + (s1 + 1/2 * s1) * mad(ndrp)
       #should the threshold fail, we provide a backup
        tresh.bup1 <- temp[2] 
        tresh.bup2 <- max(dal[3 ,which(dal[3, ] < piik[1, 2])])  #maximum of all valleys smaller than pos pesk   
        #finally, for compatibility's sake, we also reduce the original peak matrix to only contain the positive an negative populations
        piek <- piek[, c(which.min(piek[3, ]), which.max(piek[3, ]))]             
   }else{ # We have changed 'piek' to only contain two populations, we must skip the other possibilities to prevent re-analysis
   #END of THREE population case
   ##################################################################################   
   #If not three, check for exactly two populations (Standard algorithm)
     if(length(piek[1,]) == 2){
     #set the number of populations
     pops <- 2   
     #get started on the main population matrix
     piik <- piek[3, ]
     #decide on a temporarly divider for the peaks 
     kkn <- max(piek[3, ]) - diff(piek[3, ])/2
     #We find the peak base by first estimaging the standard deviation as 1/2 of the width of the peak at 0.6065 percent of the max height
     #first we find the y-location (60.65% height) for both peaks
     baes <- (dnorm(1)/dnorm(0)) * piek[2, ]
      temp <- head(krn[2, which(krn[2, ] < kkn)][order((krn[1, which(krn[2, ] < kkn)] - baes[1])^2)] , n = 10) #select 10 possible crossing candidates in the first half of the data (necessary as the data are discrete)
      temp <- c(temp[which(temp < piik[1])][order(temp[which(temp < piik[1])]-piik[1], decreasing = T)[1]],
                temp[which(temp > piik[1])][order(temp[which(temp > piik[1])]-piik[1])[1]])#take the ones nearest to the peak (one larger and one smaller)
     baes <- rbind(baes,temp)
      temp <- head(krn[2, which(krn[2, ] > kkn)][order((krn[1, which(krn[2, ] > kkn)] - baes[1, 2])^2)], n = 10) #select 10 possible crossing candidates in the second half of the data
      temp <- c(temp[which(temp < piik[2])][order(temp[which(temp < piik[2])] - piik[2], decreasing = T)[1]],
                temp[which(temp > piik[2])][order(temp[which(temp > piik[2])] - piik[2])[1]])#take the ones nearest to the peak (one larger and one smaller)
     baes <- rbind(baes,temp)
     #Expand piik to take all values
     piik <- rbind (piik, matrix(c(NA,NA,NA,NA),nrow=2,ncol=2))
     #The higher the kurtosis of the clouds the more sigmas we need to cover 99% of all the droplets, we start with standard values (4) and update below
     s1 <- 4
     s2 <- 4
     #let's populate the population matrix with the population boundaries    
     piik[2, 1] <-  piik[1, 1] - diff(baes[2, ])/2 * s1 #leftmost extreme (start of neg band)
     piik[2, 2] <-  piik[1, 2] + diff(baes[3, ])/2 * s2 #rigthmost extreme (end of pos band)
     piik[3, 1] <-  piik[1, 1] + diff(baes[2, ])/2 * s1 #first inner extreme (end of the neg band)
     piik[3, 2] <-  piik[1, 2] - diff(baes[3, ])/2 * s2 #second inner extreme (start of the pos band)
     #allright, using this rough delimination of the clouds we can calculate more accurate values
     #by iterating the update process three times (i.e. use current estimate to calculate pos and neg,
     # then recaculte sd and re-delimite bands). Three iterations is generaly enough to reach
     # stability in the first decimal place.
     for ( i  in 1 : 3 ) {
      #split in positive and negatives
       pdrp <- drp[which(drp >= piik[3, 2] & drp <= piik[2, 2])]
       ndrp <- drp[which(drp >= piik[2, 1] & drp <= piik[3, 1])]
       #calculate kurtosis of droplets
       k1 <- abs(moments(ndrp)[4])
       k2 <- abs(moments(pdrp)[4])
       #update s multiplier
       s1 <- 3.8+0.35*log(k1)+0.045*log(k1)^2   +0.75
       s2 <- 3.8+0.35*log(k2)+0.045*log(k2)^2   +0.75
       #update
       piik[2, 1] <-  piik[1, 1] - mad(ndrp) * s1 #leftmost extreme (start of neg band)
       piik[2, 2] <-  piik[1, 2] + mad(pdrp) * s2 #rigthmost extreme (end of pos band)
       piik[3, 1] <-  piik[1, 1] + mad(ndrp) * s1 #first inner extreme (end of the neg band)
       piik[3, 2] <-  piik[1, 2] - mad(pdrp) * s2 #second inner extreme (start of the pos band)
      #advance the cycle
      i <- i + 1
      }#END of the iterative update   
   #We now calculate the threshold for counting positives      
     tresh <- piik[1, 1] + (s1 + 1/2 * s1) * mad(ndrp)
     #should the threshold fail, we provide a backups
     tresh.bup1 <- min(piek[3,]) + 0.4 * abs(diff(piek[3,]))  #minimum peak + 40 % empty bandwidth
     tresh.bup2 <- min(dal[3 ,which(dal[3, ] > piik[1, 1])])  #minimum of all valleys larger than neg peak   
   }#END of the two population procedure
   ##################################################################################   
   #Last option: there is only one population (single peak) 
   #this case requres a bit more inspection to make sure we didn't overlook the some (few) positives or negatives
   if(length(piek[1,]) == 1){
     #set the number of populations
     pops <- 1     
     #Since there still may be few negatives / positives that got smoothed out. 
     #We should be able to detect them as 'outliers' from the single population
     #therefore we calculate z-scores for all droplets
     zs <- (drp - median(drp))/mad(drp)
     k1 <- abs(moments(drp)[4])
     s1 <- 4.55 + 0.35 * log(k1) + 0.045 * log(k1)^2    
     #check for outliers and check if they sare higher or lower to decide if the pop is pos or neg
     #Also make sure there is more than one outlier, otherwise none of the below makes sense
     if(any(abs(zs) > 1.5 * s1) & length(which(abs(zs) > 1.5 * s1)) >= 2){
      #we have outliers! 
      #set the number of populations
        pops <- pops + 1          
      #we now decide if they are positive or negative and produce threshold and piik matrix
      if(sum(drp[which(abs(zs) >= s1)]) > 0){ #outliers have higher fluorescence and are thus positive
         pdrp <- drp[which(abs(zs) >= 1.5 * s1)]
         ndrp <- drp[which(abs(zs) <  1.5 * s1)]      
      }else{ #outliers have lower fluorescence and are thus negative
         pdrp <- drp[which(abs(zs) <  1.5 * s1)]
         ndrp <- drp[which(abs(zs) >= 1.5 * s1)]
      }
        k1 <- abs(moments(ndrp)[4])
        k2 <- abs(moments(pdrp)[4])
        #calculate s multiplier
        s1 <- 3.8+0.35*log(k1)+0.045*log(k1)^2   +0.75
        s2 <- 3.8+0.35*log(k2)+0.045*log(k2)^2   +0.75
        #construct piik matrix      
         piik <- c(median(ndrp),median(pdrp))
         #Expand piik to take all values
         piik <- rbind (piik, matrix(c(NA,NA,NA,NA),nrow=2,ncol=2))   
         #fill in the values
         piik[2, 1] <-  piik[1, 1] -  mad(ndrp) * s1 #leftmost extreme (start of neg band)
         piik[2, 2] <-  piik[1, 2] + mad(pdrp) * s2 #rigthmost extreme (end of pos band)
         piik[3, 1] <-  piik[1, 1] + mad(ndrp) * s1 #first inner extreme (end of the neg band)
         piik[3, 2] <-  piik[1, 2] - mad(pdrp) *s2 #second inner extreme (start of the pos band)
        #we also need to add a column do 'piek' for compatibility's sake, position depends on whether outlers are pos or neg
        if(length(ndrp) > length(pdrp)){ 
          piek <- cbind(piek, c(krn[1, which.min((krn[2, ]-piik[1, 2])^2)], piek[2, ]/100, piik[1, 2]))
         }else{
          piek <- cbind(c(krn[1, which.min((krn[2, ]-piik[1, 1])^2)], piek[2, ]/100, piik[1, 1]), piek)
         } #END piek update
         #We now calculate the threshold for counting positives      
         tresh <- piik[1, 1] + (s1 + 1/2 * s1) * mad(ndrp)
     }else{ #no outliers: we give up
     fail <- 1 #notify failure
     } #end of no outlier case
   } #END of ONE population case
   } #END of the THREE population 'if - else'    
   #from here on it's the same for all cases, we start with:
   ##################################################################################   
   ##Second major control point: checks if 'tresh' exists (if not: only one population, no outliers)
   ################################################################################## 
   if(exists("tresh")){##Control LVL 2 
   #if we have a threshold we can continue, we start with a few more controls:     
   ##################################################################################   
   ##third to fifth control points: checks for stability reasons
   ################################################################################## 
   #in case of failure (tresh is NA) we fal back on the backup threshold 
   if(is.na(tresh)){##Control LVL 3
        if(!isTRUE(silent)){message("WARNING!: threshold placement failed")
                           message("Check data, arbitrary threshold was set")
                           message("---------")} 
        tresh <-  tresh.bup1
        ## if tresh is NA the piik is as well, we redo all the calculations:
        ndrp <- drp[which(drp < tresh)]
        pdrp <- drp[which(drp > tresh)]
        #construct piik matrix      
        piik <- c(median(ndrp),median(pdrp))
        #Expand piik to take all values
        piik <- rbind (piik, matrix(c(min(ndrp),max(pdrp),max(ndrp),min(pdrp)),nrow=2,ncol=2, byrow = T)) 
    }##Control LVL 3 END 
    ############################## 
    #excessive rain can make the standard deviations go haywire 
    #we build in this simple check to make sure the treshold is not placed into the positive band
    if(tresh > piik[3, 2]){ ##Control LVL 4
        if(!isTRUE(silent)){message("WARNING!: threshold inside positve cloud")
                           message("Forced lower threshold, not based on population parameters")
                           message("---------")} 
        tresh <-  tresh.bup2
    }##Control LVL 4 END 
    ##############################     
    #another problem related to standard deviations is overlapping bounderies
    #check : upper neg boundary is situated above lover pos boundary? 
    if(piik[3, 1] > piik[3, 2]){##Control LVL 5
        if(!isTRUE(silent)){message("WARNING!: population boundary overlap")
                           message("boundaries rewritten, not based on population parameters")
                           message("---------")} 
      #we can use the 'valleys' returned by findpeaks to draw new population bounderies
      #we take the valley closest to the corresponding peak as the new boundary
      piik[3, 1] <- min(dal[3 ,which(dal[3, ] > piik[1, 1])])  #minimum of all valleys larger than neg peak
      piik[3, 2] <- max(dal[3 ,which(dal[3, ] < piik[1, 2])])  #maximum of all valleys smaller than pos pesk      
    }##Control LVL 5 END    
   ##################################################################################    
   #we can now start the trivial calculations:
   ##################################################################################                   
    #count the raindrops (i.e. the population of the 'clearance band' between the positive and negative band)
     #again to avoid including negatives in the rain we use the treshold in stead of the peak base
     #rain <- length(which(drp > piik[3, 1] & drp < piik[3, 2]))
      rain <- length(which(drp > tresh & drp < piik[3, 2]))
    #split in positive and negatives so we can apply colours on the plot
     pdrp <- drp[which(drp>tresh)]
     ndrp <- drp[which(drp<=tresh)]
     tdrp <- length(drp)
    #also define 'cerntainly pos', 'rain' and 'drift'
     cdrp <- drp[which(drp > piik[3, 2] & drp < piik[2, 2])]  #c-pos is the population of the pos band only
     rdrp <- drp[which(drp > tresh & drp < piik[3, 2])]  #rain is the population of the clearance band (counted positive)
     ddrp <- drp[which(drp > piik[2, 2])]                     #drift is population above the pos band (high-Fluorescence outliers)
 ###Calculation of Performance parameters
    #1#  Resolution
     Resol <- 2 * (piik[1, 2] -piik[1, 1]) / ((piik[3, 1] - piik[2, 1]) + (piik[2, 2] - piik[3, 2]))
    #2# percentage rain
     pRain <- rain / tdrp
    #3# Compartmentalization
     Comp <- tdrp * dVol * 0.001 / sVol
 ###Calculation of Lambda and confidence bounds
     lambda <- -log(length(ndrp)/tdrp)
      lamlo <- lambda - 1.96 * sqrt((tdrp - length(ndrp))/(tdrp * length(ndrp)))   #lower 95% confidence bound for lambda 1 (cpd)
      lamhi <- lambda + 1.96 * sqrt((tdrp - length(ndrp))/(tdrp * length(ndrp)))   #upper 95% confidence bound for lambda 1 (cpd)
      lambs <- c(lambda, lamhi, lamlo)
  #total number of targets present in droplets analysed
    targets <- tdrp * lambs
  #total number of targets in sample volume
    targall <- targets * (sVol/(tdrp * dVol * 0.001))
 ###gather and name all output
      names(targets) <- c("targets", "upper", "lower")
      names(targall) <- c("targ.in.sVol", "upper", "lower")
      names(lambs) <- c("lambda", "upper", "lower")
     drops <- c(length(pdrp), length(ndrp), rain)
      names(drops) <- c("positive", "negative", "rain")
     perform <- c(Resol, pRain, Comp)
     names(perform) <- c("resolution", "%rain", "%comparted")
     names(pops) <- "populations"
 ##################################################################################                   
 ### Make plots (optional)
    if (isTRUE(plots)) {
     #plotting window:
     x11(14,7); par(mfrow=c(1, 2))
   #1#kernel density
      plot(krn[1, ] ~ krn[2, ], xlab = "FU", ylab = "density", bty = "L", type = "l", lwd = 1.5)
     #add limits to plot
      abline(h = 0)
      polygon(x = c(piik[2:3, 1], piik[3:2, 1]), y = rep(c(0, piek[2, 1]), each = 2), col = rgb(68,  140,  59, maxColorValue = 255, alpha = 0.3 * 255),  border = NA)
      abline(v = piik[1, ],  col = two.color[1], lty = 2)
      polygon(x = c(piik[2:3, 2], piik[3:2, 2]), y = rep(c(0, piek[2, 2]), each = 2), col = rgb(133,  84, 168, maxColorValue = 255, alpha = 0.3 * 255),  border = NA)
      abline(v = piik[2:3,], col = two.color[8], lty = 2)
     #add rain kernel
     if(!length(rdrp) == 0){
      temp <-             density(c(rep(0,times=length(ndrp)),rdrp,rep(0,times=length(cdrp)+length(ddrp))))$y    #add a bunch of zeros to get y-scale to fit other plot
      temp <- rbind(temp, density(c(rep(0,times=length(ndrp)),rdrp,rep(0,times=length(cdrp)+length(ddrp))))$x)
      temp <- temp[, which(temp[2,] >= min(rdrp) & temp[2,] <= max(rdrp))]#trim off the excess caused by the zeros
      polygon(x = c(min(rdrp), temp[2,], max(rdrp)), y = 1.5 * c(0, temp[1,], 0), col = rgb(245, 255,  46, maxColorValue = 255, alpha = 0.4 * 255),  border = NA)
     }
     #retrace kernel lines on top
      lines(krn[1, ] ~ krn[2, ], lwd = 1.5)
     #Add Letter
      mtext("A",side=3,cex=1.5,adj=1)
   #2#droplet view
     temp <- sample(c(1:length(drp)), size = length(cdrp), replace = F)
      plot(cdrp ~ temp, col = spectral[6], bty = "L", pch = 20 , xlab = "observation", ylab = "FU", ylim = c(min(drp), max(drp)))
     temp <- sample(c(1:length(drp)), size = length(ndrp), replace = F)
      points(ndrp ~ temp, pch = 20, col = spectral[1])
     if(!length(rdrp) == 0){
       temp <- sample(c(1:length(drp)), size = length(rdrp), replace = F)
        points(rdrp ~ temp, pch = 20, col = spectral[4])
     }
     temp <- sample(c(1:length(drp)), size = length(ddrp), replace = F)
      points(ddrp ~ temp, pch = 20, col = spectral[8])
     abline(h=tresh, col = "gray40", lty = 3)
     #Add Letter
      mtext("B",side=3,cex=1.5,adj=1)
  #2#barplots of performance
    colrs <- c(rep(spectral[6],times=7),spectral[1])
    if(Resol > 2.5){colrs[1:2] <- two.color[c(1,1)]}
    if(pRain < 0.025){colrs[3:4] <- two.color[c(1,1)]}
    if(0.3 < Comp & Comp < 0.5){colrs[5:6] <- spectral[c(3,3)]}
    if(Comp > 0.5){colrs[5:6] <- two.color[c(1,1)]}
    x11()
    barplot(matrix(c(Resol/10, rep(0, times = 9),
                     pRain*10, rep(0, times = 9),
                     Comp,     rep(0, times = 9),
                     drops[-3]/tdrp), nrow=8, ncol=4, byrow=F),
            col = colrs, bty = "L")
    lines(x = c(0,   2.4), y = c(0.25, 0.25), lty = 2)
    lines(x = c(2.5, 3.7), y = c( 0.3, 0.3),  lty = 2)
    lines(x = c(2.5, 3.7), y = c( 0.5, 0.5),  lty = 2)
    mtext(c("Resolution", "%Rain", "%compartmented","%Count"),at=c(0.7, 1.85, 3.1, 4.25),side=1)
    }
   ##################################################################################                       
   #Standard INFO messages
    if(!isTRUE(silent)){
      if(Resol < 2.5){  message("Resolution: FAILED")}else{message("Resolution: PASSED")}
      if(pRain < 0.025){message("Rain: PASSED")      }else{message("Rain: FAILED")}
      if(Comp < 0.3){message("Compartmentalization: FAILED")}else{
      if(Comp < 0.5){message("Comparmentalization: MEDIOCRE")}else{
                     message("Comparmentalization: PASSED")}}
      if(lambs[1] < 0.01){message("Confidence bounds: NOT optimal")}else{
      if(lambs[1] < 5){message("Confidence bounds: OPTIMAL")}else{
                     message("Comparmentalization: NOT optimal")}}
      }#end of messages
   ##################################################################################                       
   #ERROR messages
   }else{
    message("WARNING: only one population detected, returning NA")
   }##Control LVL 2 END
   }else{ ##Control LVL 1
    message("WARNING: something is not right, please check data")
    fail <- 1 #notify failure
   }##Control LVL 1 END
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
   } #End of failure (NA) procedure
   ##################################################################################     
   #actual output ##################################
    if(isTRUE(vec)){ #output if
        return(c(targall, lambs, perform, drops, pops))
    }else{
        return(list("targets.in.sVol" = targall, "lambda" = lambs, "performance" = perform, "droplets" = drops, "populations" = pops))
    }#output if END                           
  }##END of 'cloudy'
    
    
