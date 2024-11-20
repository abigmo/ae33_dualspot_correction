## Author: Alessandro Bigi (email: alessandro.bigi@unimore.it)
## Institution: University of Modena and Reggio Emilia, Modena, Italy
## Project: BlackAir
## Date: Oct 2022

## Suite of R functions to apply the Dual Spot® correction proposed in Drinovec et al. (doi:10.5194/amt-8-1965-2015)
## to raw attenutation of MicroAethalometer MA200 (firmware v1.4) at a desired time step

## required libraries
library(dplyr)
library(magrittr)
library(tidyr)
library(Rfast)
library(lubridate) 

### start getting the data from here ### 
input <- ma201

## list of variables which need to be averaged
## names(input)[(12, 13:15, 18:21, 24:31)] + date
ave.var <- c("date",
             "GPS.lat..ddmm.mmmmm.",
             "GPS.long..dddmm.mmmmm.",
             "GPS.speed..km.h.",
             "GPS.sat.count",
             "Battery.remaining....",
             "Accel.X",
             "Accel.Y",
             "Accel.Z",
             "Flow.total..mL.min.",
             "Flow1..mL.min.",
             "Flow2..mL.min.",
             "Sample.temp..C.",
             "Sample.RH....",
             "Sample.dewpoint..C.",
             "Internal.pressure..Pa.",
             "Internal.temp..C.")

## list of variables to be taken at the end of the measurement
## names(input)[c(1:32)[-c(12, 13:15, 18:21, 24:31)]]
## grep("Sen|ATN|Ref", names(input), value = TRUE)
end.var <- c("date",
             "Serial.number",
             "Datum.ID",
             "Session.ID",
             "Data.format.version",
             "Firmware.version",
             "App.version",
             "Date...time.local",
             "Timezone.offset..mins.",
             "Date.local..yyyy.MM.dd.",
             "Time.local..hh.mm.ss.",
             "Timebase..s.",
             "Status",
             "Tape.position",
             "Flow.setpoint..mL.min.",
             "Optical.config",
             "Readable.status",
             "UV.Sen1",
             "UV.Sen2",
             "UV.Ref",
             "UV.ATN1",
             "UV.ATN2",
             "Blue.Sen1",
             "Blue.Sen2",
             "Blue.Ref",
             "Blue.ATN1",
             "Blue.ATN2",
             "Green.Sen1",
             "Green.Sen2",
             "Green.Ref",
             "Green.ATN1",
             "Green.ATN2",
             "Red.Sen1",
             "Red.Sen2",
             "Red.Ref",
             "Red.ATN1",
             "Red.ATN2",
             "IR.Sen1",
             "IR.Sen2",
             "IR.Ref",
             "IR.ATN1",
             "IR.ATN2")


## date must be in POSIXct in UTC

a <- input %>%
    dplyr::select(all_of(ave.var)) %>%
    mutate("time.ave" = lubridate::floor_date(date, "min")) %>% ## time stamps are reported at the end of the measurement. Get the start time of each measurement to eventually group measurement together (if you want to compute BC at different timestep of the experimental sampling time)
## old code, not using lubridate: it works fine, it's just slightly more tricky to handle  ## mutate("time.ave" = as.POSIXct(floor ( as.numeric( date - 60 ) / 60 / 5) * 5 * 60, origin = "1970-01-01") + 60 * 5) %>% ## in this case we find to which '5 minute slot' each observation belongs
    group_by(time.ave) %>% ## group the data according to this variable, i.e. for dates rounded to 5 minutes
    summarise_if(.predicate = function(x) is.numeric(x),
                 .funs = ~ quantile(x = ., probs = 0.5, na.rm = TRUE)) %>% 
    mutate(date = as.POSIXct(time.ave, tz = "UTC")) %>%
    dplyr::select(all_of(ave.var)) %>%
    as.data.frame()

b <- input %>%
    dplyr::select(all_of(end.var)) %>%
    mutate("time.ave" = lubridate::floor_date(date, "min")) %>% ## time stamps are reported at the end of the measurement. Get the start time of each measurement to eventually group measurement together (if you want to compute BC at different timestep of the experimental sampling time)
    ## mutate("time.ave" = as.POSIXct(floor ( as.numeric( date - 60 ) / 60 / 5) * 5 * 60, origin = "1970-01-01") + 60 * 5) %>% ## in this case we find to which '5 minute slot' each observation belongs
    group_by(time.ave) %>% ## group the data according to this variable, i.e. for dates rounded to 5 minutes
    summarise_all(, .funs = ~ last(.)) %>%
    mutate(date = as.POSIXct(time.ave, tz = "UTC")) %>%
    dplyr::select(all_of(end.var)) %>%
    as.data.frame()

## ATN by Aethlabs is computed as
## -100 * log(input$w.Sen / input$w.Ref) + 100*log(input$w.Sen[1]/input$w.Ref[1]) ## ATN from eq. 1 in Drinovec et al, reported as variation compared to the sensor count at the start of the measurement (which ideally it should be a new spot)

## merge the two datasets
input <- merge(a,b, by = "date")
rm(a, b, end.var, ave.var)

## update the timebase if needed
input$Timebase..s. <- 300

## compute the attenuation coefficient [Mm-1]
## using Drinovec et al. eq. 15, eq. 16
## from the ATN by aethlabs at raw timebase. 10^6 converts to ng/g. 10^6 converts from m-1 to Mm-1
## this is the not-corrected attenuation coeff which, multiplied by sigma, leads to BC_not_corrected (so it corresponds as BC_NC in eq. 7 of Drinovec et al.)
input <- input %>%
    mutate_at( dplyr::vars(ends_with("ATN1")), .funs =  list( `batn` = function(x, flow, time) {
        c( diff(x), diff(x)[1] ) / 100 * 0.71*1E-5 / (flow / 10^6) / (time / 60) * 10^6 }), ## [time] is [s] and is / by 60 to convert it to min. flow is [ml/min] it's "/10^6" to take it to [m3/min], last 10^6 is to convert to Mm-1
        flow = quote(Flow1..mL.min.), time = quote(Timebase..s.)) %>% ## computing explicitly the attenutation coefficient based on the ATN1 and label it "UV.ATN1_batn"
    mutate_at( dplyr::vars(matches("batn")), ~ ifelse(. > 0, ., NA))  ## then removing the negative attenuation coefficient at each spot change

## estimate of FVRF (flow velocity ratio factor), using also fig S6
## ATN lower and upper limit for FVRF estimate
## this should be roughly independent of wavelengths
## Chakrabarty et al, 2023 AMT suggests 15 and 30
## AE33 by default should have 10 and 30
atn.f1 <- 15
atn.f2 <- 30

bla <- input %>%
    split(., f = input$Session.ID) %>%
    lapply(., function(x) dualspot.compensed.bc(x, atn.f1 = 15, atn.f2 = 30)) ## this is the core function for dual spot correction

output <- do.call(rbind.data.frame, bla)


dualspot.ae33.bc <- function(data, atn.f1, atn.f2) {

    wv <- c("UV" = 375, "Blue" = 470, "Green"= 528, "Red" = 625, "IR" = 880) ## communicated by Aethlabs
    c.ref <- 1.3 ## communicated by Aethlabs
    master.sigma <- c("UV" = 24.069, "Blue" = 19.070, "Green" = 17.028, "Red" = 14.091, "IR" = 10.120) ## MAC used by Aethlabs

    ## Estimate of not-weighted k (eq. 10 Drinovec AMT 2015)
    k.unw.estimate <- function(k = seq(1E-1, 1E-4, length.out = 1000), f1, f2, atn1, atn2, fvrf){

        ## the correction is BC_not_compensated / (1 - k * ATN1)
        ## therefore I need to check that k * ATN1 > 1
        ## I only need to check for ATN2.
        ## if it is met for ATN2, then it's met also for ATN1, since ATN1 > ATN2
        k <- k[!(k * atn1 > 1)]

        ## compute equation 10 Drinovec AMT 2015
        ## if FVRF is NA, then there is no need of the correction (no loading, ATN is very low)
        if (is.na(fvrf) == TRUE) {
            k.min  <- 0 
            return(k.min)
        }
        y <- abs(f2 / f1 * fvrf - log(1 - k * atn2) / log(1 - k * atn1))

        ## get the minimum y
        y.min <- Rfast::nth(y, 1, descending = F)

        ## get the position of the minimum y and extract the k
        k.min <- k[Rfast::min_max(y, index = TRUE)[1]]

        return(k.min)
    }

    ## function to weigh the K according to Drinovec 
    k.weight <- function(k.old, atn.ta, atn.f2, atn1, k){

        ## applying Drinovec Supplement lines 170, 174
        out <- ((atn.ta - atn1) * k.old + (atn1 - atn.f2) * k) / (atn.ta - atn.f2) 
        out[ which( atn1 < atn.f2 )] <- k.old
        out   
    }
    
    for (w in names(wv)){

        a1 <- paste0(w,".ATN1")
        a2 <- paste0(w,".ATN2")

        dd <- split(data, f = data$Tape.position)

        ## this is the intercept
        f1 <- lapply(dd, function(x) {

            y <- x[which(x[ , names(x) %in% a1] > atn.f1 & x[ , names(x) %in% a2] < atn.f2), ] ## select the ATN1 within the limits atn.f1 and atn.f2

            if (sum(complete.cases(y)) == 0) {
                message(paste0("Warning: for ", w, " there is no ATN between ", atn.f1, " and ", atn.f2, " for spot ", unique(x$Tape.position), " and Session ", unique(x$Session.ID), ". Probably ATN is very low. K will be set to zero"))
                f1 <- NA
            } 
            
            if (sum(complete.cases(y)) != 0) f1 <- lm( c( y[ , names(y) %in% a2] / y[ ,names(y) %in% a1] ) ~ y[ , names(y) %in% a1] )$coefficients[1] ## slope of the linear regression between ATN2 / ATN1 ~ ATN1 (within the range atn.f1 and atn.f2)
            return(f1)
        }) %>%
            unlist()

        ## this is the ratio of the flows
        f2 <- lapply(dd, function(x) {
            median(x[ ,"Flow1..mL.min."]) / median(x[ ,"Flow2..mL.min."])
        }) %>%
            unlist()
        
        ## fvrf is the product of f1 and f2 (eq. 12 Drinovec)
        fvrf <- f1 * f2
        
        assign(paste0("fvrf.",w), as.numeric(fvrf))
    }

    rm(f1, f2, a1, a2, w, dd, fvrf)

    ## container list with the results
    res <- NULL

    ## initialise the old k compensation factor
    k.old <- data.frame("UV" = NA, "Blue" = NA, "Green"= NA, "Red" = NA, "IR" = NA)

    ## loop through the spots
    for (tp in unique(data$Tape.position)) {

        dd <- data %>%
            dplyr::filter(Tape.position == tp)

        ## loop through the wavelengths
        for (w in names(wv)){

            k.old.w <- k.old[1, match(w, names(k.old))] ## get the k.old for the current wavelength "w"

            sigma <- master.sigma[match(w, names(wv))] ## get the sigma corresponding to the wavelength "x"

            fvrf  <- paste0("fvrf.", w) %>%
                get(.) %>% ## get the list of fvrf for all spots in this Session ID, for the wavelength "x"
                extract2( match(tp, unique(data$Tape.position)) )  ## get the fvrf of the spot number "tp"
            
            l1 <- paste0(w, ".ATN1") ## attenuation spot 1
            l2 <- paste0(w, ".ATN2") ## attenuation spot 2
            
            ## I need to use mapply in the computation of k, because the function is not vectorized
            ## I need to dynamically assign new names in mutate (see: https://stackoverflow.com/a/26003971/5375295)
            dd <- dd %>%
                dplyr::filter(Tape.position == tp) %>%
                mutate("{w}.K.ae33" := mapply(k.unw.estimate, f1 = Flow1..mL.min., f2 = Flow2..mL.min., atn1 = get(l1), atn2 = get(l2), fvrf = fvrf) ) %>% ## compute the unweighted K (Drinovec et al. eq.10)
                mutate("{w}.K.ae33" := k.weight(k.old = ifelse(!is.na(k.old.w), k.old.w, mean( tail( get(paste0(w, ".K.ae33")), 10))), ## compute the weighted K (Drinovec et al. eq.13)
                                                atn.ta = max( get(l1), na.rm = TRUE),
                                                atn.f2 = atn.f2,
                                                atn1 = get(l1),
                                                k = get(paste0(w, ".K.ae33")) )) %>%
                mutate( "{w}.BCc.ae33":= get(paste0(w, ".ATN1_batn")) / c.ref / sigma / (1 - get(paste0(w, ".K.ae33")) * get(l1)) * 1000 ) ## BC corrected in µg/m³

            ## update the value in the k.old for the following spot
            ## using the mean of the last 10 values of K (e.g. 50 minutes)  
            k.old[1, match(w, names(k.old))]  <- dd %>%
                dplyr::select( matches(paste0(w, ".K.ae33"))) %>%
                tail(., 10) %>%
                summarise_all( ~ mean(., na.rm = TRUE)) %>%
                as.numeric(.)
        }

        ## append result from the spot to a list with everything
        res[[ length(res) + 1]] <- dd
        
    }

    res <- do.call(rbind.data.frame, res)

    res %>%
        dplyr::mutate_at(dplyr::vars(matches("BC")), ~ ifelse(grepl(c("Tape advance|Start up|Optical saturation"), Readable.status) == TRUE, NA, .) )

    rm(k.old, l1, l2, tp, fvrf, sigma, w)
    rm(fvrf.IR, fvrf.Green, fvrf.Red, fvrf.Blue, fvrf.UV, atn.f1, atn.f2)

    invisible(res)
}
