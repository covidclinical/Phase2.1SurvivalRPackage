require(tidyverse)
require(data.table)
require(RcppRoll)
require(zoo)
#' Generates AKI KDIGO grades using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

aki_kdigo_grade <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_90d,2)
  if(diff > 0.3 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades using future serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in the next 7 days and next 48h
#' @noRd

aki_kdigo_grade_retro <- function(x) {
  # Instead of using the pure KDIGO definition, this looks at future Cr values, and 
  # determines whether the current value, in retrospect, represents an episode of AKI
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[7])
  baseline_48h = as.numeric(x[8])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_7d,2)
  if(diff > 0.3 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades but on serum creatinine values 7 days later using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

akd_grade_7d <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_90d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_7d = as.numeric(x[9])
  grade = 0
  ratio = round(cr_7d/baseline,2)
  diff = cr_7d - baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades but on serum creatinine values 90 days later using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

akd_grade_90d <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_90d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_90d = as.numeric(x[10])
  grade = 0
  ratio = round(cr_90d/baseline,2)
  diff = cr_90d - baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Gives the position of the lowest serum creatinine in a specified window
#' @param cr Vector containing serum creatinine values
#' @param day Vector containing time points of respective serum creatinine values
#' @param lag Specified whether to look at retrospective or future data. Default = TRUE
#' @param gap Specifies rolling window length to use. Default = 7 days
#' @noRd

pos_min <- function(cr,day,lag=TRUE,gap=7) {
  len = length(cr)
  day_pos = day
  for(i in 1:len) {
    if(lag) {
      j = max(1,i-gap)
      pos = j-1+which.min(cr[j:i])
      day_pos[i] = day[pos]
    } else {
      j = min(i+gap,len)
      pos = i-1+which.min(cr[i:j])
      day_pos[i] = day[pos]
    }
  }
  day_pos
}

#' Boolean function which returns whether a certain value is minimum/maximum
#' @param x Vector containing serial creatinine values
#' @param partial Specifies whether to count the start/end of the array as a min/max value
#' @param decreasing Specified whether to find min (TRUE) or max (FALSE)
#' @noRd

which.peaks <- function(x,partial=TRUE,decreasing=FALSE) {
  if(decreasing) {
    if(partial) {
      which(diff(c(FALSE,diff(x) > 0,TRUE)) > 0)
    } else {
      which(diff(diff(x)>0)>0) + 1
    }
  } else {
    if(partial) {
      which(diff(c(TRUE,diff(x) >= 0,FALSE)) < 0)
    } else {
      which(diff(diff(x)>=0)<0) + 1
    }
  }
}

#' Returns the time point at which the normalised serum Cr falls below a certain value
#' @param ratio Vector containing normalised serum Cr values
#' @param time_from_peak Day of interest, expressed as time from peak
#' @param target Threshold of interest, default = 1.25
#' @noRd
get_day <- function(ratio,time_from_peak,target=1.25) {
  index = purrr::detect_index(ratio,function(x) x <= 1.25)
  if(index > 0) {
    day = time_from_peak[index]
  } else {
    day = time_from_peak[length(time_from_peak)]
  }
  day
}

#' Returns a list of two data tables, aki_peak_summ (AKI events for AKI patients) and nonaki_peak_summ (Peak sCr for non-AKI),
#' each containing the following information
#' patient_num,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
#' 
#' patient_num - the patient-num corresponding to the AKI event
#' siteid - siteid (if multiple siteids were combined in same file)
#' days_since_admission - time at which peak Cr is achieved
#' day_min - time at which Cr begins to rise
#' value - peak sCr value in mg/dL
#' day_min - in the past 90 days, the day at which the lowest serum creatinine was achieved
#' day_min_retro - in the next 7 days, the day at which the lowest serum creatinine was achieved
#' min_cr_90d - the sCr value corresponding to day_min
#' min_cr_48h - lowest sCr value in the 2 days prior
#' min_cr_retro_7day - the sCr value corresponding to day_min_retro
#' min_cr_48h_retro - lowest sCr value in 2 days in future
#' min_cr_7d_final - lowest sCr value of min_cr_90d and min_cr_retro_7day
#' cr_7d - sCr value 7 days after peak
#' cr_90d - sCr value 90 days after peak
#' delta_cr - change in sCr value from baseline to peak
#' aki_kdigo - KDIGO grade corresponding to using min_cr_90d or min_cr_48h as baseline sCr
#' aki_kdigo_retro - KDIGO grade corresponding to using min_cr_retro_7day or min_cr_48h_retro as baseline sCr
#' aki_kdigo_final - Final KDIGO grade, taking the worse grade of aki_kdigo and aki_kdigo_retro. Use this value for the final AKI severity.
#' akd_7d - using the expanded acute kidney disease (AKD) definition proposed by Chawla et al, the AKD grade at 7 days post-peak
#' akd_90d - using the expanded AKD definition proposed by Chawla et al, the AKD grade at 90 days post-peak
#' 
#' @param demog Demographics table
#' @param obs Observations table

get_aki_events <- function(demog,obs){
  demographics <- demog
  observations <- obs
  message("Now proceeding to AKI detection code:")
  message("Extracting serum creatinine values...")
  # Extract serum Cr levels
  labs_cr_aki <- observations[observations$concept_code == '2160-0',] #LOINC code for Cr 2160-0
  # Remove unnecessary columns
  labs_cr_aki <- labs_cr_aki[,c("patient_num","siteid","days_since_admission","value")] ## CH
  # Filter for labs >= -90 days
  labs_cr_aki <- labs_cr_aki[labs_cr_aki$days_since_admission >= -60,]
  
  # Generate separate demographics table for patients who do not have any sCr values fulfilling 
  # the above (e.g. all the labs are before t= -90days or patient has no sCr value)
  message("Removing patients who do not have any serum creatinine values during admission...")
  pts_valid_cr <- labs_cr_aki$patient_num[labs_cr_aki$days_since_admission >= 0]
  # Keep only patients with more than one sCr value to allow determination of renal recovery
  pts_valid_cr <- pts_valid_cr[duplicated(pts_valid_cr)] 
  demographics <- demographics[demographics$patient_num %in% pts_valid_cr,]
  labs_cr_aki <- labs_cr_aki[labs_cr_aki$patient_num %in% pts_valid_cr,]
  
  # There are two possible scenarios which we have to consider when detecting each AKI event:
  # (1) AKI occurs after admission
  #   - easy to detect with the formal KDIGO definition as we only need to use older data points
  #   - we can use an approach similar to what is used in the MIMIC-III data-set: use a rolling
  #     time frame and detect the lowest Cr value to use as our baseline
  # (2) Patient presents with an AKI at the time of admission 
  #   - this makes it harder for us to determine what is the true baseline especially with limited
  #     longitudinal data
  #   - Hence one way to get around this is to generate a "retrospective" baseline (i.e. look into
  #     future Cr values and find the minimum in a rolling timeframe) and use this as a surrogate
  #     baseline
  #   - Such a workaround is the only feasible way of dealing with missing retrospective data though
  #     it is likely that we may miss quite a number of true AKIs using this method
  
  # Scenario (1): AKI occurs during admission
  # Find minimum Cr level in a rolling 90 day timeframe
  message("Generating minimum Cr level in the past 90 days")
  labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_num","days_since_admission"))
  labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_num),seq(min(days_since_admission)-90,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
  labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_num) %>% dplyr::mutate(min_cr_90d = RcppRoll::roll_min(value,91,fill=NA,na.rm=TRUE,align="right")) %>% dplyr::filter(!is.na(value)) %>% dplyr::ungroup()
  # Find minimum Cr level in a rolling 2 day timeframe (48h)
  message("Generating minimum Cr level in the past 48h")
  labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_num","days_since_admission"))
  labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_num),seq(min(days_since_admission)-2,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
  labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_num) %>% dplyr::mutate(min_cr_48h = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,align="right")) %>% dplyr::filter(!is.na(value)) %>% dplyr::ungroup()
  
  # Scenario (2): Patient presents with an AKI already on board
  # Find minimum Cr level in a rolling 7 day timeframe
  message("Generating minimum Cr level 7 days in the future")
  labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_num","days_since_admission"))
  labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_num),seq(min(days_since_admission),max(days_since_admission)+7))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
  labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_num) %>% dplyr::mutate(min_cr_retro_7day = RcppRoll::roll_min(value,8,fill=NA,na.rm=TRUE,align="left")) %>% dplyr::filter(!is.na(value)) %>% dplyr::ungroup()
  # Find minimum Cr level in a rolling 2 day timeframe (48h)
  message("Generating minimum Cr level 48h in the future")
  labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_num","days_since_admission"))
  labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_num),seq(min(days_since_admission),max(days_since_admission)+2))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
  labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_num) %>% dplyr::mutate(min_cr_48h_retro = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,align="left")) %>% dplyr::filter(!is.na(value)) %>% dplyr::ungroup()
  
  # Another outcome we are interested in is to look at acute kidney disease, AKD (in between AKI and CKD)
  # We will use the definitions proposed for AKD as described by Chawla et. al. 2017 (ref (1))
  # We are interested in renal recovery at the 7-day and 90-day timepoint
  # We will use a cutoff of recovery to 1.25x baseline Cr as recovery, as used in ref (2)
  # References:
  # 1. Chawla, L., Bellomo, R., Bihorac, A. et al. Acute kidney disease and renal recovery: consensus report
  #    of the Acute Disease Quality Initiative (ADQI) 16 Workgroup. Nat Rev Nephrol 13, 241-257 (2017). 
  #    https://doi.org/10.1038/nrneph.2017.2
  # 2. Pannu, N., James, M., Hemmelgarn, B. & Klarenbach, S. Association between AKI, Recovery of Renal 
  #    Function, and Long-Term Outcomes after Hospital Discharge. Clinical Journal of the American Society of
  #    Nephrology 8, 194-202 (2013). https://doi.org/10.2215/CJN.06480612
  
  # Generate sCr levels at +7d (cr_7d) and +90d (cr_90d) timepoints (for determining post-AKI recovery, AKD)
  labs_cr_aki <- data.table::setDT(labs_cr_aki)[,':='(cr_7d = tail(labs_cr_aki$value[labs_cr_aki$patient_num==patient_num][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_num==patient_num],days_since_admission,days_since_admission+7,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_num==patient_num][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_num==patient_num],days_since_admission,days_since_admission+90,incbounds = TRUE)],1)),by=c('patient_num','days_since_admission')][]
  
  # At this point, our table has these headers:
  # patient_num  siteid  days_since_admission  value min_cr_90d min_cr_48h  min_cr_retro_7day min_cr_48h_retro  cr_7d cr_90d
  
  # Now we have to start grading AKI severity at each time point
  # This approach is similar to how the MIMIC-III dataset generates AKI severity
  # Generate two columns using both the formal KDIGO AKI definition and the modified retrospective AKI definition
  message("Generating KDIGO severity grades for each serum Cr value")
  labs_cr_aki$aki_kdigo <- apply(labs_cr_aki,1,aki_kdigo_grade)
  labs_cr_aki$aki_kdigo_retro <- apply(labs_cr_aki,1,aki_kdigo_grade_retro)
  labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_num,days_since_admission) %>% dplyr::mutate(aki_kdigo_final = max(aki_kdigo,aki_kdigo_retro,na.rm=TRUE)) %>% dplyr::ungroup()
  
  # Generate two columns grading AKD severity at 7d and 90d (grade 0B/C is coded as 0.5)
  labs_cr_aki$akd_7d <- apply(labs_cr_aki,1,akd_grade_7d)
  labs_cr_aki$akd_90d <- apply(labs_cr_aki,1,akd_grade_90d)
  
  # Now we are going to generate the start days of each AKI
  labs_cr_aki_tmp <- labs_cr_aki
  labs_cr_aki_tmp$valid = 1
  
  # Find the day of the minimum Cr used for grading AKIs (taken as baseline)
  message("Now finding the day at which the minimum serum Cr is achieved")
  labs_cr_aki_tmp <- labs_cr_aki_tmp %>% dplyr::group_by(patient_num) %>% tidyr::complete(days_since_admission = tidyr::full_seq(days_since_admission,1)) %>% dplyr::mutate(value = zoo::na.fill(value,Inf)) %>% dplyr::ungroup()
  labs_cr_aki_tmp2 <- labs_cr_aki_tmp
  labs_cr_aki_tmp3 <- labs_cr_aki_tmp
  labs_cr_aki_tmp2 <- labs_cr_aki_tmp2 %>% split(.$patient_num) %>% purrr::map(~pos_min(.$value,.$days_since_admission)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_num')
  colnames(labs_cr_aki_tmp2)[2] <- "day_min"
  labs_cr_aki_tmp3 <- labs_cr_aki_tmp3 %>% split(.$patient_num) %>% purrr::map(~pos_min(.$value,.$days_since_admission,lag=FALSE)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_num')
  colnames(labs_cr_aki_tmp3)[2] <- "day_min_retro"
  labs_cr_aki_tmp4 <- cbind(labs_cr_aki_tmp,"day_min" = labs_cr_aki_tmp2$day_min,"day_min_retro" = labs_cr_aki_tmp3$day_min_retro)
  labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[!is.na(labs_cr_aki_tmp4$valid),]
  
  # Generate delta_cr
  message("Now identifying maxima points of serum Cr")
  labs_cr_aki_tmp4 <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_num,days_since_admission) %>% dplyr::mutate(min_cr_7d_final = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::mutate(delta_cr = value - min_cr_7d_final) %>% dplyr::ungroup()
  
  # Use the largest delta_cr to find the peak of each AKI
  labs_cr_aki_delta_maxima <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_num) %>% dplyr::filter(delta_cr %in% delta_cr[which.peaks(delta_cr,decreasing=FALSE)]) %>% dplyr::ungroup()
  labs_cr_aki_delta_maxima$delta_is_max = 1
  labs_cr_aki_delta_maxima <- labs_cr_aki_delta_maxima %>% dplyr::rename(delta_maxima = delta_cr) %>% dplyr::select(patient_num,days_since_admission,delta_maxima,delta_is_max)
  labs_cr_aki_tmp4 <- merge(labs_cr_aki_tmp4,labs_cr_aki_delta_maxima,by=c("patient_num","days_since_admission"),all.x=TRUE)
  
  # Filter for KDIGO grades > 0
  message("Now generating tables of all AKI events")
  labs_cr_aki_tmp5 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final > 0,]
  labs_cr_aki_tmp5[is.na(labs_cr_aki_tmp5)] <- 0
  # Filter for maxima of delta_cr (which should give us the peaks)
  labs_cr_aki_tmp5 <- labs_cr_aki_tmp5[labs_cr_aki_tmp5$delta_is_max > 0,]
  
  # Filter and reorder columns to generate our final table of all AKI events
  aki_peak_summ <- labs_cr_aki_tmp5 %>% dplyr::select(patient_num,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
  
  aki_peak_summ <- aki_peak_summ %>% dplyr::distinct(patient_num,days_since_admission,.keep_all=TRUE)
  aki_peak_summ <- aki_peak_summ[aki_peak_summ$days_since_admission >= 0,]
  
  # Final headers for aki_peak_summ:
  # patient_num,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
  # days_since_admission - time at which peak Cr is achieved
  # day_min - time at which Cr begins to rise
  
  # Generate a separate table (for reference) of all creatinine peaks not fulfilling KDIGO AKI criteria
  message("Now generating a separate table for serum Cr peaks which do not reach AKI definitions")
  labs_cr_nonaki <- labs_cr_aki_tmp4[!(labs_cr_aki_tmp4$patient_num %in% aki_peak_summ$patient_num),]
  labs_cr_nonaki[is.na(labs_cr_nonaki)] <- 0
  labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$delta_is_max > 0,]
  labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_num,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
  labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$days_since_admission >= 0,]
  # Generate the highest Cr peak for non-AKI peaks detected
  nonaki_peak_summ <- labs_cr_nonaki %>% dplyr::group_by(patient_num) %>% dplyr::slice(which.max(delta_cr)) %>% dplyr::ungroup()
  return(list(aki_peak_summ,nonaki_peak_summ))
}