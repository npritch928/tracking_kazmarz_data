library(xtable)
library(tidyverse)
#Functions in this file are 
# read_it - Function to read in results files
# run_output - Function to create a single dataframe given a directory
# moving_av_phase - Runs a two stage moving average
# estimators - calculates the true and estimated moving average for a dataframe
# confidence - vectorised function to compute credible intervals
# confidence_add - function that adds credible interval to dataframe 


# Function to read individual datafiles 

read_it <-function(file,direct){
  parts <- as.vector(str_split(file,"_",simplify = T))
  if (parts[4] == "Achlio") {
    parts <- parts[-5]
  }
  tib_out <- read_csv(str_c(direct,file)) %>%
    mutate(Iteration = 1:nrow(.),
           Matrix_Type = parts[1],
           Rows = as.numeric(parts[2]),
           Cols = as.numeric(parts[3]),
           Sampling_Method = parts[4],
           Sample_size = as.numeric(parts[5]),
           Width = as.numeric(str_split(parts[6],"\\.",simplify = T)[1]))%>%
    select(Iteration,Matrix_Type,Rows,Cols,Sampling_Method,Sample_size,Width,everything())
  return(tib_out)
}

read_it_rep <- function(file,direct){
 parts <- as.vector(str_split(file,"_|\\.",simplify = T))
  if (parts[4] == "Achlio") {
    parts <- parts[-5]
  }

  tib_out <- read_csv(str_c(direct,file)) %>%
    mutate(Iteration = 1:nrow(.),
           Replication = parts[7],
           Matrix_Type = parts[1],
           Rows = as.numeric(parts[2]),
           Cols = as.numeric(parts[3]),
           Sampling_Method = parts[4],
           Sample_size = as.numeric(parts[5]),
           Width = as.numeric(str_split(parts[6],"\\.",simplify = T)[1]))%>%
    select(Iteration,Matrix_Type,Rows,Cols,Sampling_Method,Sample_size,Width,everything()) %>%
    mutate(Obs_Res = Obs_Res * case_when(Sampling_Method == "Uniform"~ Rows / Sample_size,TRUE~1))
  return(tib_out)
} 

confidence_big<-function(rho, iota, var, width, smin, alpha,ncols,  eta = 1){
  
  #tstar <-(log(2/(1-alpha)) * iota * omega * 2 * c)
  #tstart <- ifelse( tstar > sqrt(tstar/omega), tstar, sqrt(tstar/omega))
  #tstart <- sqrt(2 * log(2/alpha) * iota * (9 * ncols)^2 / (4 * smin^2) *(1+log(width))/(eta*width))
  tstart <- sqrt(2 * log(2/alpha) * iota * var *(1+log(width))/(eta*width))
  up <- rho + tstart
  low <- rho - tstart
  #low <- ifelse(low <= 0, min(1e-16,min(rho)*10^-2), low)
  return(tibble("eta" = eta, Lower = low, Upper = up, sigma = sqrt( (9 * ncols)^2 / (4 * smin^2) *(1+log(width))/(eta*width))))
}

stopping_big <- function(ncols, ms, var, width, iota, nu, delt1, delt2, xi1, xi2, eta=1){
  #c <-  2 * ms^2 * width *eta  / ( (9 * ncols)^2 * (1 + log(width)))   
  c <-  2 * width *eta  / ( (var) * (1 + log(width)))  
  t1 <- c*(delt1 - 1)^2 * nu^2/(sqrt(iota)*log(1/xi1))
  t2 <- c*(1 - delt2)^2 * nu^2/(sqrt(iota)*log(1/xi2))
  thres <- ifelse(t1 > t2, t2, t1)
  #stopped <- ifelse(iota <= thres, T, F)
  return(thres)
}

read_big <- function(filena, srow, scol, smin,var,wid, alpha = .95, nu = 100, delt1 = 1.1, delt2 = .9, xi1=.01, xi2=.01, eta=1){
  gou <- read_csv(filena) %>%
    mutate(ncol = scol, nrow = srow)%>%
    mutate(it = rep(1:(nrow(.)%/%4),4)) %>%
    pivot_wider(names_from = "type", values_from = "value") %>%
    new_mov_av(wid)
  #CI <- conf_calc(gou$Mov_Obs_Res,gou$Mov_Obs_Res2,vari,gou$Width,alpha)
return(
  gou %>%
    #mutate(Upper= CI$Upper,Lower= CI$Lower)%>%#confidence_big(Mov_Obs_Res,Mov_Obs_Res2,var,Width,smin,alpha,scol,eta=eta)) %>%
    mutate(threshold = stopping_big(scol,smin,var,Width,Mov_Obs_Res2,nu, delt1, delt2, xi1, xi2,eta)) %>%
    mutate(nu = nu, delt1 = delt1, delt2 = delt2, xi1 = xi1, xi2 = xi2, eta = eta, stop = rho < nu)
  )
}

conf_calc<-function(rho,iota,var,width,alpha,eta=1){
  tstart <- sqrt(2 * log(2/alpha) * iota * var *(1+log(width))/(eta*width))
  up <- rho + tstart
  low <- rho - tstart
  return(tibble(Lower = low, Upper = up))
}

stopping_graph_big<-function(results){
  name <- str_c("Error from stopping criterion")
  
  g <- results %>%
    filter(sqrt(iota) <= threshold) %>%
    summarise(n = n(), 
              t1 = sum(stop & rho > nu * delt1), 
              t2 = sum(!stop & rho <= nu * delt2), 
              correct = sum(stop & rho <= nu * delt1 | !stop & rho > nu * delt2))%>%
    ungroup()%>%
    summarise(total = sum(n), 
              t1 = sum(t1)/total, 
              correct = sum(correct)/total, 
              t2 = sum(t2)/total) 
  return(g)
}

#Function that creates a full dataset given a directory
#Used for experiments 1-7
#Assumes that director has only csv files inside of it
run_output <- function(direct){
  if(!str_detect(direct,"/$")){
    direct <- str_c(direct,"/")
  }
  files <- dir(direct)
  data_list<- map(files, read_it, direct)
  return(bind_rows(data_list))
}
run_output_rep <- function(direct){
  if(!str_detect(direct,"/$")){
    direct <- str_c(direct,"/")
  }
  files <- dir(direct)
  data_list<- map(files, ~read_it_rep(direct,.x))
  return(bind_rows(data_list))
}
#Moving av two phases operates on assumption that in convergence stage mean 
#of miving average greater than true observations
Moving_av_phase<-function(width,values,i,group_it,start)
{
  if(is.na(start) & values[i] <= values[ifelse(i - group_it + 1 == 1, i ,i-1)] ){
    
    start <- NA
    return(c(start,values[i],1))
    
  }else if(is.na(start) & values[i] > values[ifelse(i - group_it + 1 == 1, i ,i-1)]){
    
    start <- i
    return(c(start,values[i],1))
    
  }else{

    if(i - start >= width ){
      
      me<-mean(values[(i - width + 1):i])
      return(c(start,me, width))
      
    }else{
      
      me<-mean(values[start:i])
      return(c(start,me,i - start + 1))
      
    }
    
  }
  
}
Moving_max_phase<-function(width,values,i,group_it,start)
{
  if(is.na(start) & values[i] <= values[ifelse(i - group_it + 1 == 1, i ,i-1)] ){
    
    start <- NA
    return(c(start,values[i],1))
    
  }else if(is.na(start) & values[i] > values[ifelse(i - group_it + 1 == 1, i ,i-1)]){
    
    start <- i
    return(c(start,values[i],1))
    
  }else{
    
    if(i - start >= width ){
      
      me<-max(values[(i - width + 1):i])
      return(c(start,me, width))
      
    }else{
      
      me<-max(values[start:i])
      return(c(start,me,i - start + 1))
      
    }
    
  }
  
}

#Moving Average calculation
Moving_av <- function(width,values,i)
{
  if(i > width){
    mean(values[(i - width + 1):i])
  }else{
    mean(values[1:i])
  }
}
#Function that gets the delayed error 
error_conf <- function(error,iteration,width){
  len <- length(error)
  new_err <- numeric(len)
  for(i in 1:len){
    j <- iteration[i]
    w <- width[i]
    new_err[i] <- error[j - w + 1]
  }
  return(new_err)
}

#function that calculates rho_k and iot_k at every 
#iteration as well as their true values
estimators<-function(results,width){
  
  widths <- unique(results$Width)
  
  if(length(widths) > 1 ){
    results <- results %>%
      filter(Width == widths[1])
  }
  if (sum("Replication" %in% colnames(results))>=1){
    results <- results %>% 
      arrange(Matrix_Type,Sampling_Method,Sample_size, Rows, Replication, Iteration)
  }else if(sum("Matrix" %in% colnames(results))==1){
    results <- results %>% 
      arrange(Matrix_Type, Iteration)
  }else if(sum("Matrix" %in% colnames(results))==0){
    results <- results %>% 
      arrange(Matrix_Type, Iteration)
  }else{
   results <- results %>% 
      arrange(Matrix_Type,Sampling_Method,Sample_size, Rows, Iteration) 
  }
  
  len <- nrow(results)
  #Define the residual values that will be used
  Obs_res <- results$Obs_Res
  True_res <- results$True_Res
  Obs_res2 <- Obs_res^2
  True_res2 <- True_res^2
  #Assumes every expirement runs for same num of iterations
  iterations <- results$Iteration
  #Define vectors for saving rhos and iotas and true values
  Mov_Obs_res <- numeric(len)
  Mov_True_res <- numeric(len)
  Mov_Obs_res2 <- numeric(len)
  Mov_True_res2 <- numeric(len)
  Widths <- numeric(len)
  
  rhos_obs <- c(NA,NA,NA)
  rhos_true <- c(NA,NA,NA)
  iota_obs <- c(NA,NA,NA)
  iota_true <- c(NA,NA,NA)
  
  for (i in 1:len) {
    if(iterations[i] == 1){
      j <- i
      rhos_obs[1] <- NA
      iota_obs[1] <- NA
    }
    rhos_obs <- Moving_av_phase(width,Obs_res,i,j,rhos_obs[1])
    rhos_true <- Moving_av_phase(width,True_res,i,j,rhos_obs[1])
    iota_obs <- Moving_av_phase(width,Obs_res2,i,j,iota_obs[1])
    iota_true <- Moving_av_phase(width,True_res2,i,j,iota_obs[1])
    #rhos_obs <- Moving_av(width,Obs_res,i)
    #rhos_true <- Moving_av(width,True_res,i)
    #iota_obs <- Moving_av(width,Obs_res2,i)
    #iota_true <- Moving_av(width,True_res,i)
    Mov_Obs_res[i] <- rhos_obs[2]
    Mov_True_res[i] <- rhos_true[2]
    Mov_Obs_res2[i] <- iota_obs[2]
    Mov_True_res2[i] <- iota_true[2]
    Widths[i] <- rhos_obs[3]
  }
  return( results%>%
            mutate(Width = Widths, 
               Mov_Obs_Res = Mov_Obs_res, 
               Mov_True_Res = Mov_True_res, 
               Mov_Obs_Res2 = Mov_Obs_res2,
               Mov_True_Res2 = Mov_True_res2)
        )
}

estimatorsBig<-function(results,width){
  
  widths <- unique(results$Width)
  
  results <- results %>% 
      arrange(Iteration)
  
  len <- nrow(results)
  #Define the residual values that will be used
  Obs_res <- results$rho
  Obs_res2 <- results$iota
  #Assumes every expirement runs for same num of iterations
  iterations <- results$Iteration
  #Define vectors for saving rhos and iotas and true values
  Mov_Obs_res <- numeric(len)
  Mov_Obs_res2 <- numeric(len)
  Widths <- numeric(len)
  
  rhos_obs <- c(NA,NA,NA)
  iota_obs <- c(NA,NA,NA)
  
  for (i in 1:len) {
    if(iterations[i] == 1){
      j <- i
      rhos_obs[1] <- NA
      iota_obs[1] <- NA
    }
    rhos_obs <- Moving_av_phase(width,Obs_res,i,j,rhos_obs[1])
    iota_obs <- Moving_av_phase(width,Obs_res2,i,j,iota_obs[1])
    #rhos_obs <- Moving_av(width,Obs_res,i)
    #rhos_true <- Moving_av(width,True_res,i)
    #iota_obs <- Moving_av(width,Obs_res2,i)
    #iota_true <- Moving_av(width,True_res,i)
    Mov_Obs_res[i] <- rhos_obs[2]
    Mov_Obs_res2[i] <- iota_obs[2]
    Widths[i] <- rhos_obs[3]
  }
  return( results%>%
            mutate(Width = Widths, 
                   rho = Mov_Obs_res, 
                   iota = Mov_Obs_res2)
  )
}

#Function that computes the confidence bounds 
confidence<-function(rho, iota, width, Sampling_meth, Sample_size, alpha, etav, nrows=256, sigma2= 0, mega = 0){
  eta = etav#ifelse(etav == F, case_when(Sampling_meth == "Achlio" | Sampling_meth == "gauss" ~2,#13
                                   #Sampling_meth == "FJLT" ~ 2),etav) #188
  if(max(sigma2) == 0){
    c <- case_when(Sampling_meth == "Achlio" | Sampling_meth == "gauss"~ (1+log(width)) / (Sample_size*width*1.1*eta),#.23467),
                 Sampling_meth == "FJLT" ~ (1 + log(width))/ (Sample_size*width*.03125*eta),
              TRUE ~ ((1+log(width))*nrows^2/(4*Sample_size^2))  / (width*eta))
    #
    omega <- case_when(Sampling_meth == "Achlio" | Sampling_meth == "gauss"~.46,#.1127,
                  Sampling_meth == "FJLT" ~ 1/16,
                  TRUE ~ 0)
  }else{
    c <- (sigma2 * (1+log(width))) / (width*eta)
    omega <- mega
  }
  #tstar <-(log(2/(1-alpha)) * iota * omega * 2 * c)
  #tstart <- ifelse( tstar > sqrt(tstar/omega), tstar, sqrt(tstar/omega))
  tstart <-ifelse(sqrt(2* log(2 / (1-alpha))* c * iota)> sqrt(iota)*2*log(2/(1-alpha))* omega /( eta *width), sqrt(2* log(2 / (1-alpha))* c * iota), sqrt(iota)*2*log(2/(1-alpha))* omega / (eta*width))
  up <- rho + tstart
  low <- rho - tstart
  #low <- ifelse(low <= 0, min(1e-16,min(rho)*10^-2), low)
  return(tibble("eta" = eta, Lower = low, Upper = up))
}

#functions that adds confidence bounds to the dataframe assumes you ran estimators first
confidence_add <- function(results,alpha,etav = 1, sigma2 = sigma2){
  if(!"Mov_Obs_Res2" %in% colnames(results)){
    return(cat("ERROR:Make sure you run the function estimators() first."))
  }
  etav <- rep(etav,nrow(results))
  if (sum("Replication" %in% colnames(results))>=1){
    results <- results %>% 
      arrange(Matrix_Type,Sampling_Method,Sample_size, Rows, Replication, Iteration)
  }else{
   results <- results %>% 
      arrange(Matrix_Type,Sampling_Method,Sample_size, Rows, Iteration) 
  }
  nr <- max(results$Rows)
  lowup_est <- confidence(results$Mov_Obs_Res,
             results$Mov_Obs_Res2,
             results$Width,
             results$Sampling_Method,
             results$Sample_size,
             alpha,
             etav, 
             nrows = nr,
             sigma2 = sigma2)
  er <- results$Error #error_conf(results$Error,results$Iteration,results$Width)
  if(sum("A_Norm" %in% colnames(results))>=1){
    lowup_true <- confidence( results$Mov_Obs_Res,
                               er^2 * results$A_Norm^4,
                               results$Width,
                               results$Sampling_Method,
                               results$Sample_size,
                               alpha,
                               etav,
                               nrows = nr)
    results <- results %>%
      bind_cols(lowup_true) %>%
      mutate(Lower_true = Lower, Upper_true = Upper)%>%
      select(-Upper,-Lower,-eta)
  }
  return(
    bind_cols(results, lowup_est)%>%
      mutate(Lower_est = Lower, Upper_est = Upper)%>%
      select(-Upper,-Lower)
  )
  
}

#Function designed to indicate when you should stop
stopping <- function(Sampling_meth, Sample_size, width, iota, nu, delt1, delt2, xi1, xi2, eta,nrows=256, sigma2= 0, mega = 0){
  if(max(sigma2) == 0){
    c <- case_when(Sampling_meth == "Achlio" | Sampling_meth == "gauss"~
              (eta*Sample_size*width * 0.23467)/(1+log(width)) ,
              Sampling_meth == "FJLT" ~ (eta*Sample_size*width *.03125)/(1 + log(width)),
              Sampling_meth == "Uniform"~ (eta * width * (4* Sample_size^2)/nrows^2)/(1 + log(width)))
  #
  #256/(2* Sample_size)
  omega <- case_when(Sampling_meth == "Achlio" | Sampling_meth == "gauss"~
                  .1127,
                  Sampling_meth == "FJLT" ~ 1/16,
                  Sampling_meth == "Uniform" ~ 0)
  }else{
    c<- (eta * width)/ (sigma2 * (1 + log(width)))
    omega <- mega
  }
  #t1 <- c * min((delt1 - 1)^2 * nu^2, (delt1 - 1) * nu / omega)/log(2/xi1)
  #t2 <- c * min((1 - delt2)^2 * nu^2, (1 - delt2) * nu / omega)/log(2/xi2)
  a1 <- c*(delt1 - 1)^2 * nu^2/(2 *sqrt(iota)*log(1/(xi1)))
  b1 <- (delt1 - 1) * nu /(2*log(1/(xi1))*omega)
  a2 <- c*(1 - delt2)^2 * nu^2/(2*sqrt(iota)*log(1/xi2))
  b2 <- (1-delt2) * nu / (2*log(1/xi2)*omega)
  t1 <- ifelse(a1 < b1, a1, b1)
  t2 <- ifelse(a2 < b2, a2, b2)
  #t1 <- (min(c*(delt1 - 1)^2 * nu^2/(sqrt(iota)*log(2/xi1)), (delt1 - 1) * nu / (log(2/xi1)*omega)))^2
  #t2 <- (min(c*(1 - delt2)^2 * nu^2/(sqrt(iota)*log(2/xi2)), (1-delt2) * nu / (log(2/xi2)*omega)))^2
  thres <- ifelse(t1 > t2, t2, t1)
  return(thres)
}

#Function that adds columns indicating according to stopping criteria where
#the algorithm should be stopped versus the true average residuals

stopping_add <- function(results, nuv, delt1v, delt2v, xi1v, xi2v){
  nr <- max(results$Rows)
  st <- results %>% 
    mutate(stop = Mov_Obs_Res <= nuv, 
           threshold = stopping(Sampling_Method, Sample_size, Width, Mov_Obs_Res2, nuv, delt1v, delt2v, xi1v, xi2v, eta, nrows = nr),
           nu = nuv, 
           delt1 = delt1v, 
           delt2 = delt2v, 
           xi1 = xi1v,
           xi2 = xi2v)
  
  return(st)
}
