#Function that reads in the initial dataset and adds required information
read_it_rep <- function(file, direct){
  parts <- as.vector(str_split(file,"_",simplify = T))
  if (parts[4] == "Achlio") {
    parts <- parts[-5]
  }
  tib_out <- read_csv(str_c(direct,file)) %>%
    mutate(Iteration = it,
           #Most basic matrix info is in the filename
           Matrix_Type = parts[1],
           Rows = as.numeric(parts[2]),
           Cols = as.numeric(parts[3]),
           Sampling_Method = parts[4],
           Sample_size = as.numeric(parts[5]),
           Width = as.numeric(str_split(parts[6],"\\.",simplify = T)[1]),
           #accidentally recorded norm of values not norm squared so correct
           Obs_Res = Obs_Res^2, 
           True_Res = True_Res^2
           )%>%
    select(Iteration, Matrix_Type, Rows, Cols, Sampling_Method, Sample_size, Width, everything()) %>%
    mutate(Obs_Res = Obs_Res * case_when(Sampling_Method == "Uniform"~ Rows / Sample_size, TRUE~1),
           Obs_Res2 = Obs_Res^2,
           True_Res2 = True_Res^2)
  return(tib_out)
}
#Get the confidence intervals
confidence_rep<- function(results, width, alpha, eta=1, sigma2=0){
  #Get 1st repetition to compute the rhos
  rep1 <- results %>%
    filter(reps == 1, ob == 1)
  l = nrow(rep1)
  Mov_Obs_res = numeric(l)
  Mov_True_res = numeric(l)
  Mov_Obs_res2 = numeric(l)
  for(i in 1:l){
    Mov_Obs_res[i] = Moving_av(width, rep1$Obs_Res, i)
    Mov_True_res[i] = Moving_av(width, rep1$True_Res, i)
    Mov_Obs_res2[i] = Moving_av(width, rep1$Obs_Res2, i)
  }
  rep1 <- rep1 %>%
    mutate( Mov_Obs_Res = Mov_Obs_res,
            Mov_True_Res = Mov_True_res,
            Mov_Obs_Res2 = Mov_Obs_res2) %>%
    filter(it >= width) %>% 
    mutate(rel_err = abs(Obs_Res - Mov_Obs_Res)/Mov_Obs_Res, 
           rel_err_t = abs(Obs_Res - True_Res)/True_Res) 
  rep1 <- rep1 %>%
    #estimators(width) %>%
    confidence_add(alpha, eta, sigma2=sigma2) %>% 
    rename(rho = Mov_True_Res, orho = Mov_Obs_Res, oiota = Mov_Obs_Res2) %>%
    select(it, reps, Sampling_Method, Rows, Sample_size, Matrix_Type, Upper_est, Lower_est, rho, orho, oiota, rel_err, rel_err_t) %>%
    rename(Iteration = it)
  
  #Define the rho 
  #rep1 <-  %>%
    #select(Iteration, Matrix_Type, Sampling_Method, Mov_True_Res) %>%
    
  return(
    gt<-results %>% 
      group_by(Matrix_Type, Sampling_Method, Rows, Sample_size, it, reps) %>%
      #Take mean over the columns to get the moving averages
      summarise(Mov_True_Res = mean(True_Res), 
                Mov_Obs_Res = mean(Obs_Res),
                Mov_Obs_Res2 = mean(Obs_Res2),
                Error = mean(Error)) %>%
      rename(Iteration = it) %>%
      mutate(Width = width) %>%
      #Shift to account for the fact that the observations at I1 are actually i15
      mutate(Iteration = Iteration + width - 1) %>%
      filter(Iteration <= 300) %>%
      left_join(rep1, by = c("Iteration","reps","Sampling_Method", "Rows", "Sample_size","Matrix_Type")) %>%
      group_by(Sampling_Method, Matrix_Type, Iteration, Sample_size) %>%
      #expand the estimates to all observations
      mutate(orho = max(orho, na.rm=T),
             rho = max(rho, na.rm=T),
             oiota = max(oiota, na.rm=T),
             sigma2 = max(sigma2, na.rm = T),
             rel_err = max(rel_err, na.rm = T),
             rel_err_t = max(rel_err_t, na.rm = T),
             Upper_est = max(Upper_est, na.rm = T),
             Lower_est = max(Lower_est, na.rm = T)) 
    
  )
}

stopping_add_rep <- function(results, nuv, delt1v, delt2v, xi1v, xi2v, eta, sigma2=sigma2){
  if(!("obs_res" %in%colnames(results))){
    results <- results %>%
      mutate(obs_res = Mov_Obs_Res)
  }
  return(
    results %>% 
      mutate(eta = eta) %>%
      mutate(stop = obs_res <= nuv, 
             threshold = stopping(Sampling_Method, Sample_size, Width, Mov_Obs_Res2, nuv, delt1v, delt2v, xi1v, xi2v, eta, sigma2 = sigma2, nrow = Rows),
             nu = nuv, 
             delt1 = delt1v, 
             delt2 = delt2v, 
             xi1 = xi1v,
             xi2 = xi2v)
  )
}
#Function adding interval and criterion to all entries of the repeated samples
get_es_ci_st <- function(file,direct,width = 30, alph = .05, nuv=1, delt1v=1.1, delt2v=.9, xi1v=.05, xi2v=.05,eta=1, Sigma2=0){
  return(
    read_it_rep(file, direct) %>%
      confidence_rep(width,alph,sigma2=Sigma2,eta = eta) %>% 
      rename(obs_res = Mov_Obs_Res) %>% 
      stopping_add_rep(nuv, delt1v, delt2v, xi1v, xi2v,eta,sigma2=Sigma2) %>%
      rename(Mov_Obs_Res = obs_res)
  )
}

#Function that reads in all the datasets of the repeated samples
read_all_rep <- function(direct,width = 30, alph = .95, nuv=1e-5, delt1v=1.1, delt2v=.9, xi1v=.05, xi2v=.05,eta=1, sigma2 = 0){
  if(!str_detect(direct,"/$")){
    direct <- str_c(direct,"/")
  }
  files <- dir(direct)
  data_list<- map(files, get_es_ci_st, direct, width, alph, nuv, delt1v, delt2v, xi1v, xi2v, eta, Sigma2=sigma2)
  return(bind_rows(data_list))
}

#function to measure stopping errors
stopping_graph_rep<-function(results){
  name <- str_c("Error from stopping criterion")
  
  g <- results %>%
    #filter(sqrt(Mov_Obs_Res2) <= threshold) %>%
    group_by(Matrix_Type, Sampling_Method, Sample_size) %>%
    summarise(n = n(), 
              #Compute Late stopping
              t2 = sum(Mov_Obs_Res <= nu & Mov_True_Res > nu * delt1 & sqrt(Mov_Obs_Res2) <= threshold), 
              n2 = sum(Mov_True_Res > nu * delt1),
              #Compute early stopping
              t1 = sum(Mov_Obs_Res > nu & Mov_True_Res <= nu * delt2 & sqrt(Mov_Obs_Res2) <= threshold),
              n1 = sum(Mov_True_Res <= nu * delt2), 
              #Compute power i.e. probability you don't stop when you are supposed to stop
              p = sum(Mov_Obs_Res < nu & Mov_True_Res < nu & sqrt(Mov_Obs_Res2) > threshold),
               n3 = sum(Mov_Obs_Res < nu & Mov_True_Res < nu )) #%>%
    # ungroup()%>%
    # summarise(total = sum(n),
    #           n1 = sum(n1),
    #           t1 = sum(t1),
    #           r1 = t1/n1,
    #           n2 = sum(n2),
    #           t2 = sum(t2),
    #           r2 = t2/n2,
    #           n3 = sum(n3),
    #           p = sum(p),
    #           rp = p/n3)
  return(g)
}
#COmpute stopping failures without mov_obs_res2 threshold
stopping_graph_rep_nc<-function(results){
  name <- str_c("Error from stopping criterion")
  
  g <- results %>%
    group_by(Matrix_Type, Sampling_Method, Sample_size) %>%
    summarise(n = n(), 
              #Count how often you stop too early
              t1 = sum(stop & Mov_True_Res > nu * delt1), 
              #Total number of times you could stop too late
              n1 = sum(Mov_True_Res > nu * delt1),
              #Count how often you stop too late
              t2 = sum(!stop & Mov_True_Res <= nu * delt2),#, 
              #Total number of times you could stop too late
              n2 = sum(Mov_True_Res <= nu * delt2),
              correct = sum(stop & Mov_True_Res <= nu * delt1 | !stop & Mov_True_Res > nu * delt2))%>%view()
    # ungroup()%>%
    # summarise(total = sum(n), 
    #           n1 = sum(n1),
    #           t1 = sum(t1),
    #           r1 = t1/n1,
    #           n2 = sum(n2),
    #           t2 = sum(t2),
    #           r2 = t2/n2, 
    #           correct = 1 - r2 - r1) 
  return(g)
}

