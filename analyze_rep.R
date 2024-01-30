library(xtable)
library(tidyverse)
source("analyze.R")
#Functions in this file are 
# read_it_rep - Function to read in replicated results files
# confidence_rep - Function that adds credible interval to dataframe 
# get_all_obs - Function that gets moving average of each replication
# and that corresponds to the original interval

# Function to read individual datafiles 

read_it_rep <-function(direct,file){
  parts <- as.vector(str_split(file,"_",simplify = T))
  if (parts[4] == "Achlio") {
    parts <- parts[-5]
  }
  tib_out <- read_csv(str_c(direct,file)) %>%
    mutate(Iteration = it,
           Matrix_Type = parts[1],
           Rows = as.numeric(parts[2]),
           Cols = as.numeric(parts[3]),
           Sampling_Method = parts[4],
           Sample_size = as.numeric(parts[5]),
           Width = as.numeric(str_split(parts[6],"\\.",simplify = T)[1]),
           Obs_Res = Obs_Res^2, True_Res = True_Res^2)%>%
    select(Iteration,Matrix_Type,Rows,Cols,Sampling_Method,Sample_size,Width,everything()) %>%
    mutate(Obs_Res = Obs_Res * case_when(Sampling_Method == "Uniform"~ Rows / Sample_size, TRUE~1),
           Obs_Res2 = Obs_Res^2)
  return(tib_out)
}

get_rhos <- function(file,direct){
  tib <- read_it_rep(direct,file)
  tib1 <- confidence_rep(tib,tib$Width[1],.95)
  tib <- get_all_obs_rho(tib,tib1,15)
  return(
    tib %>% 
      group_by(Matrix_Type,Sampling_Method,it) %>%
      #summarise(rho = mean(Mean_Obs_Res), obs_res = max(Mov_Obs_Res)) %>%
      mutate(error = abs(rho - obs_res), rel_err = abs(rho - obs_res)/rho) %>%
      rename(Iteration = it)%>%
      drop_na()
      #group_by(Matrix_Type,Sampling_Method,it,reps) %>%
      #mutate(rho = mean(Obs_Res))%>%
      #summarise(rho = mean(rho), obs_res = max(Mov_Obs_Res)) %>%
      #ungroup() %>%
      #group_by(Matrix_Type,Sampling_Method,it) %>%
      #summarise(rho = mean(rho), obs_res = max(obs_res))
  )
}

run_rhos <- function(direct){
  if(!str_detect(direct,"/$")){
    direct <- str_c(direct,"/")
  }
  files <- dir(direct)
  data_list<- map(files, get_rhos, direct)
  return(bind_rows(data_list))
}
# Function that returns confidence interval for 1st rrep. 
confidence_rep<- function(results, width, alpha,eta=1,sigma2=0){
  rep1 <- results %>%
    filter(reps == 1, ob == 1) 
  
  rep1 <- estimators(rep1, width)
  
  rep1 <- confidence_add(rep1, alpha, eta,sigma2=sigma2)
  
  return(rep1)
}

# Function that returns confidence interval for 1st rrep. 
confidence_rep_const<- function(results, width, alpha,eta=1){
  rep1 <- results %>%
    filter(reps == 1, ob == 1) 
  
  rep1 <- estimators(rep1, width)
  
  rep1 <- confidence_add(rep1, alpha,eta)
  
  return(rep1)
}
#Functions that make credible interval at each iteration for each observation
estimators_rep_eo<-function(results){
  return(
    results%>% 
      group_by(Sampling_Method, Matrix_Type,Sample_size, Rows, Width, it,reps) %>% 
      mutate(Mov_Obs_Res = mean(Obs_Res), 
             Mov_True_Res = mean(True_Res),
             sigma = abs(Obs_Res - Mov_Obs_Res) / Mov_Obs_Res,
             sigma_t = abs(Obs_Res - True_Res) / True_Res) %>%
       summarise(Mov_Obs_Res = mean(Obs_Res), 
                 Mov_Obs_Res2 = mean(Obs_Res2), 
                 sigma = (max(sigma) / 2)^2,
                 sigma_t = (max(sigma_t) / 2)^2) %>%
      rename(Iteration = it) 
  )
}

confidence_add_eo <- function(results,alpha,etav = 1, nrow = 256, sigma2=0){
  if("eta" %in% colnames(results)){
    results <- results %>%
      select(-eta)
  }
  etav <- rep(etav,nrow(results))
  results <- results %>% 
      arrange(Matrix_Type,Sampling_Method,Sample_size, Rows, Iteration, reps) 
  lowup_est <- confidence(results$Mov_Obs_Res,
                          results$Mov_Obs_Res2,
                          results$Width,
                          results$Sampling_Method,
                          results$Sample_size,
                          alpha,
                          etav,
                          nrows = nrow,
                          sigma2 = sigma2)
  
  return(
    bind_cols(results, lowup_est)%>%
      mutate(Lower_est = Lower, Upper_est = Upper, Iteration = Iteration + Width - 1)%>%
      select(-Upper,-Lower)
  )
  
}


confidence_rep_eo<- function(results, width, alpha,eta=1,sigma2=0){
  
  rep1 <- estimators_rep_eo(results) 
  Nrows <- max(rep1$Rows)
  rep1 <- confidence_add_eo(rep1, alpha, eta, nrow=Nrows, sigma2= sigma2)
  rhos <- results %>%
    group_by(it, ob) %>%
    summarise(Mean_Obs_Res = mean(Obs_Res)) %>%
    ungroup() %>%
    group_by(it) %>%
    summarise(rho = mean(Mean_Obs_Res)) %>%
    rename(Iteration = it) %>%
    mutate(Iteration = Iteration + width - 1)
  
  
  
  return(left_join(rhos, rep1, by = c("Iteration")))
}


get_es_ci_st <- function(file,direct,width = 15, alph = .95, nuv=1, delt1v=1.1, delt2v=.9, xi1v=.01, xi2v=.01,eta=1, sigma2=0){
  return(
    read_it_rep(direct,file) %>%
      confidence_rep_eo(width,alph,sigma2=sigma2) %>% 
      rename(obs_res = Mov_Obs_Res) %>% 
      stopping_add_rep(nuv, delt1v, delt2v, xi1v, xi2v,eta,sigma2=sigma2) %>%
      rename(Mov_Obs_Res = obs_res) %>% view()
  )
}

read_all_rep <- function(direct,width = 15, alph = .95, nuv=1e-5, delt1v=1.1, delt2v=.9, xi1v=.01, xi2v=.01,eta=1, sigma2 = 0){
  if(!str_detect(direct,"/$")){
    direct <- str_c(direct,"/")
  }
  files <- dir(direct)
  data_list<- map(files, get_es_ci_st, direct, width, alph, nuv, delt1v, delt2v, xi1v, xi2v, eta, sigma2=sigma2)
  return(bind_rows(data_list))
}

get_all_obs <- function(results, rep1, width){
  widths <- rep1 %>%
    filter(Width < width, it >= width) %>%
    pull(Width)
  results <- results %>%
    mutate(it  = it + width - 1)
  k <- width
  for(i in widths){
    results <- results %>%
      filter(it == k & ob > k - i & ob <= k | it != k)
    k <- k + 1
  }
  set_to_join <- rep1 %>%
    select(it,Matrix_Type,Sampling_Method,Sample_size,Width,Mov_Obs_Res,Lower_est,Upper_est) 
  avs <- results %>% 
    filter(it <= 500)%>%
    group_by(it, reps) %>%
    summarise(Mean_Obs_Res = mean(Obs_Res), Mean_True_Res = mean(True_Res),Mean_Obs_Res2 = mean(Obs_Res^2), Mean_True_Res2 = mean(True_Res^2))
  return(left_join(avs, set_to_join, by = c("it")) %>%
    mutate(Bounds = ifelse(Mean_True_Res > Upper_est | Mean_True_Res < Lower_est, "yes", "no"),
           Upperd = Upper_est - Mov_Obs_Res,
           Lowerd = Lower_est - Mov_Obs_Res,
           Meand = Mean_True_Res - Mov_Obs_Res)
  )
}

get_all_obs_us <- function(results, rep1, width){
  results <- results %>%
    mutate(grouping = reps %% 1) %>%
    group_by(it, ob, grouping) %>%
    summarise(Mean_Obs_Res = mean(Obs_Res))
  widths <- rep1 %>%
    filter(Width < width, it >= width) %>%
    pull(Width)
  results <- results %>%
    mutate(it  = it + width - 1) %>% 
    filter(it <= 300) %>%
    group_by(it,grouping) %>% 
    summarise(rho = mean(Mean_Obs_Res))
  set_to_join <- rep1 %>%
    select(it,Matrix_Type,Sampling_Method,Sample_size,Width,Mov_Obs_Res,Mov_Obs_Res2,Lower_est,Upper_est) 
  
  return(left_join(results, set_to_join, by = c("it"))%>%
           mutate(Bounds = ifelse(rho > Upper_est | rho < Lower_est, "yes", "no"),
                  Upperd = Upper_est - Mov_Obs_Res,
                  Lowerd = Lower_est - Mov_Obs_Res,
                  Meand = rho - Mov_Obs_Res)
  )
}

stopping_graph_rep<-function(results){
  name <- str_c("Error from stopping criterion")
  
  g <- results %>%
    #filter(sqrt(Mov_Obs_Res2) <= threshold) %>%
    group_by(Matrix_Type, Sampling_Method, Sample_size) %>%
    summarise(n = n(), 
              t1 = sum(Mov_Obs_Res < nu & Mov_True_Res > nu * delt1 & sqrt(Mov_Obs_Res2) <= threshold), 
              n1 = sum(Mov_True_Res > nu * delt1),# & orho < nu),
              t2 = sum(Mov_Obs_Res > nu & Mov_True_Res <= nu * delt2 & sqrt(Mov_Obs_Res2) <= threshold),
              n2 = sum(Mov_True_Res <= nu * delt2), #& orho > nu),
               correct = sum(stop & rho <= nu * delt1 | !stop & rho > nu * delt2)) %>%
    ungroup()%>%
    summarise(total = sum(n),
              n1 = sum(n1),
              t1 = sum(t1),
              r1 = t1/n1,
              #correct = sum(correct)/total,
              n2 = sum(n2),
              t2 = sum(t2),
              r2 = t2/n2)
   return(g)
}


get_all_obs_rho <- function(results, rep1, width){
  results <- results %>%
    group_by(it, ob) %>%
    summarise(Mean_Obs_Res = mean(Obs_Res))
  widths <- rep1 %>%
    filter(Width < width, it >= width) %>%
    pull(Width)
  results <- results %>%
    mutate(it  = it + width - 1)
  results <- results %>% filter(it <= 300) %>% group_by(it) %>% summarise(rho = mean(Mean_Obs_Res))
  set_to_join <- rep1 %>%
    select(it,Matrix_Type,Sampling_Method,Mov_Obs_Res,Mov_Obs_Res2,Sample_size,Width) %>%
    rename(obs_res = Mov_Obs_Res)
  return(left_join(results, set_to_join, by = c("it"))
  )
}

stopping_add_rep <- function(results, nuv, delt1v, delt2v, xi1v, xi2v,eta, sigma2 = 0){
  if(!("obs_res" %in%colnames(results))){
    results <- results %>%
      mutate(obs_res = Mov_Obs_Res)
  }
  st <- results %>% 
    mutate(eta = eta) %>%
    mutate(stop = obs_res <= nuv, 
           threshold = stopping(Sampling_Method, Sample_size, Width, Mov_Obs_Res2, nuv, delt1v, delt2v, xi1v, xi2v, eta, sigma2 = sigma2),
           nu = nuv, 
           delt1 = delt1v, 
           delt2 = delt2v, 
           xi1 = xi1v,
           xi2 = xi2v)
  
  return(st)
}

make_graphs <- function(results){
  mat_name <- results$Matrix_Type[1]
  samp_name <- results$Sampling_Method[1]
  name <- str_c("Distribution of 1000 repetitions for a\n", mat_name ," matrix")
  g<- ggplot(results,aes(x = it))+
    geom_point(aes(y=Meand, col = Bounds),size = 1, alpha = .1, )+
    geom_ribbon(aes(ymin = Lowerd, ymax = Upperd),fill = "black",col = "black", alpha = .2)+
    #geom_vline(aes(yintercept = 0))+
    scale_x_log10()+
    #scale_y_log10()+
    scale_color_manual(values = c("forestgreen","red"), labels = c("Remains Within","Exceeds"))+
    theme_bw()+
    ylab("Residual Norm Squared")+
    xlab("Iteration")+
    ggtitle(name)+ 
    theme(legend.position="none")+
    theme(text = element_text(size = 18))
  return(g)
}
# Takes output from get_all_obs and adds stopping criteria information

stopping_criteria <- function(results, nu, delt1, delt2, xi1, xi2, eta = 1){
  Sampling_meth <- results$Sampling_Method[1]
  width <- results$Width
  Sample_size <- results$Sample_size
  #Calculate Appropiate constants
  if(Sampling_meth == "Achlio" | Sampling_meth == "gauss"){
    c <- 2.12258 * (1+log(width)) / (Sample_size*width)
    omega <- .1
  }else{
    c <-  32 *  (1 + log(width))/ (Sample_size*width)
    omega <- 1/16
  }

  
  #Calculate the thresholds
  t1 <- 1/c * eta * min((delt1 - 1)^2 * nu^2, (delt1 - 1) * nu / omega)/(2 * log(1/xi1))
  t2 <- 1/c * eta *  min((1 - delt2)^2 * nu^2, (1 - delt2) * nu / omega)/(2 * log(1/xi2))
  thres <- ifelse(t1 > t2, t2, t1)
  # Apply the thresholds to determine if satisfied
  return( results %>%
            ungroup()%>%
            mutate(thres = thres) %>%
    mutate(stopped = ifelse(Mean_Obs_Res2 <= thres, T, F),
           l_nu = ifelse(Mean_True_Res <= nu, T, F),
           nu = nu, 
           delt1 = delt1, 
           delt2 = delt2, 
           xi1 = xi1,
           xi2 = xi2)%>%
      mutate(error = case_when(stopped == T & l_nu == T ~ "Correct",
                               stopped == F & l_nu == T ~ "Type 1",
                               stopped == T & l_nu == F ~ "Type 2"))
  )
    
}

#Function designed to take raw stopping criteria information and 
#create a graph
graph_stop <- function(stop){
  g <- stop %>%
    ggplot(aes(x = error, fill = error))+
    geom_bar()
  
  return(g)
}

stop %>%
  filter(Mean_Obs_Res2 <= thres,l_nu == T) %>%
  group_by(reps) %>%
  summarise(cond_t1 = mean(Mean_Obs_Res > delt1 * nu))%>%#mean(Mean_Obs_Res > delt1 * nu))
summarise(mean(cond_t1))

stop %>%
  filter(Mean_Obs_Res2 <= thres, l_nu == F) %>%
  group_by(reps) %>%
  summarise(cond_t2 = mean(Mean_Obs_Res < delt2 * nu))%>%
  summarise(mean(cond_t2))#mean(Mean_Obs_Res < delt2 * nu))

stop %>%
  filter(Mean_Obs_Res2 <= thres, l_nu == F) %>%
  group_by(reps) %>%
  summarise(cond_t2 = mean(Mean_Obs_Res > delt2 * nu & Mean_Obs_Res < delt1 * nu))%>%
  summarise(mean(cond_t2))

stop %>%
  summarise(satisfied = mean(l_nu == T),outside = mean(l_nu == F))

stop %>%
  filter(Mean_Obs_Res2 <= thres) %>%
  mutate(Error= case_when(l_nu == T & Mean_Obs_Res > delt1 * nu ~ "Type 1",
                          l_nu == F & Mean_Obs_Res < delt2 * nu ~ "Type 2",
                          TRUE ~ "No Error")) %>%
  ggplot(aes(x = reps, y= Mean_Obs_Res))+
  geom_point(aes(col = Error))+
  geom_line(aes(y = nu))+
  geom_line(aes(y = delt1 * nu, col = l_nu))+
  geom_line(aes(y = delt2 * nu))+
  facet_wrap(~l_nu)
