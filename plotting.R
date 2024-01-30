DIRECT <- "~/Documents/Research/metrics-linearsystems-v2/code"
setwd(DIRECT)


library(xtable)
library(tidyverse)
library(scales)
source("analyze.R")
source("plotter.R")
source("analyze_rep.R")
source("Repeat_analysis.R")
results <- run_output("./Old_computer/Uresults25")

ests <- map(seq(1,500,5),~estimators(results,.x))
tvs <- bind_rows(ests)

tvs_sub<-tvs%>% filter(Iteration>10, Width > 1,Sampling_Method == "gauss",Matrix_Type %in% c("golub","phillips","baart"))
lmer(log10(tvs_sub$Error)~tvs_sub$Matrix_Type+log10(tvs_sub$True_Res)+(1|tvs_sub$Width)) %>% summary
lmer(log10(tvs_sub$Error)~tvs_sub$Matrix_Type+log10(tvs_sub$Mov_True_Res)+(1|tvs_sub$Width)) %>% summary
weighted_model <- lm(log10(tvs_sub$Error)~0+ tvs_sub$Matrix_Type* log10(tvs_sub$Mov_True_Res),weights = 1/tvs_sub$Width)
unweighted_model <- lm(log10(tvs_sub$Error)~0+tvs_sub$Matrix_Type*log10(tvs_sub$True_Res))

tvs_sub%>%add_predictions(weighted_model) %>%add_residuals(weighted_model) %>%filter(Iteration >5, Matrix_Type == "baart")%>%ggplot(aes(x = Mov_True_Res, y =log10(Error),color=Matrix_Type))+geom_point()+geom_line(aes(x = Mov_True_Res, y = pred),color="blue")+scale_x_log10()+theme(legend.position = "none")
tvs_sub%>%add_predictions(unweighted_model) %>%add_residuals(unweighted_model) %>%filter(Iteration >5, Matrix_Type == "baart")%>%ggplot(aes(x = True_Res, y =log10(Error),color=Matrix_Type))+geom_point()+geom_line(aes(x = True_Res, y = pred),color="blue")+scale_x_log10()+theme(legend.position = "none")
summary(unweighted_model)
summary(weighted_model)
tvs_sub%>%add_predictions(weighted_model) %>%add_residuals(weighted_model)%>%
  rename(resid1 = resid)%>%add_residuals(unweighted_model)%>%filter(Matrix_Type == "baart")%>%
  ggplot(aes(x = pred,y =resid))+
  geom_point()+
  geom_point(aes(y=resid1/Width),color = "blue")+
  theme(legend.position = "none")

tvs_sub %>% filter(Matrix_Type == "baart")%>%ggplot(aes(x = Iteration, y = Error))+ geom_point()
est_100 <- estimators(results, 100)


estimators_100 <- confidence_add(est_100, .95)

estimators_100 <- stopping_add(estimators_100, 100, 1.1,.9,.01,.01)

estimators_100 <- estimators_100 %>%
  drop_na()

est_100 %>% filter(Matrix_Type == "phillips", Sampling_Method == "Uniform") %>% filter(Iteration > 3) %>%
  select(Iteration,Error,Mov_True_Res,True_Res, Obs_Res, Mov_Obs_Res) %>% write_csv("../article/csvs/phillips.csv")
#Look at big system solve example
results_big <- run_output("./big.csv")

results_big <- estimators(results_big, 100)

results_big <- confidence_add(results_big, .95,F)

results_big <- stopping_add(results_big, 1, 1.1,.9,.01,.01)

#Big experiment Results
width <- 1000
read_csv("./big_sre.csv") %>%
   filter(value > 1) %>%
   group_by(type) %>%
   mutate(Iteration = row_number()) %>%
   ungroup() %>%
    pivot_wider(id_cols = "Iteration",names_from = "type", values_from = "value")%>%
   mutate(rho = as.numeric(rho), iota = rho^2) %>% 
  mutate(Sample_Size = 50, 
         Nrows = 132651, 
         Width = width,
         rho = as.numeric(rho) * (Nrows/Sample_Size),  
         iota = as.numeric(iota) * (Nrows/Sample_Size), 
         Obs_Res = rho ) %>%
   estimatorsBig(width) %>%
  mutate(confidence(rho,iota,width,"uniform",50,.95,16, nrow = Nrows, sigma2 = 0)) %>%
  mutate(rho = (rho^(1)), Upper = Upper^1, Lower = sign(Lower) * abs(Lower)^(1)) %>%
  rename(it = Iteration, Mov_Obs_Res = rho) %>% 
  mutate(ratio = Obs_Res/Mov_Obs_Res) %>% 
  filter(it %% 25 == 0) %>%
  #drop_na() %>% ggplot(aes(x = it, y = Obs_Res)) + geom_line()+geom_line(aes(x = it, y = Mov_Obs_Res), color = "red")+geom_line(aes(x = it, y = Upper), color = "green")+scale_y_log10()
  write_csv("../article/csvs/big.csv")
  
  # ggplot(aes(x= Iteration))+
  # geom_line(aes(y =Upper))+
  # geom_line(aes(y = Lower))+
  # geom_line(aes(y  = rho), color = "red")
 
 res_big<-read_big("big3.csv",20,100^3,0.02,0.107807,nu=200,eta = 1)
 res_big<-read_big("big3.csv",20,100^3,0.02,.07,500,nu=200,eta = 1)
 res_big100<-read_big("big31.csv",20,100^3,0.02,.172,300,nu=400,eta = 1)
 res_big100 <- res_big100 %>% mutate(conf_calc(Mov_Obs_Res,Mov_Obs_Res2,vari,Width,alpha))
 res_big100 %>% stopping_graph_big()
 res_big100 %>% summarise(mean(Mov_True_Res> Upper | Mov_True_Res < Lower))
 res_big100 %>%
   filter( it < 3000) %>%
   ggplot(aes(x = it))+
   geom_line(aes(y=Lower))+
   geom_line(aes(y=Upper))+
   #geom_line(data = res_big100,aes(y=Lower),color = "blue",alpha = .5)+
   #geom_line(data = res_big100,aes(y=Upper),color = "blue",alpha = .5)+
   geom_line(aes(y=Mov_Obs_Res), col = "red")+
   geom_line(aes(y=Mov_True_Res), col = "green")+
   scale_y_log10()
    #geom_hline(aes(yintercept = 100))
 res_big100 %>%
   filter(it %% 10 == 1, it < 3000) %>%
   write_csv("../article/csvs/big.csv")
 
 compute_r_from_rho <- function(widths,rhos){
   ln = length(rhos)
   res <- numeric(ln)
   #mw <- max(widths)
   for(i in 1:ln){
     if(widths[i] == 1){
       #print(i)
       res[i] = rhos[i]
      }else if(widths[i-1] != widths[i]){
        #print(i)
        res[i] = rhos[i] * widths[i] - rhos[i-1]*widths[i-1]
      }else{
        #print(i)
        res[i] = rhos[i] * widths[i] - rhos[i-1]*widths[i-1] + res[i - widths[i] - 1] 
     }
   }
   return(res)
 }
 
 get_var_rel_err<-function(data,nvals){
   res_o <- compute_r_from_rho(data$width,data$rho)
   res_e <- compute_r_from_rho(data$width,data$trho)
   var(abs(res_o[1:nvals] - res_e[1:nvals])/ res_e[1:nvals])
 }
 
 new_mov_av<-function(data,wi){
   data %>% mutate(Iteration = it) %>%
     mutate(True_Res = compute_r_from_rho(width,trho))%>%
      mutate(Obs_Res = compute_r_from_rho(width,rho)) %>%
                  #mutate(iota = compute_r_from_rho(width,iota)) %>%
      estimators(wi)%>%
     return()
 }
 
 #Consider Everything
 all_rho <- read_all_rep("./Uresults_rep/", nuv = 1, alph = .05, xi1v = .01, xi2v  =.01, eta = 5)
 all_rho <- read_all_rep("./urep_results_rho1024/", nuv = 1, alph=.05, xi1v = .01, xi2v = .01, eta = 100000)
 #all_rho <- read_all_rep("./Uresults_rep/", nuv = 1, sigma2 = 0)
 all_rho %>%
   mutate(rel_err = min(rel_err, 512/25)) %>%
   group_by(Matrix_Type) %>%
   mutate(serr = max(rel_err)^2/4) %>%# filter(reps==1, Matrix_Type == "baart") %>% view()
   mutate(Upper_est = Upper_est * sqrt(2^2) / (Rows/(2*Sample_size)), 
                    Lower_est = Lower_est * sqrt(2^2) / (Rows/(2*Sample_size))) %>% 
   group_by(Matrix_Type) %>%
   #ungroup() %>%
   summarise(fails = mean(rho > Upper_est | rho < Lower_est)) %>% 
   arrange(desc(fails)) %>% view()
   filter(Matrix_Type == "lotkin") %>% ggplot(aes(x = Iteration, y = Upper_est))+geom_line(color = "black")+geom_point(aes(y=Mov_Obs_Res),alpha = .1, color = "green")+ geom_point(aes(y=Mov_True_Res),alpha = .1, color = "blue") + scale_y_log10()+geom_line(aes(y=orho), color = "red")+geom_line(aes(y=rho), color = "violet")
 # Examine the stopping criteria with no eta
 all_rho %>%
   #filter(Mov_Obs_Res2 <= threshold) %>%
   stopping_graph_rep() %>% view()
 all_rho %>%
   #filter(Mov_Obs_Res2 <= threshold) %>%
   stopping_graph_rep_nc() #%>% view()
 # Get overall failure rate with no eta
 all_rho %>%
   #group_by(Matrix_Type) %>%
   ungroup() %>%
   summarise(fails = mean(Mov_True_Res > Upper_est | Mov_True_Res < Lower_est)) %>% arrange(desc(fails)) %>% view()
 #Code to make graphs
 group_by(Iteration) %>%
   mutate(Upper_est = Upper_est/sqrt(21), outside = Mov_True_Res > Upper_est ) %>%
   ggplot(aes(x = Iteration, y = Upper_est)) +
   geom_point(aes(y = Mov_True_Res, color = outside), alpha = .1) +
   geom_line(aes(y = rho), color = "blue") +
   geom_line(aes(y=orho), color = "red")+
   geom_line() +
   geom_line(aes(y = max(Lower_est,1e-8))) +
   scale_y_log10()
 
 #Get Information for Wilkinson Graph no eta
 all_rho %>%
   filter(Matrix_Type == "rosser", Sampling_Method == "Uniform") %>%
   filter(Mov_True_Res <= Upper_est & Mov_True_Res >= Lower_est, row_number() %in% sample(1:nrow(.), 2709, replace = F)) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(min_up =  Upper_est, max_low = Lower_est, rho = Mov_True_Res) %>%
   #summarise(rho = max(rho),orho = max(orho), min_up = -1*max(-Upper_est)/sqrt(17), max_low = -1 *  min(-Lower_est)/sqrt(17)) %>%
   #filter(Mov_True_Res <= Upper_est | Mov_True_Res >= Lower_est) %>%
   #summarise(maxrho = max(Mov_True_Res), minrho = min(Mov_True_Res), 
  #           min_up = max(Upper_est), max_low = max(min(Lower_est),1e-7)) %>%
   write_csv("../article/csvs/wilkrepinside1.csv")
 
 

 all_rho %>%
   filter(Matrix_Type == "rosser", Sampling_Method == "Uniform") %>%
   filter(Mov_True_Res >= Upper_est | Mov_True_Res <= Lower_est) %>%
   group_by(Iteration) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(rho = Mov_True_Res) %>%
   write_csv("../article/csvs/wilkrepoutside1.csv")  
 set.seed(123)
 all_rho %>%
   filter(Matrix_Type == "phillips", Sampling_Method == "Uniform") %>%
   group_by(Iteration) %>%
   #summarise(rho = max(rho),orho = max(orho), min_up = -1*max(-Upper_est)/sqrt(17), max_low = -1 *  min(-Lower_est)/sqrt(17)) %>%
   filter(Mov_True_Res <= Upper_est & Mov_True_Res >= Lower_est, row_number() %in% sample(1:nrow(.), 2684, replace = F)) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(min_up =  Upper_est, max_low = Lower_est, rho = Mov_True_Res) %>%
   #summarise(maxrho = max(Mov_True_Res), minrho = min(Mov_True_Res), 
   #         min_up = max(Upper_est), max_low = max(min(Lower_est),1e-23)) %>%
   write_csv("../article/csvs/rohrepinside1.csv")

 
 all_rho %>%
   filter(Matrix_Type == "phillips", Sampling_Method == "Uniform") %>%
   filter(Mov_True_Res >= Upper_est | Mov_True_Res <= Lower_est)  %>%
   filter(row_number() %in% sample(1:nrow(.), 23, replace = F)) %>%
   group_by(Iteration) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(rho = Mov_True_Res)%>%
   write_csv("../article/csvs/rohrepoutside1.csv")  
 
 
 
 #Do Experiment to find the appropiate eta Values
 
 gauss <- all_rho %>% filter(Sampling_Method == "gauss")
 FJLT <- all_rho %>% filter(Sampling_Method == "FJLT")
 achlio <- all_rho %>% filter(Sampling_Method == "Achlio")
 uniform <- all_rho %>%ungroup() %>% filter(Sampling_Method == "Uniform")
 alph = .95
 # Begin by testing the correct eta for the rho values
 find_eta <- function(data,start,en,alph){
   test_va = 1 -alph
   fr <- 0
    for(eta1 in start:en){
      old <- c(fr,eta1 - 1)
      fr <- data %>%
      #confidence_add_eo(alph,eta) %>%
      ungroup() %>%
      mutate(Upper_est = orho + ((Upper_est - orho)/ sqrt(eta1)), Lower_est = orho + (Lower_est - orho)/sqrt(eta1)) %>%
      summarise(me = mean(Mov_True_Res > Upper_est | Mov_True_Res < Lower_est)) %>%
      pull(me)
      if(fr > test_va){
        return(matrix(append(old,c(fr,eta1)),ncol = 2, byrow = T))
      }
    }
   return(matrix(append(old,c(fr,eta1)),ncol = 2, byrow = T))
 }
 
 gauss_et<- find_eta(gauss,5,20,alph)
 gauss_et
 ach_et<- find_eta(achlio,5,20,alph)
 ach_et
 FJLT_et <- find_eta(FJLT,185,190,alph)
 FJLT_et
 Uniform_et <- find_eta(uniform,30,50,alph)
 Uniform_et
 # Use Gauss = Achlio = 5, FJLT = 188
 #Do Stopping experiment again with new eta values
 # overall_eta <- gauss %>% 
 #   mutate(obs_res = Mov_Obs_Res ) %>%
 #   confidence_add_eo(alph,26) %>%
 #   mutate(Iteration = Iteration - 14) %>%
 #   stopping_add_rep(100, 1.1,.9,.01,.01,26)%>%
 #   bind_rows(achlio %>% 
 #               mutate(obs_res = Mov_Obs_Res ) %>%
 #               confidence_add_eo(alph,26) %>%
 #               mutate(Iteration = Iteration - 14) %>%
 #               stopping_add_rep(1, 1.1,.9,.01,.01,26)) %>%
 #   bind_rows(FJLT %>% 
 #               mutate(obs_res = Mov_Obs_Res ) %>%
 #               confidence_add_eo(alph,188) %>%
 #               mutate(Iteration = Iteration - 14) %>%
 #               stopping_add_rep(1, 1.1,.9,.01,.01,188))
 overall_eta <- uniform %>% 
   mutate(obs_res = Mov_Obs_Res ) %>%
   #confidence_add_eo(alph,21) %>%
   mutate(Upper_est = orho + (Upper_est - orho)/sqrt(41), Lower_est = orho + (Lower_est - orho)/sqrt(41)) %>%
   mutate(Iteration = Iteration - 30) %>%
   #stopping_add_rep(1e-10, 1.1,.9,.01,.01) %>%
   mutate(threshold = threshold * (60))

 # Examine the stopping criteria with no eta
 view(overall_eta %>%
   stopping_graph_rep())
 overall_eta %>%
   stopping_graph_rep_nc()
 # Get overall failure rate with no eta
 overall_eta %>%
   group_by(Matrix_Type) %>%
   #ungroup() %>%
   summarise(mean(Mov_True_Res > Upper_est | Mov_True_Res < Lower_est)) %>% view()
 #Get Information for Wilkinson Graph no eta
 overall_eta %>%
   filter(Matrix_Type == "rosser", Sampling_Method == "Uniform") %>%
   group_by(Iteration) %>%
   filter(Mov_True_Res <= Upper_est & Mov_True_Res >= Lower_est, row_number() %in% sample(1:nrow(.), 2586, replace = F)) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(min_up =  Upper_est, max_low = Lower_est, rho = Mov_True_Res) %>%
   #filter(Mov_True_Res <= Upper_est | Mov_True_Res >= Lower_est) %>%
   #summarise(maxrho = max(Mov_True_Res), minrho = min(Mov_True_Res), 
    #         min_up = max(Upper_est), max_low = max(min(Lower_est),1e-7)) %>% 
   write_csv("../article/csvs/wilkrepinside1_eta.csv")
 
 overall_eta %>%
   filter(Matrix_Type == "rosser", Sampling_Method == "Uniform") %>%
   filter(Mov_True_Res >= Upper_est | Mov_True_Res <= Lower_est ) %>%
   filter(row_number() %in% sample(1:nrow(.), 120, replace = F)) %>%
   group_by(Iteration) %>%
   select(Iteration, Mov_Obs_Res, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(rho = Mov_True_Res) %>%
   group_by(Iteration) %>%
   write_csv("../article/csvs/wilkrepoutside1_eta.csv")  
 
 overall_eta %>%
   filter(Matrix_Type == "phillips", Sampling_Method == "Uniform") %>%
   group_by(Iteration) %>%
   filter(Mov_True_Res <= Upper_est & Mov_True_Res >= Lower_est, row_number() %in% sample(1:nrow(.), 2467, replace = F)) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(min_up =  Upper_est, max_low = Lower_est, rho = Mov_True_Res) %>%
   #filter(Mov_True_Res <= Upper_est | Mov_True_Res >= Lower_est) %>%
   #summarise(maxrho = max(Mov_True_Res), minrho = min(Mov_True_Res), 
  #           min_up = max(Upper_est), max_low = max(min(Lower_est),1e-23)) %>%
   write_csv("../article/csvs/rohrepinside1_eta.csv")
 
 overall_eta %>%
   filter(Matrix_Type == "phillips", Sampling_Method == "Uniform") %>%
   filter(Mov_True_Res >= Upper_est | Mov_True_Res <= Lower_est) %>%
   filter( row_number() %in% sample(1:nrow(.), 241, replace = F)) %>%
   group_by(Iteration) %>%
   select(Iteration, Mov_True_Res, Upper_est, Lower_est) %>%
   rename(rho = Mov_True_Res) %>%
   write_csv("../article/csvs/rohrepoutside1_eta.csv")  
 
 overall_eta %>% 
   filter(Sampling_Method == "Uniform", Matrix_Type == "foxgood") %>% 
   mutate(abs_diff = abs(rho - Mov_Obs_Res)) %>% group_by(Iteration, rho) %>% 
   summarise(obs_res = max(abs_diff),mean_diff = median(abs_diff), min_diff = min(abs_diff)) %>%
   write_csv("../article/csvs/iotarho.csv")
 
 overall_eta %>%
   ungroup() %>%
   mutate(rel_err = abs(obs_res - rho)/rho) %>%
   group_by(Iteration) %>%
   #rename(Iteration = it)%>%
   summarise("5th" = min(rel_err),"50th" = quantile(rel_err,.5),"95th" = max(rel_err) ) %>%
   write_csv("../article/csvs/convergedrho.csv")
 
 
 
 
 
 # mutate(Upper_est = Upper_est, outside = Mov_True_Res > Upper_est ) %>%
 #   ggplot(aes(x = Iteration )) +
 #   geom_point(aes(y = Mov_True_Res, color = outside), alpha = .1) +
 #   geom_line(aes(y = rho), color = "blue") +
 #   geom_line(aes(y= orho), color = "red")+
 #   geom_line(aes(y = Upper_est)) +
 #   geom_line(aes(y = max(Lower_est,1e-8))) +
 #   scale_y_log10()
 
 
 
 
 
 
  
#Consider results for repeated samples
rep_wilk <- read_it_rep("./results_rep/","wilkinson_256_256_gauss_25_15.csv")
rep_golub <- read_it_rep("./results_rep/","golub_256_256_gauss_25_15.csv")
rep_rohess <- read_it_rep("./results_rep/","rohess_256_256_gauss_25_15.csv")
rep1_wilk <- confidence_rep(rep_wilk, 15, .95,eta=1)
rep1_golub <- confidence_rep(rep_golub, 15, .95)
rep1_rohess <- confidence_rep(rep_rohess, 15, .95)
sum_wilk <- get_all_obs_us(rep_wilk, rep1_wilk, 15)

#Consider repeated samples with how often uncertainty sets fail to cover rho
rep_wilk <- read_it_rep("./results_rep/","wilkinson_256_256_gauss_25_15.csv")
ci_wilk <- confidence_rep_eo(rep_wilk,15,.95)
ci_wilk %>%
  summarise(mean(rho > Upper_est | rho< Lower_est))
#Save the values inside
ci_wilk %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho <= min_up | rho >= max_low) %>%
  write_csv("../article/csvs/wilkrepinside1.csv")
#save the values outside
ci_wilk %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho > min_up | rho< max_low) %>%
  write_csv("../article/csvs/wilkrepoutside1.csv") 


rep_golub <- read_it_rep("./results_rep/","golub_256_256_gauss_25_15.csv")
ci_golub <- confidence_rep_eo(rep_golub,15,.95)
ci_golub %>%
  summarise(mean(rho > Upper_est | rho< Lower_est))
#Save the values inside
ci_golub %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho <= min_up | rho >= max_low) %>%
  write_csv("../article/csvs/golubrepinside1.csv")
#save the values outside
ci_golub %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho > min_up | rho< max_low) %>%
  write_csv("../article/csvs/golubrepoutside1.csv") 

rep_rohess <- read_it_rep("./results_rep/","rohess_256_256_gauss_25_15.csv")
ci_rohess <- confidence_rep_eo(rep_rohess,15,.95)
ci_rohess %>%
  summarise(mean(rho > Upper_est | rho< Lower_est))
#Save the values inside
ci_rohess %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho <= min_up | rho >= max_low) %>%
  write_csv("../article/csvs/rohrepinside1.csv")
#save the values outside
ci_rohess %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho > min_up | rho< max_low) %>%
  write_csv("../article/csvs/rohrepoutside1.csv")















#Do the Same Analysis with etas = 13



rep_wilk <- read_it_rep("./results_rep/","wilkinson_256_256_gauss_25_15.csv")
ci_wilk <- confidence_rep_eo(rep_wilk,15,.95, eta =  5)
ci_wilk %>%
  summarise(mean(rho > Upper_est | rho< Lower_est))
#Save the values inside
ci_wilk %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho <= min_up | rho >= max_low) %>%
  write_csv("../article/csvs/wilkrepinside1_eta.csv")
#save the values outside
ci_wilk %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho > min_up | rho< max_low) %>%
  write_csv("../article/csvs/wilkrepoutside1_eta.csv") 


rep_golub <- read_it_rep("./results_rep/","golub_256_256_gauss_25_15.csv")
ci_golub <- confidence_rep_eo(rep_golub,15,.95, eta =  5)
ci_golub %>%
  summarise(mean(rho > Upper_est | rho< Lower_est))
#Save the values inside
ci_golub %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho <= min_up | rho >= max_low) %>%
  write_csv("../article/csvs/golubrepinside1_eta.csv")
#save the values outside
ci_golub %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho > min_up | rho< max_low) %>%
  write_csv("../article/csvs/golubrepoutside1_eta.csv") 

rep_rohess <- read_it_rep("./results_rep/","rohess_256_256_gauss_25_15.csv")
ci_rohess <- confidence_rep_eo(rep_rohess,15,.95, eta =  5)
ci_rohess %>%
  summarise(mean(rho > Upper_est | rho< Lower_est))
#Save the values inside
ci_rohess %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho <= min_up | rho >= max_low) %>%
  write_csv("../article/csvs/rohrepinside1_eta.csv")
#save the values outside
ci_rohess %>% 
  group_by(Iteration) %>%
  summarise(rho = max(rho), min_up = -1*max(-Upper_est), max_low = -1 *  min(-Lower_est)) %>%
  filter(rho > min_up | rho< max_low) %>%
  write_csv("../article/csvs/rohrepoutside1_eta.csv")


#Get repeat Rho for the first graph in experiments

rhos <- run_rhos("./rep_results_rho")

rhos %>% filter(Sampling_Method == "gauss",Matrix_Type == "rohess") %>%
  write_csv("../article/csvs/iotarho.csv")

#ungroup()%>%
#group_by(it)%>%
#summarise(rho = mean(rho), obs_res = max(Mov_Obs_Res))%>%
#ggplot()+
#geom_line(aes(x = Iteration, y = rho))+#abs(rho - obs_res)))+
#scale_y_log10()+
#geom_line(aes(x = Iteration, y=obs_res), color ="red", linetype = "dashed")

rhos %>%
  ungroup() %>%
  mutate(rel_err = abs(obs_res - rho)/rho) %>%
  group_by(Iteration) %>%
  #rename(Iteration = it)%>%
  summarise("5th" = min(rel_err),"50th" = quantile(rel_err,.5),"95th" = max(rel_err) ) %>%
  write_csv("../article/csvs/convergedrho.csv")

rhos %>% 
  filter(Matrix_Type == "golub") %>%
  select(Iteration,Matrix_Type,Sample_size,Mov_Obs_Res,rho)%>%
  ggplot()+
  geom_line(aes(x = Iteration, y = Mov_Obs_Res))+
  geom_line(aes(x = Iteration, y=rho), color ="red", linetype = "dashed")

#look at Stopping information in experiments
st<-stopping_add_rep(rhos, 1, 1.1,.9,.01,.01,1)
stopping_graph_rep(st) %>%
  write_csv("../article/csvs/stopping_graph.csv")
st_eta<-stopping_add_rep(rhos, 1, 1.1,.9,.01,.01,13)
stopping_graph_rep(st_eta)

#Experiment to find optimal eta for each sampling method








#Make csv for points inside curve
sum_wilk %>%
  filter(Meand <= Upperd & Meand >= Lowerd) %>%
  group_by(it) %>%
  summarise(max = max(Meand), min = min(Meand), up = max(Upperd), low = min(Lowerd))%>%
  write_csv("../article/csvs/wilkrepinside.csv")

sum_wilk %>%
  filter(Meand > Upperd | Meand < Lowerd) %>%
  select(it,Meand) %>%
  write_csv("../article/csvs/wilkrepoutside.csv")
  
sum_golub <- get_all_obs_us(rep_golub, rep1_golub, 15)
sum_golub %>%
  filter(Meand <= Upperd & Meand >= Lowerd) %>%
  group_by(it) %>%
  summarise(max = max(Meand), min = min(Meand), up = max(Upperd), low = min(Lowerd))%>%
  write_csv("../article/csvs/golubrepinside.csv")

sum_golub %>%
  filter(Meand > Upperd | Meand < Lowerd) %>%
  select(it,Meand) %>%
  write_csv("../article/csvs/golubrepoutside.csv")


sum_rohess <- get_all_obs_us(rep_rohess, rep1_rohess, 15)

sum_rohess %>%
  filter(Meand <= Upperd & Meand >= Lowerd) %>%
  group_by(it) %>%
  summarise(max = max(Meand), min = min(Meand), up = max(Upperd), low = min(Lowerd))%>%
  write_csv("../article/csvs/rohrepinside.csv")

sum_rohess %>%
  filter(Meand > Upperd | Meand < Lowerd) %>%
  select(it,Meand) %>%
  write_csv("../article/csvs/rohrepoutside.csv")

sum_rohess %>% ungroup%>%summarise(mean(Bounds == "yes"))
sum_golub %>% ungroup%>%summarise(mean(Bounds == "yes"))
sum_wilk %>% ungroup%>%summarise(mean(Bounds == "yes"))

# Repeat of stopping


#Make Plot for rhos
pr<-plot_rho(estimators_100,rel=T)

plot_rho(estimators_100)
pr %>% 
  #filter(Iteration %% 20 == 0) %>% 
  write_csv("../article/csvs/convergedrho.csv")


ggsave("./plots/rho.png",width = 6, height =3)

plot_rho(estimators_100, rel = T)
ggsave("../plots/rho_rel.png",width = 6, height =3)
#Make Plot for iota

pi <- plot_iota(estimators_100, rel = T)
pi %>% 
  #filter(Iteration %% 20 == 0) %>% 
  write_csv("../article/csvs/convergediota.csv")

ggsave("./plots/iota.png",width = 6, height =3)

plot_iota(estimators_100, rel = T)
ggsave("./plots/rho_rel.png",width = 6, height =3)
#the rhoviota
estimators_100 %>%
  filter(Matrix_Type == "rohess", Sampling_Method == "gauss") %>%
  select(Iteration,Matrix_Type, Sampling_Method, Sample_size, Mov_Obs_Res,Mov_Obs_Res2)%>%
  rename( "iota"=Mov_Obs_Res2 ,  "rho"=Mov_Obs_Res ) %>%
  write_csv("../article/csvs/iotarho.csv")
#Make plot for true and estimated Confidence Interval
Matrix_types <-  c("golub","wilkinson")
plot_int(estimators_100,Matrix_types,"gauss",100, type = "both")
plot_int(estimators_100,Matrix_types,"gauss",100, type = "est")
ggsave("../submission/plots/ci_plots.png",width = 6, height =3)
#Make full results table
tab_ci <- confidence_table(estimators_100, width = 100)
tab_ci[,7:11] <- lapply(tab_ci[,7:11], function(x)
  paste0("\\cellcolor{", colortab(x), "}", ifelse(is.na(x),0,x)))
#Add color gradient to the table
tab_ci%>%
  filter(`Sampling Method` != "SRHT")%>%
  select(-Sample_Size) %>%
  rename(`$>10$` = `>10`,`$>1e-5$` = `>1e-5`,`$>1e-10$` = `>1e-10`,`$>0$` = `>0`)%>%
  xtable(type = "latex",digits=c(0,0,0,0,0,0,4,4,4,4,4))%>%
  print(tabular.environment="longtable",file = "../metrics-leastsquares/GenCoD/Article/plots/ci_table_noeta.tex",
        sanitize.text.function = identity, include.rownames=FALSE)

#Make results table for specific Matrix type
MAT = "gravity"

tab_ci %>%
  filter(`Sampling Method` != "SRHT")%>%
  filter(`Matrix Type` == MAT)%>%
  select(-eta, -Sample_Size, -Width, -Rows, -`Matrix Type`)%>%
  rename(`$>10$` = `>10`,`$>1e-5$` = `>1e-5`,`$>1e-10$` = `>1e-10`,`$>0$` = `>0`)%>%
  xtable(type = "latex",digits=c(0,0,4,4,4,4,4))%>%
  print(tabular.environment="longtable",file = "../submission/plots/exp_table.tex",
        sanitize.text.function = identity, include.rownames=FALSE)

#Make bar plot for stopping point errors 
stopping_graph(estimators_100) %>%
  write_csv("../article/csvs/stopping_graph.csv")

ggsave("./plots/stop_true.png",width = 4, height = 2 )
plot_stop_no_delt(estimators_100)
ggsave("../plots/no_delt_stop.png",width = 4, height = 2 )
#Make plots from large matrix tests
#plot_int(results_big,"ls","gauss",100, type = "est")

plot_int(results_big,"rohess","gauss",100, type = "est")
ggsave("../submission/plots/rohess_big.png",width = 4, height = 2 )
#Make Plots for repeated Sample Test
make_graphs(sum_wilk)
ggsave("../submission/plots/rep_wilk.png",width = 6, height =3 )
make_graphs(sum_golub)
ggsave("../submission/plots/rep_golub.png",width = 6, height =3 )
make_graphs(sum_rohess)
ggsave("../submission/plots/rep_rohess.png",width = 6, height = 3 )

#Do an eta analysis
eta_res <- tibble( "Sampling Method" = NULL, "Width" = NULL,        
 "Matrix Type" = NULL,    "Rows"  = NULL, "Sample_Size" = NULL,
 "eta" = NULL, ">10" = NULL, ">1e-5" = NULL, ">1e-10" = NULL,
 ">0" = NULL,  "Total" = NULL)
results25 <- run_output("./results25")
results50 <- run_output("./results50")
results200 <- run_output("./results200")
for(i in c(15,30,100)){
  etas <- seq(1,1000)
  
  est_25 <- estimators(results25, i)
  eta_list_25 <- map(etas,~ eta_test(est_25,.95,i,.x))
  eta_res_25 <- bind_rows(eta_list_25)

   
   est_100_100 <- estimators(results, i)
   eta_list_100 <- map(etas,~ eta_test(est_100_100,.95,i,.x))
   eta_res_100 <- bind_rows(eta_list_100)
   
   
   est_100_50 <- estimators(results50, i)
   eta_list_50 <- map(etas,~ eta_test(est_100_50,.95,i,.x))
   eta_res_50 <- bind_rows(eta_list_50)

   
   est_100_200 <- estimators(results200, i)
   eta_list_200 <- map(etas,~ eta_test(est_100_200,.95,i,.x))
   eta_res_200 <- bind_rows(eta_list_200)
   
  eta_res <- rbind(eta_res,eta_res_25,eta_res_50,eta_res_100, eta_res_200)
}



eta_res%>%
  mutate(Sample_Size = round(Sample_Size / Rows,2))%>%
  drop_na()%>%
  filter(`>0`==0 & `>1e-10` < .7)%>%
  group_by(Width,`Sampling Method`,`Sample_Size`,eta)%>%
  summarise(Total = mean(Total))%>%
  mutate(eta = eta,Width = as.factor(Width), `Sample_Size` = as.factor(`Sample_Size`))%>%
  filter(Total <= .05)%>%
  group_by(`Sampling Method`, `Sample_Size`,Width)%>%
  mutate(max = max(Total))%>%
  filter(Total == max)%>%
  mutate(min = min(eta))%>%
  filter(eta == min)%>%
  mutate(eta = eta)%>%
  filter(`Sampling Method` != "SRHT" )%>%
  pivot_wider(id_cols = c(Width,`Sample_Size`),names_from = `Sampling Method`,values_from = eta)%>%
  rename("Percent Sampled" = `Sample_Size`)%>%
  xtable(type = "latex", digits=c(0,0,3,0,0,0))%>%
  print(tabular.environment="longtable",file = "../metrics-leastsquares/GenCoD/Article/plots/eta_tab.tex",
        sanitize.text.function = identity,include.rownames=FALSE)

eta_res%>%
  filter(`Sampling Method` != "SRHT")%>%
  mutate(Total = `>10`+`>1e-5`+`>1e-10`)%>%
  group_by(Width,Sample_Size,eta,`Sampling Method`)%>%
  summarise(Total = mean(Total))%>%
  drop_na()  %>%
  mutate(Width = as.factor(Width), Sample_Size = as.factor(Sample_Size))%>%
  ggplot(aes(x = eta, y = Total, col = Sample_Size))+
  geom_line()+
  geom_hline(yintercept = .05, color = "red", linetype = "dotted")+
  facet_wrap(`Sampling Method`~Width,nrow = 1)+
  xlab(expression(eta))+
  ylab("Proportion exceeding CI")+
  #scale_x_log10()+
  scale_color_manual(values = c("darkorange","purple","navy"),name = "Sample Dimension")+
  theme_bw()+
  ggtitle("Average proportion of iterations where true moving average exceeds a\n95% credible interval")+
  theme(text = element_text(size = 9))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("../submission/plots/eta.png", width = 6, height = 4)



#Repeat with the appropiate Eta values
estimators_100_eta <- confidence_add(est_100, .95,5)

plot_int(estimators_100_eta,Matrix_types,"gauss",100, type = "est")
ggsave("../submission/plots/ci_plots_eta.png",width = 6, height =3)
#Make full results table
tab_ci_eta <- confidence_table(estimators_100_eta,100 )
tab_ci_eta[,7:11] <- lapply(tab_ci_eta[,7:11], function(x)
  paste0("\\cellcolor{", colortab(x), "}", ifelse(is.na(x),0,x)))
#Add color gradient to the table
tab_ci_eta%>%
  filter(`Sampling Method` != "SRHT")%>%
  select(-Sample_Size) %>%
  rename(`$>10$` = `>10`,`$>1e-5$` = `>1e-5`,`$>1e-10$` = `>1e-10`,`$>0$` = `>0`)%>%
  xtable(type = "latex",digits=c(0,0,0,0,0,0,4,4,4,4,4))%>%
  print(tabular.environment="longtable",file = "../submission/plots/results_table_eta.tex",
        sanitize.text.function = identity, include.rownames=FALSE)

#Make results table for specific Matrix type
MAT = "gravity"

tab_ci_eta %>%
  filter(`Sampling Method` != "SRHT")%>%
  filter(`Matrix Type` == MAT)%>%
  select(-eta, -Sample_Size, -Width, -Rows, -`Matrix Type`)%>%
  rename(`$>10$` = `>10`,`$>1e-5$` = `>1e-5`,`$>1e-10$` = `>1e-10`,`$>0$` = `>0`)%>%
  xtable(type = "latex",digits=c(0,0,4,4,4,4,4))%>%
  print(tabular.environment="longtable",file = "../submission/plots/exp_table_eta.tex",
        sanitize.text.function = identity, include.rownames=FALSE)

#Make bar plot for stopping point errors
estimators_100_eta <- stopping_add(estimators_100_eta, 100, 1.1,.9,.01,.01)
stopping_graph(estimators_100_eta)
ggsave("../submission/plots/stop_true_eta.png",width = 4, height = 2 )
plot_stop_no_delt(estimators_100_eta)
ggsave("../submission/plots/no_delt_stop_eta.png",width = 4, height = 2 )

rep1_wilk_eta <- confidence_rep(rep_wilk, 15, .95, 13)
rep1_golub_eta <- confidence_rep(rep_golub, 15, .95, 13)
rep1_rohess_eta <- confidence_rep(rep_rohess, 15, .95, 13)
sum_wilk_eta <- get_all_obs_us(rep_wilk, rep1_wilk_eta, 15)
sum_golub_eta <- get_all_obs_us(rep_golub, rep1_golub_eta, 15)
sum_rohess_eta <- get_all_obs_us(rep_rohess, rep1_rohess_eta, 15)

sum_wilk_eta %>%
  filter(Meand <= Upperd & Meand >= Lowerd) %>%
  group_by(it) %>%
  summarise(max = max(Meand), min = min(Meand), up = max(Upperd), low = min(Lowerd))%>%
  write_csv("../article/csvs/wilkrepinsideeta.csv")

sum_wilk_eta %>%
  filter(Meand > Upperd | Meand < Lowerd) %>%
  select(it,Meand) %>%
  write_csv("../article/csvs/wilkrepoutsideeta.csv")

sum_golub_eta %>%
  filter(Meand <= Upperd & Meand >= Lowerd) %>%
  group_by(it) %>%
  summarise(max = max(Meand), min = min(Meand), up = max(Upperd), low = min(Lowerd))%>%
  write_csv("../article/csvs/golubrepinsideeta.csv")

sum_golub_eta %>%
  filter(Meand > Upperd | Meand < Lowerd) %>%
  select(it,Meand) %>%
  write_csv("../article/csvs/golubrepoutsideeta.csv")

sum_rohess_eta %>%
  filter(Meand <= Upperd & Meand >= Lowerd) %>%
  group_by(it) %>%
  summarise(max = max(Meand), min = min(Meand), up = max(Upperd), low = min(Lowerd))%>%
  write_csv("../article/csvs/rohessrepinsideeta.csv")

sum_rohess_eta %>%
  filter(Meand > Upperd | Meand < Lowerd) %>%
  select(it,Meand) %>%
  write_csv("../article/csvs/rohessrepoutsideeta.csv")


make_graphs(sum_wilk_eta)
ggsave("../submission/plots/rep_wilk_eta.png",width = 6, height =3 )
make_graphs(sum_golub_eta)
ggsave("../submission/plots/rep_golub_eta.png",width = 6, height =3 )
make_graphs(sum_rohess_eta)
ggsave("../submission/plots/rep_rohess_eta.png",width = 6, height = 3 )

#Calculate failure rates
sum_wilk_eta %>%
  ungroup()%>%
  summarise(failure_rate = mean(Meand > Upperd | Meand < Lowerd))
sum_golub_eta %>%
  ungroup()%>%
  summarise(failure_rate = mean(Meand > Upperd | Meand < Lowerd))
sum_rohess_eta %>%
  ungroup()%>%
  summarise(failure_rate = mean(Meand > Upperd | Meand < Lowerd))

#Create Intro Plot with lines
estimators_100%>%
     filter(Sampling_Method == "gauss", Matrix_Type == "randsvd", True_Res > 0, Iteration == 1 | Iteration %% 20 == 0) %>% 
  write_csv("../metrics-leastsquares/GenCoD/Article/csvs/randsvd_intro.csv")


estimators_100%>%
  filter(Sampling_Method == "gauss", Matrix_Type == "randsvd", True_Res > 0, Iteration == 1 | Iteration %% 20 == 0) %>%
  mutate(Bin = case_when(Error < 1150 & Error > 1100 ~ "1100 - 1150",
                         Error < 1100 & Error > 1050 ~ "1050 - 1100",
                         Error < 1050 & Error > 1000 ~ "1000 - 1050",
                         Error < 1000 & Error > 950 ~ "950 - 1000",
                         Error < 950 & Error > 900 ~ "900 - 950",
                         Error < 900 & Error > 850 ~ "850 - 900",
                         Error < 850 & Error > 800 ~ "800 - 850",
                         TRUE ~ "other"))%>%
  group_by(Bin)%>%
  summarise(O5 = quantile(Obs_Res,.05), O95 = quantile(Obs_Res,.95), T5 = quantile(True_Res, .05), T95 = quantile(True_Res,.95))

estimators_100%>%
  filter(Sampling_Method == "gauss", Matrix_Type == "randsvd", True_Res > 0, Iteration == 1 | Iteration %% 20 == 0) %>%
  mutate(Bin = case_when(Obs_Res < 40 & Obs_Res > 1 ~ "1 - 40",
                         Obs_Res < 1 & Obs_Res > .5 ~ ".5 - 1",
                         Obs_Res < .5 & Obs_Res > .1  ~ ".1 - .5",
                         Obs_Res < .1 & Obs_Res > .05 ~ ".05 - .1",
                         Obs_Res < .05 & Obs_Res > .01 ~ ".01 - .05",
                         Obs_Res < .01 & Obs_Res > .005 ~ ".005 - .01",
                         Obs_Res < .005 & Obs_Res > .001 ~ ".001 - .005",
                         TRUE ~ "other"))%>%
  group_by(Bin)%>%
  summarise(E5 = quantile(Error,.05), E95 = quantile(Error,.95))

estimators_100%>%
  filter(Sampling_Method == "gauss", Matrix_Type == "randsvd", True_Res > 0, Iteration == 1 | Iteration %% 20 == 0) %>%
  mutate(Bin = case_when(True_Res < 40 & True_Res > 1 ~ "1 - 40",
                         True_Res < 1 & True_Res > .5 ~ ".5 - 1",
                         True_Res < .5 & True_Res > .1  ~ ".1 - .5",
                         True_Res < .1 & True_Res > .05 ~ ".05 - .1",
                         True_Res < .05 & True_Res > .01 ~ ".01 - .05",
                         True_Res < .01 & True_Res > .005 ~ ".005 - .01",
                         True_Res < .005 & True_Res > .001 ~ ".001 - .005",
                         TRUE ~ "other"))%>%
  group_by(Bin)%>%
  summarise(E5 = quantile(Error,.05), E95 = quantile(Error,.95))


#Make a better Error plot 
lower = seq(0,150,by = 1)
upper = seq(-1,151,by = .5)
p = estimators_100%>%
  filter(Sampling_Method == "gauss", Matrix_Type == "randsvd", True_Res > 0, Iteration == 1 | Iteration %% 20 == 0) %>%
  mutate(error_bins = cut(Error,breaks = upper)) %>%
  group_by(error_bins)%>%
  filter(n() != 1, Width == 100)%>%
  unique() %>%
  ungroup()%>%
  #mutate(lower_bin = str_extract([0-9]{1,3}))%>%
  rename("Sketch Res" = Obs_Res,  "True Res"= True_Res, "Mov True" = Mov_True_Res) %>%
  pivot_longer(cols = c( `True Res`, `Mov True`), names_to = "Type", values_to = "Residual")%>%
  #boxplot(data=.,  Gradient~error_bins*Type,log = "y")
  ggplot()+
    geom_boxplot(aes(x = error_bins, y= Residual, col = Type))+
    scale_y_log10()+
    scale_color_manual(values = c("red","blue"))+
    theme_classic()+
    xlab("Absolute Error Bins")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = c(.95,.2),
        legend.box.background = element_rect(color = "white"),
        legend.title = element_blank())+
  labs(col = "")+
  ggtitle("Box plots for sketched and true residuals at different error ranges")
    
tikz(file = "~/Documents/Research_stuff/Linear_Systems/metrics-linearsystems-v2/article/plots/errorhist.tex", 
     width = 5, height = 3)
print(p)
dev.off()
#      
#      pivot_longer(cols =  c(True_Res,Mov_True_Res), names_to = "Method", values_to = "norm2") %>%
#      mutate(val = case_when(Method == "True_Res"~2,
#                                                        Method == "Mov_True_Res"~1,
#                                                        TRUE~3), Method = reorder(Method,-val)) %>%
#      ggplot(aes(x = Iteration, y = norm2, color = Method))+
#      geom_line()+
#      scale_x_log10()+
#      scale_y_log10()+
#      theme_bw()+
#      scale_color_manual(labels = c("True Residual","Moving Average"), values = c("green","blue"))+
#      ylab("Residual Norm Squared")+
#      ggtitle("Norm squared residuals compared to moving average of \nthose residuals at different iterations")+
#   theme(text = element_text(size = 18))+ 
#   theme(legend.position="bottom")
# ggsave("../submission/plots/Moler_true_res_gauss.png",width = 6, height = 4)
#Create  Intro plot with errors
estimators_100%>%
     filter(Sampling_Method == "gauss", Matrix_Type == "randsvd",True_Res >0 )
  
#     mutate(True_Res.x = True_Res, Mov_True_Res.x = Mov_True_Res, Obs_Res.x = Obs_Res) %>%
#      pivot_longer(cols =  c(True_Res,Mov_True_Res), names_to = "Method", values_to = "norm2") %>%
#      mutate(val = case_when(Method == "True_Res"~2,
#                                                        Method == "Mov_True_Res"~1,
#                                                        TRUE~3), Method = reorder(Method,-val)) %>%
#      ggplot(aes(x = Error, y = norm2, color = Method))+
#      geom_point()+
#      geom_point(aes(x = Error, y = Obs_Res.x),color = "violet", alpha = .2)+ 
#      geom_point(aes(x = Error, y = True_Res.x), color  = "green", alpha = .2)+
#      geom_point(aes(x = Error, y = Mov_True_Res.x),color = "blue", alpha = .2)+
#      scale_x_log10()+
#      scale_y_log10()+
#      theme_bw()+
#      scale_color_manual(labels = c("True Residual","Moving Average"), values = c("green","blue"))+
#      ylab("Residual Norm Squared")+
#   ggtitle("Norm squared residuals compared to moving average of \nthose residuals at different errors")+
#   theme(text = element_text(size = 18))+ 
#   theme(legend.position="bottom")
# ggsave("../submission/plots/Moler_true_res_gauss_point.png",width = 6, height = 4)
