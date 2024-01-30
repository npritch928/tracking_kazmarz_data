library(xtable)
library(tidyverse)
library(cowplot)
#Creates Plot of Highest lowest and average difference between actual and predicted rho
#Requires input of results after estimators has been called

plot_rho <- function(results, rel = F){
  if(!rel){
    g <- results %>%
      ungroup()%>%
      group_by(Matrix_Type, Sampling_Method)%>%
      mutate(Conv = ifelse(min(True_Res) < 1e-08,"Converged","Failed to Converge")) %>%
      ungroup()%>%
      mutate(diff_rho = abs(Mov_Obs_Res - Mov_True_Res))%>%
      group_by(Iteration,Conv)%>%
      summarise(`95th` = quantile(diff_rho,.95), `5th` = quantile(diff_rho,.05), `50th` = median(diff_rho))
      # pivot_longer(cols = c(`95th`,`5th`,`50th`), names_to = "Percentile", values_to = "Absolute_diff")%>%
      # ungroup()%>%
      # mutate(per = as.numeric(str_extract(Percentile,"[0-9]+")), Percentile = reorder(Percentile,per))%>%
      # ggplot()+
      # geom_line(aes(x = Iteration, y = Absolute_diff, col = Percentile))+
      # facet_wrap(~Conv)+
      # scale_y_log10()+
      # scale_color_manual(values = c("blue","purple","red"))+
      # ylab("Absolute Difference")+
      # ggtitle("Absolute difference between true and estimated\nmoving average of squared residuals")+
      # theme_bw()+
      # theme(text = element_text(size = 9))+
      # theme(legend.position = "bottom")
  }else{
    g <- results %>%
      ungroup()%>%
      group_by(Matrix_Type, Sampling_Method)%>%
      mutate(Conv = ifelse(min(True_Res) < 1e-08,"Converged","Failed to Converge")) %>%
      ungroup()%>%
      filter(Mov_True_Res != 0) %>%
      mutate(diff_rho = abs(Mov_Obs_Res - Mov_True_Res)/Mov_True_Res)%>%
      group_by(Iteration,Conv)%>%
      summarise(`95th` = quantile(diff_rho,.95), `5th` = quantile(diff_rho,.05), `50th` = median(diff_rho))
      # pivot_longer(cols = c(`95th`,`5th`,`50th`), names_to = "Percentile", values_to = "Absolute_diff")%>%
      # ungroup()%>%
      # mutate(per = as.numeric(str_extract(Percentile,"[0-9]+")), Percentile = reorder(Percentile,per))%>%
      # ggplot()+
      # geom_line(aes(x = Iteration, y = Absolute_diff, col = Percentile))+
      # facet_wrap(~Conv)+
      # scale_y_log10()+
      # scale_color_manual(values = c("blue","purple","red"))+
      # ylab("Relative Difference")+
      # ggtitle("Relative difference between true and estimated\nmoving average of squared residuals")+
      # theme_bw()+
      # theme(text = element_text(size = 9))+
      # theme(legend.position = "bottom")
  }
  return(g)
}

#Creates Plot of Highest lowest and average difference between actual and predicted iota
#Requires input of results after estimators has been called

plot_iota <- function(results, rel = F){
  if(! rel){
    g <- results %>%
      ungroup()%>%
      group_by(Matrix_Type, Sampling_Method)%>%
      mutate(Conv = ifelse(min(True_Res) < 1e-08,"Converged","Failed to Converge")) %>%
      ungroup()%>%
      mutate(diff_iota = abs(Mov_Obs_Res2 - Mov_True_Res2))%>%
      group_by(Iteration,Conv)%>%
      summarise(`95th` = quantile(diff_iota,.95), `5th` = quantile(diff_iota,.05), `50th` = median(diff_iota))
      # pivot_longer(cols = c(`95th`,`5th`,`50th`), names_to = "Percentile", values_to = "Absolute_diff")%>%
      # ungroup()%>%
      # mutate(per = as.numeric(str_extract(Percentile,"[0-9]+")), Percentile = reorder(Percentile,per))%>%
      # ggplot()+
      # geom_line(aes(x = Iteration, y = Absolute_diff, col = Percentile))+
      # facet_wrap(~Conv)+
      # scale_y_log10()+
      # scale_color_manual(values = c("blue","purple","red"))+
      # ylab("Absolute Difference")+
      # ggtitle("Absolute difference between true and estimated\nmoving average of residuals to the fourth power")+
      # theme_bw()+
      # theme(text = element_text(size = 9))+
      # theme(legend.position = "bottom")
  }else{
    g <- results %>%
      ungroup()%>%
      group_by(Matrix_Type, Sampling_Method)%>%
      mutate(Conv = ifelse(min(True_Res) < 1e-08,"Converged","Failed to Converge")) %>%
      ungroup()%>%
      mutate(diff_iota = abs(Mov_Obs_Res2 - Mov_True_Res2)/Mov_True_Res2)%>%
      mutate(diff_iota = ifelse(is.na(diff_iota),0,diff_iota)) %>%
      filter(diff_iota != 0)%>%
      group_by(Iteration,Conv)%>%
      summarise(`95th` = quantile(diff_iota,.95), `5th` = quantile(diff_iota,.05), `50th` = median(diff_iota))#%>%
      #pivot_longer(cols = c(`95th`,`5th`,`50th`), names_to = "Percentile", values_to = "Absolute_diff")%>%
      #ungroup()%>%
      #mutate(per = as.numeric(str_extract(Percentile,"[0-9]+")), Percentile = reorder(Percentile,per))%>%
      #ggplot()+
      #geom_line(aes(x = Iteration, y = Absolute_diff, col = Percentile))+
      #facet_wrap(~Conv)+
      #scale_y_log10()+
      #scale_color_manual(values = c("blue","purple","red"))+
      #ylab("Relative Difference")+
      #ggtitle("Relative difference between true and estimated\nmoving average of residuals to the fourth power")+
      #theme_bw()+
      #theme(text = element_text(size = 9))+
      #theme(legend.position = "bottom")
  }
  
  return(g)
}
#Function designed to plot the true and estimated confidence intervals
#type can be either both, true, or est
plot_int <- function(results, Matrix_type, Sampling_meth, width, type = "both"){
  name <- str_c( "Credible intervals for estimators of width ",width, " for ",ifelse(Sampling_meth == "gauss","Gaussian",Sampling_meth), ' sampling')
  if(type == "both"){
    g <- results %>%
      filter(Sampling_Method == Sampling_meth, 
             Matrix_Type %in% Matrix_type)%>%
      mutate(Sampling_Method = str_to_title(Sampling_Method))%>%
      pivot_longer(cols = c(Mov_Obs_Res, Mov_True_Res, True_Res), 
                   names_to = "Residual Type",
                   values_to = "Value") %>%
      mutate(           Ordering = case_when(`Residual Type` == "Mov_Obs_Res" ~ 3,
                                             `Residual Type` == "Mov_True_Res" ~ 2,
                                             `Residual Type` == "True_Res" ~ 1),
                        `Residual Type` = case_when(`Residual Type` == "Mov_Obs_Res" ~ "Average Observed",
                                         `Residual Type` == "Mov_True_Res" ~ "Average True",
                                         `Residual Type` == "True_Res" ~ "True Residual"),
                        `Residual Type` = reorder(`Residual Type`,Ordering))%>%
      filter(Ordering == 1)%>%
      select(-Ordering)%>%
      unite(Estimated,Lower_est,Upper_est)%>%
      unite(True,Lower_true,Upper_true)%>%
      pivot_longer(cols = c(Estimated,True), 
                   names_to = "Interval Type",
                   values_to = "int") %>%
      separate(int,into = c("Lower","Upper"),convert=T,sep ="_") %>%
      ggplot()+
      geom_ribbon(aes(x = Iteration, ymin = Lower, ymax = Upper, fill = `Interval Type`), alpha = .2)+
      geom_line(aes(x = Iteration, y = Value, col = `Residual Type`))+
      facet_wrap(~Matrix_Type)+
      #scale_y_log10()+
      scale_x_log10()+
      ylab("Residual Norm Squared")+
      scale_fill_manual(values = c("black","purple"))+
      ggtitle(name)+
      theme_bw()+
      theme(text = element_text(size = 9))+
      theme(legend.position = "bottom")
  } else if(type == "est"){
    g <- results %>%
      filter(Sampling_Method == Sampling_meth, 
             Matrix_Type %in% Matrix_type)%>%
      mutate(Sampling_Method = str_to_title(Sampling_Method))%>%
      pivot_longer(cols = c(Mov_Obs_Res, Mov_True_Res, True_Res), 
                   names_to = "Residual Type",
                   values_to = "Value") %>%
      mutate(           Ordering = case_when(`Residual Type` == "Mov_Obs_Res" ~ 3,
                                             `Residual Type` == "Mov_True_Res" ~ 2,
                                             `Residual Type` == "True_Res" ~ 1),
                        `Residual Type` = case_when(`Residual Type` == "Mov_Obs_Res" ~ "Average Observed",
                                                    `Residual Type` == "Mov_True_Res" ~ "Average True",
                                                    `Residual Type` == "True_Res" ~ "True Residual"),
                        `Residual Type` = reorder(`Residual Type`,Ordering))%>%
      filter(Ordering != 1)%>%
      select(-Ordering)%>%
      ggplot()+
      geom_ribbon(aes(x = Iteration, ymin = Lower_est, ymax = Upper_est), fill = "black", alpha = .3)+
      geom_line(aes(x = Iteration, y = Value, col = `Residual Type`))+
      scale_color_manual(values = c("red","blue"))+
      facet_wrap(~Matrix_Type)+
      #scale_y_log10()+
      scale_x_log10()+
      ylab("Residual Norm Squared")+
      ggtitle(name)+
      theme_bw()+
      theme(text = element_text(size = 9))+
      theme(legend.position = "bottom")
    
  }else {
    g <- results %>%
      filter(Sampling_Method == Sampling_meth, 
             Matrix_Type %in% Matrix_type)%>%
      mutate(Sampling_Method = str_to_title(Sampling_Method))%>%
      pivot_longer(cols = c(Mov_Obs_Res, Mov_True_Res, True_Res), 
                   names_to = "Residual Type",
                   values_to = "Value") %>%
      filter(Ordering == 1) %>%
      mutate(           Ordering = case_when(`Residual Type` == "Mov_Obs_Res" ~ 3,
                                             `Residual Type` == "Mov_True_Res" ~ 2,
                                             `Residual Type` == "True_Res" ~ 1),
                        `Residual Type` = case_when(`Residual Type` == "Mov_Obs_Res" ~ "Average Observed",
                                                    `Residual Type` == "Mov_True_Res" ~ "Average True",
                                                    `Residual Type` == "True_Res" ~ "True Residual"),
                        `Residual Type` = reorder(`Residual Type`,Ordering))%>%
      filter(Ordering != 1)%>%
      select(-Ordering)%>%
      ggplot()+
      geom_ribbon(aes(x = Iteration, ymin = Lower_true, ymax = Upper_true), fill = "purple", alpha = .2)+
      geom_line(aes(x = Iteration, y = Value, col = `Residual Type`))+
      facet_wrap(~Matrix_Type)+
      #scale_y_log10()+
      scale_x_log10()+
      ylab("Residual Norm Squared")+
      ggtitle(name)+
      theme_bw()+
      theme(text = element_text(size = 9))+
      theme(legend.position = "bottom")
    
  }
    return(g)
}

#Function counting number of times interval is exceeded by true moving 
#average 
count_outside <- function(results){
  p<-results %>% 
    mutate(Out_int = ifelse((Mov_True_Res > Upper_est | Mov_True_Res < Lower_est),T,F), 
           Res_Lev = case_when(Mov_True_Res > 10 ~ "g10",
                               Mov_True_Res <= 10 & Mov_True_Res > 1e-5 ~ "n5",
                               Mov_True_Res <= 1e-5 & Mov_True_Res > 1e-10 ~"n10",
                               Mov_True_Res <= 1e-10 ~ "z10"))%>%
    group_by(Matrix_Type, Rows, Sampling_Method, eta, Res_Lev)%>%
    summarise(Iteration = n(), Count_Outside = sum(Out_int)) %>%
    ungroup() %>% 
    group_by(Matrix_Type, Rows, Sampling_Method, eta) %>%
    mutate(Total = sum(Count_Outside)/sum(Iteration), Count_Outside = Count_Outside / sum(Iteration))%>%
    pivot_wider(id_cols = c(Matrix_Type, Rows, Sampling_Method, eta, Total), 
                names_from = Res_Lev,
                values_from = Count_Outside,
                values_fill = 0)
  cols = c("g10","n5","n10","z10")
  cn <- cols[!cols %in% colnames(p)]
  for (i in seq_along(cn)) {
    p<-p%>%mutate(!!cn[i] := 0)
  }
  p<-p%>%mutate( "Sample_Size" =  results$Sample_size[1]) %>% 
    mutate(g10 = round(g10,4), n5 = round(n5,4), n10 = round(n10,4), z10 = round(z10,4), Total = round(Total,4))%>%
    rename(">10" = g10,
           ">1e-5" = n5,
           ">1e-10" = n10,
           ">0" = z10) %>%
    select(Matrix_Type, Rows, Sampling_Method, Sample_Size, eta, `>10`,`>1e-5`,`>1e-10`,`>0`,Total)
  return(p)
}
#Function that does a spectrum of red coloring for table
colortab<-function(x){
  ifelse(is.na(x),"white",as.character(
    cut(x,c(-.9999,.05,.075,.1,.5,Inf),
        labels = c("white","red1","red2","red3","red4"),
        right = T))
  )
}
#Function to generate a table that records the number of times the confidence 
#interval fails to cover the true moving average
confidence_table<-function(results, width = 100)
{
  tab_ci <- count_outside(results)
  tab_ci <- tab_ci%>%
    mutate(Width = width) %>%
    select(Sampling_Method,Width,everything())%>%
    ungroup()%>%
    arrange(Sampling_Method,Matrix_Type,Width)%>%
    rename("Sampling Method" = Sampling_Method, "Matrix Type" = Matrix_Type)

  return(tab_ci)
}
#Function that combines functions to create the eta table
eta_test<-function(results,alpha,width,eta){
  res <- confidence_add(results,alpha,eta)
  ct <- confidence_table(res,width)
  return(ct)
}
#Make a bar graph showing the distribution of the outcomes of stopping
# Assumes that stopping_add has been run prior 
stopping_graph<-function(results){
  name <- str_c("Error from stopping criterion")
  
  g <- results %>%
    filter(Mov_Obs_Res2 <= threshold) %>%
    group_by(Matrix_Type, Sampling_Method, Sample_size, Rows) %>%
    #mutate(class = case_when(stop & Mov_Obs_Res > nu * delt1 ~ "t1",
     #                        !stop & Mov_Obs_Res < nu * delt2 ~ "t2",
      #                       stop & Mov_Obs_Res <= nu * delt1 | !stop & Mov_Obs_Res >= nu * delt2~"correct"))
    summarise(n = n(), 
              t1 = sum(stop & Mov_True_Res > nu * delt1), 
              t2 = sum(!stop & Mov_True_Res <= nu * delt2), 
              correct = sum(stop & Mov_True_Res <= nu * delt1 | !stop & Mov_True_Res > nu * delt2))%>%
    ungroup()%>%
    summarise(total = sum(n), 
              t1 = sum(t1)/total, 
              correct = sum(correct)/total, 
               t2 = sum(t2)/total) 
    # pivot_longer(cols = c(t1,t2,correct),
    #              names_to = "Error",
    #              values_to = "Proportion")%>%
    # mutate(Error = case_when(
    #   Error == "t1" ~ "Type I",
    #   Error == "t2" ~ "Type II",
    #   Error == "correct" ~ "No Error"
    # ))%>%
    # ggplot()+
    # geom_col(aes(x = Error, y = Proportion, fill = Error),show.legend = F)+
    # theme_bw()+
    # scale_fill_manual(values = c("green", "red", "blue"))+
    # ggtitle(name)
  return(g)
}

plot_stop_no_delt<-function(results){
  thr <- results %>%
    mutate(c1 = ifelse(Mov_Obs_Res < nu, T,F), correct = ifelse(Mov_True_Res < nu, T,F))%>%
    filter(c1 == T | correct == T) %>% 
    filter(Mov_Obs_Res2 < threshold) %>%
    group_by(Matrix_Type,Sampling_Method) %>% 
    mutate(min= min(Iteration))%>%
    filter(min == Iteration)%>%
    mutate(Error = case_when(c1 == T & correct == T~"Correct",
                             c1 == F & correct == T ~"Type I",
                             c1 == T & correct == F~"Type II"))%>%
    group_by(Error)%>%
    summarise(n = n()) %>% 
    ungroup() %>% 
    mutate(total = sum(n), Rate = n/total)%>%
    mutate(Type = "Thresholding")
  
  nthr <- results %>%
    mutate(c1 = ifelse(Mov_Obs_Res < nu, T,F), correct = ifelse(Mov_True_Res < nu, T,F))%>%
    filter(c1 == T | correct == T) %>% 
    group_by(Matrix_Type,Sampling_Method) %>% 
    mutate(min= min(Iteration))%>%
    filter(min == Iteration)%>%
    mutate(Error = case_when(c1 == T & correct == T~"Correct",
                             c1 == F & correct == T ~"Type I",
                             c1 == T & correct == F~"Type II"))%>%
    group_by(Error)%>%
    summarise(n = n()) %>% 
    ungroup() %>% 
    mutate(total = sum(n), Rate = n/total)%>%
    mutate(Type = "No Thresholding")
  
  g<-rbind(thr,nthr)%>%
    ggplot(aes(x = Error, y = Rate, fill = Error)) + 
    geom_col()+
    theme_bw()+
    facet_wrap(~Type)+
    scale_fill_manual(values = c("green", "red", "blue"))+
    ggtitle("Stopping criteria error rates")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
  
  return(g)
}
