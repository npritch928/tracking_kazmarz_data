iter_plot = function(Sampling_meth,stride,Matrix_types, Matrix_size, Samp_size, Average_Width, results)
{
  name <- str_c(ifelse(Sampling_meth == "gauss","Gauss",Sampling_meth), " Sampling For Matrix of size ",  Matrix_size,"\n moving average of width ", Average_Width)
  p <- results%>%
    filter(Matrix_Type %in% Matrix_types,
           Width == Average_Width,
           Sampling_Method == Sampling_meth, 
           Sample_size == Samp_size, 
           Iteration %% stride == 0 || Iteration == 1, True_Res > 5e-16)%>%
    pivot_longer(cols = c("Obs_Res","Mov_Obs_Res","Mov_True_Res","True_Res","Error"),
                 names_to = "Method", 
                 values_to = "Norm")%>%
    filter(!(Method %in% c("Error","Obs_Res")))
  p$Method<-ifelse(p$Method == "Mov_True_Res","zTrue_Res_Mov",p$Method)
  p$Method<-ifelse(p$Method == "Mov_Obs_Res","zMov_Obs_res",p$Method)
  g<- p%>%ggplot()+
    geom_line(aes(x = Iteration, y = Norm, col = Method))+
    facet_wrap(.~Matrix_Type)+
    scale_y_log10()+
    scale_x_log10()+
    theme_bw()+
    ggtitle(name)+
    scale_color_manual(labels = c( "True Residual","Average Observed Residual","Average True Residual"),
                       values = c("green","red","blue"))+
    ylab("Norm - Squared")
  return(g)
}
# function to compute the moving average 
Moving_av <- function(width,values,i)
{
  if(i > width){
    mean(values[(i - width + 1):i])
  }else{
    mean(values[1:i])
  }
}
#Moving av two phases operates on assumption that in convergence stage mean 
#of miving average greater than true observations
Moving_av_phase<-function(width,values,i,start)
{
  if(is.na(start) & values[i] <= values[ifelse(i == 1, 1 ,i-1)] ){
    
    start <- NA
    return(c(start,values[i]))
    
  }else if(is.na(start) & values[i] > values[i-1]){
    start <- i
    return(c(start,values[i]))
    
  }else{
    
    if(i > width + start){
      me<-mean(values[(i - width + 1):i])
    }else{
      me<-mean(values[start:i])
    }
    return(c(start,me))
  }
  
}


# Produces a confidence interval
confidence_band <- function(Sampling_meth, delta, Sample_size, Matrix_size, Obs_res, width, Conf_Lev,eta = 4 * width)
{
  Obs_four <- Obs_res^2
  if(Sampling_meth == "Achlio" | Sampling_meth == "gauss"){
    c<- 2.12258 * (1+log(width)) / (Sample_size*width)
    omega <- .1127
  }else{
    c <- 32 * (1 + log(width))/ (Sample_size*width)
    omega <- 1/16
  }
  iterates <- length(Obs_res)
  point <- numeric(iterates)
  vrnce <- numeric(iterates)
  #Assuming even matrix has the same number of iterates
  if(Sample_size != 1){
    for(i in seq_len(iterates)){
      point[i] <- Moving_av(width,Obs_res,i)
      vrnce[i] <- c * Moving_av(width,Obs_four,i)
    }
  }else{
    NULL
  }
  tstar <-(-log((1-Conf_Lev)/2) * 2 * vrnce / (omega*eta))
  tstart <- ifelse(tstar > 1, tstar, sqrt(omega*tstar))
  up <- point + tstart
  low <- point - tstart
  low <- ifelse(low <= 0,min(1e-16,min(point)*10^-2), low)
  return(tibble(Upper = up, Lower = low))
  #return(tibble(Upper = up))
}
#Createsconfidence  bands using 2 phase moving average
confidence_band2 <- function(Sampling_meth, delta, Sample_size, Matrix_size, Obs_res, Width, Conf_Lev,eta)
{
  
  Obs_four <- Obs_res^2
  if(Sampling_meth == "Achlio" | Sampling_meth == "gauss"){
    c<- 2.12258 * (1+log(width)) / (Sample_size*width)
    omega <- .1
  }else{
    c <- 32 *  (1 + log(width))/ (Sample_size*width)
    omega <- 1/16
  }
  iterates <- length(Obs_res)
  point <- numeric(iterates)
  vrnce <- numeric(iterates)
  if(Sample_size != 1){
    d <- c(NA,NA)
    for(i in seq_len(iterates)){
      width <- ifelse(is.na(d[1]),1,Width)
      d <- Moving_av_phase(width,Obs_res,i,d[1])
      point[i] <- d[2]
      vrnce[i] <- c * Moving_av_phase(width,Obs_four,i,d[1])[2]
    }
  }else{
    NULL
  }
  tstar <-(-log((1-Conf_Lev)/2) * 2 * vrnce / (omega*eta))
  tstart <- ifelse(tstar > 1, tstar, sqrt(omega*tstar))
  up <- point + tstart
  low <- point - tstart
  low <- ifelse(low <= 0, min(1e-16,min(point)*10^-2), low)
  return(tibble(Upper = up, Lower = low))
}

#Produces the plot of the confidence interval
Conf_plot <- function(Sampling_meth, epsilon, Sample_siz, Matrix_types,Matrix_size, width, Conf_Lev, stride,results,eta)
{
  Matrix_types <- Matrix_types[order(Matrix_types)]
  name <- str_c(ifelse(Sampling_meth == "gauss","Gauss",Sampling_meth), " sampling for matrix of size ",  Matrix_size,",\n moving average of width ", width, ",\n and ", 100*Conf_Lev, "% Confidence Interval")
  p <- results%>%
    filter(Sampling_Method == Sampling_meth, Sample_size == Sample_siz, Matrix_size == Rows,
           Matrix_Type %in% Matrix_types, Width == width, Iteration %% stride == 0 || Iteration == 1)%>%
    arrange(Matrix_Type)

  cb <- tibble()
  for (i in Matrix_types) {
    Obs_res <- p%>%select(Matrix_Type,Obs_Res) %>% 
      filter(Matrix_Type == i)  %>% 
      pull(Obs_Res)
    cb <- rbind(cb,
                confidence_band(Sampling_meth, 
                                epsilon, Sample_siz, 
                                Matrix_size, Obs_res, width, Conf_Lev,eta)
                )
    
  }
  
  p <- p%>%
      cbind(cb)
  
  p <- p%>%pivot_longer(cols = c("Obs_Res","Mov_Obs_Res","Mov_True_Res","True_Res","Error"),
                 names_to = "Method", 
                 values_to = "Norm")%>%
    filter(!(Method %in% c("Error","Obs_Res")))
  p$Lower <- ifelse(p$Method == "Mov_Obs_Res", p$Lower,p$Norm)
  p$Upper <- ifelse(p$Method == "Mov_Obs_Res", p$Upper,p$Norm)
  p$Method <- ifelse(p$Method == "Mov_True_Res","zTrue_Res_Mov",p$Method)
  p$Method <- ifelse(p$Method == "Mov_Obs_Res","zMov_Obs_res",p$Method)
  
  g<- p%>%ggplot(aes(x = Iteration, y = Norm))+
    geom_line(aes( col = Method))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = .3, linetype ="dotted")+
    facet_wrap(.~Matrix_Type)+
    scale_y_log10()+
    scale_x_log10()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1))+
    ggtitle(name)+
    scale_color_manual(labels = c( "True Residual","Average Observed Residual","Average True Residual"),
                       values = c("green","red","blue"))+
    scale_fill_manual(labels = c( "True Residual","Average Observed Residual","Average True Residual"),
                       values = c("green","red","blue"))+
   ylab("Norm - Squared")
  return(g)
}
#Builds a  two  phase confidence plot for demonstration purposes of effectiveness
Conf_plot2 <- function(Sampling_meth, epsilon, Sample_siz, Matrix_types,Matrix_size, width, Conf_Lev, stride,results,eta){
  name <- str_c(ifelse(Sampling_meth == "gauss","Gauss",Sampling_meth), " sampling for matrix of size ",  Matrix_size,",\n two phase moving average of width 1 and then ", width, ",\n and ", 100*Conf_Lev, "% Confidence Interval")
  e_col <- tibble(Lower = NULL, Upper = NULL, True_Mov_av_2 = NULL, Obs_Mov_av_2 = NULL)
  Matrix_types <- Matrix_types[order(Matrix_types)]
  results <- results %>%
    filter(Matrix_Type %in% Matrix_types)%>%
    arrange(Matrix_Type)
  for (j in Matrix_types) {
    pl <- results%>%
      filter(Sampling_Method == Sampling_meth, Sample_size == Sample_siz, Matrix_size == Rows,
             Matrix_Type == j , Width == width)
    it <- nrow(pl)
    Obs_res <- pl$Obs_Res
    True_res <- pl$True_Res
    True_Mov_av_2 <- numeric(it)
    Obs_Mov_av_2 <- numeric(it)
    d1 <- c(NA,NA)
    d2 <- c(NA,NA)
    for(i in seq_along(1:it)){
      d1 <- Moving_av_phase(width,Obs_res,i,d1[1])
      d2 <- Moving_av_phase(width,True_res,i,d2[1])
      True_Mov_av_2[i] <- d2[2]
      Obs_Mov_av_2[i] <- d1[2]
    }
  
    cb <- confidence_band2(Sampling_meth, epsilon, Sample_siz, Matrix_size, Obs_res, width, Conf_Lev,eta)
    e_col <- rbind(e_col,cbind(cb,True_Mov_av_2,Obs_Mov_av_2))
  }
  p <- results%>%
    filter(Sampling_Method == Sampling_meth, Sample_size == Sample_siz, Matrix_size == Rows,
           Matrix_Type %in% Matrix_types , Width == width)
  p <- p%>%
    cbind(e_col)

  p <- p%>%pivot_longer(cols = c("Obs_Res","Obs_Mov_av_2","True_Mov_av_2","True_Res","Error"),
                        names_to = "Method", 
                        values_to = "Norm")%>%
    filter(!(Method %in% c("Error","Obs_Res")))
  p$Lower <- ifelse(p$Method == "Obs_Mov_av_2", p$Lower,p$Norm)
  p$Upper <- ifelse(p$Method == "Obs_Mov_av_2", p$Upper,p$Norm)
  p$Method <- ifelse(p$Method == "True_Mov_av_2","zTrue_Res_Mov",p$Method)
  p$Method <- ifelse(p$Method == "Obs_Mov_av_2","zMov_Obs_res",p$Method)
  
  g<- p%>%ggplot(aes(x = Iteration, y = Norm))+
    geom_line(aes( col = Method))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = .3, linetype ="dotted")+
    facet_wrap(.~Matrix_Type)+
    scale_y_log10()+
    scale_x_log10()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1))+
    ggtitle(name)+
    scale_color_manual(labels = c( "True Residual","Average Observed Residual","Average True Residual"),
                       values = c("green","red","blue"))+
    scale_fill_manual(labels = c( "True Residual","Average Observed Residual","Average True Residual"),
                      values = c("green","red","blue"))+
    ylab("Norm - Squared")
  return(g)
}
#Function used to classify the true  residuals for the creation of the table
classify<-function(value){
  if(value > 10){
    return("g10")
  }else if(value <= 10 & value > 1e-5){
    return("n5")
  }else if(value <= 1e-5 & value > 1e-10){
    return("n10")
  }else if(value <= 1e-10 & value > 1e-15){
    return("n15")
  }else if(value <= 1e-15){
    return("l15")
  }
}

count_outside <- function(Sampling_meth, epsilon, width, Sample_siz, Matrix_types,Matrix_size, Conf_Lev,results, eta = 4*width){
  Matrix_types <- Matrix_types[order(Matrix_types)]
  p <- results%>%
    filter(Sampling_Method == Sampling_meth, Sample_size == Sample_siz, Matrix_size == Rows,
           Matrix_Type %in% Matrix_types, Width == width)%>%
    arrange(Matrix_Type)

  #Use adaptive window when width is specified to be 1.
  if(width != 1){
    cb <- tibble()
    for (i in Matrix_types) {
      Obs_res <- p%>%select(Matrix_Type,Obs_Res) %>% 
        filter(Matrix_Type == i)  %>% 
        pull(Obs_Res)
      cb <- rbind(cb,
                  confidence_band(Sampling_meth, 
                                  epsilon, Sample_siz, 
                                  Matrix_size, Obs_res, width, Conf_Lev, eta)
      )
      
    }
  
  }else{
    Obs_res <- p$Obs_Res
    cb <- confidence_band2(Sampling_meth, epsilon, Sample_siz, Matrix_size, Obs_res, width, Conf_Lev)
  }
  p <- p%>%
    cbind(cb)%>%
    mutate(Out_int = ifelse((Mov_True_Res > Upper | Mov_True_Res < Lower),T,F), 
           Res_Lev = case_when(Mov_True_Res > 10 ~ "g10",
                               Mov_True_Res <= 10 & Mov_True_Res > 1e-5 ~ "n5",
                               Mov_True_Res <= 1e-5 & Mov_True_Res > 1e-10 ~"n10",
                               Mov_True_Res <= 1e-10 ~ "z10"))%>%
    group_by(Matrix_Type, Width, Rows, Sampling_Method, Res_Lev)%>%
    summarise(Iteration = n(), Count_Outside = sum(Out_int)/Iteration) %>%
    pivot_wider(id_cols = c(Matrix_Type, Width, Rows, Sampling_Method), 
                names_from = Res_Lev,
                values_from = Count_Outside,
                values_fill = 0)
    cols = c("g10","n5","n10","z10")
    cn <- cols[!cols %in% colnames(p)]
    for (i in seq_along(cn)) {
      p<-p%>%mutate(!!cn[i] := 0)
    }
    p<-p%>%mutate(Total = sum(c(g10,n5 ,n10,z10),na.rm = T), "eta" = eta, "Sample_Size" =  Sample_siz) %>% 
      mutate(g10 = round(g10,4), n5 = round(n5,4), n10 = round(n10,4), z10 = round(z10,4), Total = round(Total,4))%>%
      rename(">10" = g10,
           ">1e-5" = n5,
           ">1e-10" = n10,
           ">0" = z10) %>%
      select(Matrix_Type, Width, Rows, Sampling_Method, Sample_Size, eta, `>10`,`>1e-5`,`>1e-10`,`>0`,Total)
  return(p)
}

CI_rel_true <- function(Sampling_meth, epsilon, width, Sample_siz, Matrix_types,Matrix_size, Conf_Lev,results){
  p <- results%>%
    filter(Sampling_Method %in% Sampling_meth, Sample_size == Sample_siz, Matrix_size == Rows,
           Matrix_Type %in% Matrix_types, Width == width)
  Obs_res <- p$Obs_Res
  cb <- confidence_band(Sampling_meth, epsilon, Sample_siz, Matrix_size, Obs_res, width, Conf_Lev)
  p <- p%>%
    cbind(cb)%>%
    mutate(Out_int = ifelse(True_Res > Upper | True_Res < Lower,T,F), Int_Width = Upper - Lower, distance = abs(Mov_Obs_Res - Mov_True_Res), relative_distance = abs(Int_Width - distance)/distance)%>%
    group_by(Iteration)%>%
    summarise(relative_distance = mean(relative_distance))
  g <- p %>%ggplot(aes(x = Iteration, y = relative_distance))+
    geom_line(aes( col = Sampling_meth))+
    #facet_wrap(.~Matrix_Type)+
    scale_y_log10()+
    scale_x_log10()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1))+
    #ggtitle(name)+
    scale_color_manual(labels = c( "True Residual","Average True Residual"),
                       values = c("green","blue"))+
    scale_fill_manual(labels = c( "True Residual","Average True Residual"),
                      values = c("green","blue"))+
    ylab("Norm - Squared")
  return(g)
}