draw_plot <- function(df,y_name,y_lab){
  y_name <- sym(y_name)
  p2 <- ggplot(data = df , aes(y = !!y_name, x = dt)) +
    geom_point()  + ylab(y_lab) + xlab("Date") +
    geom_line(aes(y = predC, x = dt, color="counterfactual",linetype = "counterfactual"), size=1)+
    geom_line(aes(y = pred, x = dt, color="fitted values",linetype = "fitted values"), size=1)+
    scale_color_manual("Legend", values = c("counterfactual"="red","fitted values"="black")) +
    scale_linetype_manual("Legend",values = c("counterfactual"=2,"fitted values"=1)) +
    theme(axis.text.x = element_text(angle = 90))+
    geom_vline(xintercept = as.Date("2020-03-01", "%Y-%m-%d"), color = "blue") +
    geom_vline(xintercept = as.Date("2020-09-01", "%Y-%m-%d"), color = "green") +
    geom_vline(xintercept = as.Date("2021-01-01", "%Y-%m-%d"), color = "yellow") +
    theme(legend.key=element_blank())
  return(p2)
}


est_df <- function(model){
   df <- data.frame(matrix(NA,8,2))
   CI <- confint(model)
   sum_df <- data.frame(summary(model)[[12]])
   sum_df <- sum_df[,c(1,4)]
   est <- round(sum_df$Estimate,digits = 3)
   lower <- round(CI[,1], digits=2)
   upper<- round(CI[,2],digits=2)
   df[,1] <- paste0(est," (",lower,", ",upper,")")
   df[,2] <- round(sum_df[,2],digits=2)
   ind1 <- which(df[,2]<0.001,arr.ind = TRUE)
   df[ind1,2] <- "P<0.001"
   df <- df[c(1,6:8),]
   rownames(df) <- c("Interecpt","Time","COVID-19 period","Time: COVID-19 period")
   colnames(df) <- c("Estimate (95% CI range)", "p-value")
   return(df)
 }

est_df_RR <- function(model){
   df <- data.frame(matrix(NA,8,2))
   CI <- confint(model)
   sum_df <- data.frame(summary(model)[[12]])
   sum_df <- sum_df[,c(1,4)]
   est <- round(exp(sum_df$Estimate),digits = 2)
   lower <- round(exp(CI[,1]), digits=2)
   upper<- round(exp(CI[,2]),digits=2)
   df[,1] <- paste0(est," (",lower,", ",upper,")")
   df[,2] <- round(sum_df[,2],digits=2)
   ind1 <- which(df[,2]<0.001,arr.ind = TRUE)
   df[ind1,2] <- "P<0.001"
   df <- df[c(1,6:8),]
   rownames(df) <- c("Interecpt","Time","COVID-19 period","Time: COVID-19 period")
   colnames(df) <- c("Estimate (95% CI range)", "p-value")
   return(df)
}

create_monthly_data <- function(df_merge,outcome){
   # Daily suicide data
   df <- df_merge %>% group_by(dt) %>%
      summarise(daily = sum(!!sym(outcome)),
                insured = sum(ninsured),
                proportion = daily / insured,
                percent = daily / insured * 100)
   # Monthly suicide data
   df_monthly <- df %>% group_by(year = lubridate::year(dt), month = lubridate::month(dt)) %>% 
      summarise(month_total = sum(daily),
                insured = mean(insured)) %>% 
      mutate(percent = month_total / insured * 100)
   # Add date column
   df_monthly$dt <- lubridate::ymd(paste(df_monthly$year, df_monthly$month, "01"))
   # Add time column (a sequence from 1 to 98)
   df_monthly$time <- seq(nrow(df_monthly))
   # Add the COVID-19 indicator
   df_monthly$indicator <- ifelse((df_monthly$year==2020 & df_monthly$month>2)| df_monthly$year==2021,1,0)
   return(df_monthly)
}

predictions <- function(model,df_monthly){
   # Compute predicted number of deaths from this model
   predA <- predict(model,type="response")
   df_monthly$pred <- predA
   # Compute counterfactual values
   new_df_monthly <- df_monthly
   new_df_monthly$indicator <- 0
   pred_counter <- predict(model,type="response",newdata = new_df_monthly)
   df_monthly$predC <- pred_counter
   df_monthly$predC[1:86] <- NA
   return(df_monthly)
}

RR <- function(model){
   # COVID-19 relative risk
   beta_COVID <- model$coefficients["indicator"]
   beta_interaction <- model$coefficients["time:indicator"]
   RR <- exp(beta_COVID+beta_interaction*87:98)
   meanRR <- mean(RR)
   geomRR <- prod(RR)^(1/12) #almost the same
   cov_mat <- vcov(model)
   b <- sum(87:98)/12
   var_lin <- cov_mat["indicator","indicator"]+b^2*cov_mat["time:indicator","time:indicator"]+b*2*cov_mat["indicator","time:indicator"]
   lin_est <- beta_COVID+b*beta_interaction
   CI_lin <- c(lin_est-1.96*sqrt(var_lin),lin_est+1.96*sqrt(var_lin))
   CI_RR <- exp(CI_lin)
   ret_vec <-c(geomRR,CI_RR)
   names(ret_vec) <- c("RR", "2.5%", "97.5%")
   return(ret_vec)
}

infer <- function(model,df_monthly,y_lab){
   # Extract regression estimates, including confidence intervals and p-values (eTable1)
   small_df <- est_df(model)
   # Add predictions
   df_monthly_new <- predictions(model,df_monthly)
   # Plot 
   p <- draw_plot(df_monthly_new,"month_total",y_lab) #+ ggtitle("Monthly Count of Suicide Attempts")
   ## Relative risk (RR) and confidence interval
   RR <- round(RR(model_main),digits=2)
return(list(small_df=small_df,df_monthly=df_monthly_new,p=p,RR=RR))
}