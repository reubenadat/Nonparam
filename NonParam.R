#########################################################
# Job: NonParametric survival Method
# Name: Reuben Adatorwovor
# Date: September 04 2020
#
# #########################################################
rm(list = ls())
#library('VineCopula')
require('copula')
require(survival)
#library(tikzDevice)
library(latex2exp)
library(numDeriv)
library(boot)
#require("kmconfband")
#parameters##
 theta = 4/3 # the dependence parameter
 alpha1 = 5
 alpha2 = 2
 lambda1 = 1.2
 lambda2 = 2.1
 ###################################################
 #####
 direct = "/Users/adat/Documents/UNC/Dissertation/ProjectF/Rcode"
source ("/Users/adat/Documents/UNC/Dissertation/Source_func/Smart_func.r")
   set.seed(223)
   n.obs = 2500
   num_replicates = 1000
    gen_datagum000_100 = list()
   for (jj in seq(num_replicates)) {
   #clayton copula
   gumbel_cop = gumbelCopula(param = theta, dim = 2)
   bvd = mvdc(copula = gumbel_cop, margins = c("weibull","weibull"),paramMargins = list(list( shape = alpha1, scale = lambda1),list( shape = alpha2, scale = lambda2)) )
    all_data = smart_df(rMvdc(n.obs,bvd) )
    all_data = smart_df(cbind(T = pmin(all_data$X1, all_data$X2), C = runif(n.obs, 0, 11.7) ) )
    all_data$E = ifelse(all_data$T < all_data$C, 1, 0) 
    ### observed time to event is the min of (T, C)
	data_curr = transform(all_data, X = pmin(T, C) )
   	gen_datagum000_100[[jj]] = data_curr
   }  
   

   #  all_data$T = pmin(all_data$X1, all_data$X2)
    # all_data$C = runif(n.obs, 0, 7.1)
  ########### observed time to event is the min of (T, C)#############
     #    data_curr = transform(all_data, X = pmin(T, C), XX1 = pmin(X1, C), XX2 = pmin(X2, C) )
      #   data_curr$E = (ifelse(data_curr$T < data_curr$C, 1, 0)  ) #1==yes, 0==no
       #  data_curr$E1 = (ifelse(data_curr$XX1 < data_curr$C, 1, 0)  ) #1==yes, 0==no
# data_curr$E2 = (ifelse(data_curr$XX2 < data_curr$C, 1, 0)  ) #1==yes, 0==no
 
 #        data_curr = data_curr[,c("X","XX1","XX2","E","E1","E2")] #remove X1 and X2
  #   gen_bootgumbel000_2500[[jj]] = data_curr
   # datat1[[jj]] = as.numeric(table(data_curr$E)[2]/dim(data_curr)[1])  
   
   
   
   ##############Create Censoring Variable #####################################
#my_survival = function(currdata){
#test example code:
#what_rep = 1
#curr_data = gen_datagum000_100[[what_rep]][seq(sample_size),]
#data_cu = smart_df(curr_data)
#data_cur  = cbind(T = pmin(data_cu$X1, data_cu$X2), C = runif(n.obs, 0, 11.7) )
### observed time to event is the min of (T, C)
#data_curr = transform(data_cur, X = pmin(T, C) )
#data_curr$E = (ifelse(data_curr$T > data_curr$C, 0, 1)  )#1==yes, 0==no
##################Create true estimates###################
################################Save Data###################################
saveRDS(gen_datagum000_100,file = file.path(direct,"gen_datagum000_100.rds"))
 ################################Call data###################################
gen_datagum000_100 = readRDS(file.path(direct,"gen_datagum000_100.rds"))
#############################################################################
#############################################################################

func = function(dat,qq){
         df = dat[c("X","E")]
         #tt = quantile(gen_datagum000_100[[rr]]$X,qq)
         tt = quantile(dat$X,qq)
         #tt = quantile(dat$X)[c(2,3,4)] 
         
         S2 = with(df, Surv(X,E))
         surv_fit1 = survfit(formula = S2 ~ 1, data = df, conf.type = "log-log", type = "kaplan-meier")

         surv = stepfun(surv_fit1$time, c(1, surv_fit1$surv),ties = "ordered")
         #std.err = stepfun(surv_fit1$time, c(1, surv_fit1$std.err),ties = "ordered")
         #var_est = (std.err(tt))^2
         
         weib2 = 1 - pweibull(tt, alpha2, scale = lambda2, lower.tail = TRUE, log.p = FALSE)
         surv_est = exp(-( (-log(surv(tt)))^(theta) - (-log(weib2))^(theta) )^(1/theta) )
         as.numeric(surv_est)
  }

###########################################################################################
 func_stat = function(dat,qq) {
         len_qq = length(qq)
         s1_hat = rep(NA,len_qq)
         for(jj in seq(len_qq)){
                 s1_hat[jj] = func(dat = dat,qq = qq[jj])
                 #s1_hat[jj] = exp(-( (-log((survest(jj) )))^theta - ( -log(weib2[jj]))^theta  )^(1/theta) )
         }
         s1_hat
  }
#########################################################################################
 wrap_func_stat = function(dat,index,qq){
         func_stat(dat[index,],qq = qq)
  } 
#########################################################################################
 true_s1 = function(dat,qq){
         tt = quantile(dat$X,qq)
         weib1 = 1 - pweibull(tt, alpha1, scale = lambda1, lower.tail = TRUE, log.p = FALSE)
     as.numeric(weib1)
 }
################ initialize the parameters############
sample_size = N = n.obs;  qq = c(0.25,0.50,0.75)
all_estimates = c() ;  unc_params = c()

 #############################################################################


for (rr in sample_size ){
	## rr = 1
	cat(paste0("N = ",sample_size,"\n"))
   		for (what_rep in seq(num_replicates)) {
   			if (what_rep %% 5 == 0) cat(paste0(what_rep," "))
   			tt = quantile(gen_datagum000_100[[what_rep]]$X,qq)
   			propC = 1 - mean(gen_datagum000_100[[what_rep]]$E)
   			one_rep = boot(data = gen_datagum000_100[[what_rep]], statistic = wrap_func_stat, qq =qq, R = 500)
   			one_s1_hat =  one_rep$t0
   			one_var = diag(var(one_rep$t))
   			tru_s1 = true_s1(gen_datagum000_100[[what_rep]], qq )
   			unc_params = c(tt, one_s1_hat,tru_s1, one_var, propC)
         all_estimates = rbind(all_estimates, smart_df(N, t1 = unc_params[1], t2 = unc_params[2],t3 = unc_params[3], s1_hat.t1 = unc_params[4], s1_hat.t2 = unc_params[5], s1_hat.t3 = unc_params[6], s1.t1 = unc_params[7], s1.t2 = unc_params[8], s1.t3 = unc_params[9],var.t1 = unc_params[10],var.t2 = unc_params[11],var.t3 = unc_params[12], propC = unc_params[13] ) )
       
   			}	
   	  all_estimates
}
#######################################################################
 summary = c()
   # for (N in sample_size) {
          # tmp_dist = summary(all_estimates[which(all_estimates$N== N  ),])
 tmp_df = smart_df(t = apply(all_estimates[which(all_estimates$N == N), c("t1","t2","t3")],2, mean) )
 tmp_df$s1_hat = apply(all_estimates[which(all_estimates$N == N), c("s1_hat.t1","s1_hat.t2","s1_hat.t3")],2,mean)
 tmp_df$s1 = apply(all_estimates[which(all_estimates$N == N),c("s1.t1","s1.t2","s1.t3")],2,mean)
 tmp_df$bias = (tmp_df$s1_hat - tmp_df$s1)
 tmp_df$mod_var = apply(all_estimates[which(all_estimates$N == N),c("var.t1","var.t2","var.t3")],2,median)
 tmp_df$emp_var = apply(all_estimates[which(all_estimates$N == N), c("s1_hat.t1","s1_hat.t2","s1_hat.t3")],2, var)
 summary = rbind(summary,tmp_df)
#  }
 summary
 
###############################################################
 CP_s1.t1 = sum(all_estimates$s1_hat.t1 - qnorm(0.975) * sqrt(all_estimates$var.t1) < tmp_df$s1[1] & tmp_df$s1[1] < all_estimates$s1_hat.t1 + qnorm(0.975) * sqrt(all_estimates$var.t1) )/dim(all_estimates)[1]
 
 CP_s1.t2 = sum(all_estimates$s1_hat.t2 - qnorm(0.975) * sqrt(all_estimates$var.t2) < tmp_df$s1[2] & tmp_df$s1[2] < all_estimates$s1_hat.t2 + qnorm(0.975) * sqrt(all_estimates$var.t2) )/dim(all_estimates)[1]
 
 CP_s1.t3 = sum(all_estimates$s1_hat.t3 - qnorm(0.975) * sqrt(all_estimates$var.t3) <  tmp_df$s1[3] & tmp_df$s1[3] < all_estimates$s1_hat.t3 + qnorm(0.975) * sqrt(all_estimates$var.t3) )/dim(all_estimates)[1]
 
 CP = rbind(CP_s1.t1,CP_s1.t2,CP_s1.t3)
 
 summary_all = cbind(N,summary, CP)
 
 summary_all[,c(2:8)] = lapply(summary_all[c(2:8)],round,5)
 
  #apply(summary_all[,c(1:14)],1,function(x) cat(paste0(paste(x ,collapse = " & "))," \\\\\n"))

 apply(summary_all[,c(1:8)],1,function(x) cat(paste0(paste(x ,collapse = " & "))," \\\\\n"))
summary(all_estimates["propC"])













 surv(dat$X)
  qq = c(0.25,0.50,0.75)
  func(dat,qq)
  head(func(dat))
 
 ###################################################################
 
 func_stat = function(dat,qq) {
         # s1_hat =  sapply(xx, function(x) uniroot(func, interval= c(0, 1) )$root)
         #dat = all_reps1; qq = 0.25
         # dat = all_reps1; qq = c(0.25,.5)
         len_qq = length(qq)
         s1_hat = rep(NA,len_qq)
         for(jj in seq(len_qq)){
         #	jj =1
                 s1_hat[jj] = exp(-( (-log(surv(dat$X)))^(theta) - (-log(dat$weib2))^(theta) )^(1/theta) )
    dat = na.omit(dat)

                 #uniroot(func,dat = dat,qq = qq[jj], interval= c(0, 1.4) )$root
                 #s1_hat[jj] = exp(-( (-log((survest(jj) )))^theta - ( -log(weib2[jj]))^theta  )^(1/theta) )
         }
         s1_hat
 }
 
 
  func_stat(dat,qq)
 
 
 
 
 data_curr = func(data_curr)

 data_curr$weib1 = true_est(data_curr$X)
 ############Variance of the joint survival################
 surv_var = function(t){
	surv_err =  stepfun(surv_fit1$time, c(1, surv_fit1$std.err))
	surv_err(t)*surv_err(t)
}
###########################################################
func_deriv = function(theta, dat){
	
}
###########################################################


 func = function(dat,pars){
         df = data_curr[c("X","E")]
         tt = data_curr$X#quantile(dat$X,qq)
         #tt = quantile(dat$X)[c(2,3,4)] 
         S2 = with(df, Surv(X,E))
         surv_fit1 = survfit(formula = S2 ~ 1, data = df, conf.type = "log-log", type ="kaplan-meier")
         
  		}
  		
  		
  		
 surv_var( 0.0704 )
 0.0704    992       1  0.99899 0.00101     0.992866
 
 #func = function(dat){
 #	df = data_curr[c("X","E")]
  #       tt = data_curr$X#quantile(dat$X,qq)
         #tt = quantile(dat$X)[c(2,3,4)] 
   #      S2 = with(df, Surv(X,E))
   #      surv_fit1 = survfit(formula = S2 ~ 1, data = df, conf.type = "log-log", type ="kaplan-meier")
    #     data_curr$surv =  stepfun(surv_fit1$time, c(1, surv_fit1$surv))
    #     data_curr$weib2 = 1 - pweibull(tt, alpha2, scale = lambda2, lower.tail = TRUE, log.p = FALSE)
    #     data_curr$out_est = exp(-( (-log(surv(tt)))^(theta) - (-log(weib2))^(theta) )^(1/theta) )
 #}
 
 func = function(dat,qq){
         df = dat[c("X","E")]
         tt = quantile(dat$X,qq)
         #tt = quantile(dat$X)[c(2,3,4)] 
         S2 = with(df, Surv(X,E))
         surv_fit1 = survfit(formula = S2 ~ 1, data = df, conf.type = "log-log", type ="kaplan-meier")
         surv =  stepfun(surv_fit1$time, c(1, surv_fit1$surv))
         weib2 = exp(-(tt/lambda2)^alpha2)
         #surv_tt = surv(tt)
         out_est = exp(-( (-log(surv(tt)))^(theta) - (-log(weib2))^(theta) )^(1/theta))
         as.numeric(out_est)
}
 ########################################################
 func_stat = function(dat,qq) {
        len_qq = length(qq)
         s1_hat = rep(NA,len_qq)
         for(jj in seq(len_qq)){
                 s1_hat[jj] = quantile(out_est,qq[jj])
                          }
         s1_hat
 }
 ########################################################
 #s1_hat = func_stat(dat = all_reps1, qq = c(0.25,0.50,.75) )
 wrap_func_stat = function(dat,index,qq){
         func_stat(dat[index,],qq = qq)
 }
 ####################################################################
 
 ###################################################################
 true_s1 = function(dat,qq){
         tt = quantile(dat$X,qq)
         weib1 = exp(-(tt/lambda1)^alpha1)
     as.numeric(weib1)
 }
 

	
############################################# #################################

#data_curr$est1 = survest(data_curr$X) 
#bias = data_curr$est1
par(mfrow =c(2,2))
plot.ci(surv_fit1) 
 plot(data_curr[order(data_curr$X),][,c("X","weib1")], type = "l",col ="red",ylab = 'Survival Probability', xlab = 'Event Times', main ='Estimated Survival Probability from 1000 Samples  with 0 Dependence')
 lines(data_curr[order(data_curr$X),][,c("X","s1_hat")],type="l", col ='blue')
  lines(data_curr[order(data_curr$X),][,c("X","weib2")],type="l", col ='green')
    lines(data_curr[order(data_curr$X),][,c("X","weib1")],type="l", col ='orange')
 legend('bottomleft', legend = TeX(sprintf( c("$\\S_{T_1}(X)$", "$\\widehat{S_{T_1}(X)}$") )), col = c("red","blue"), lwd=1)
#surv_hat(q_time)
#true_est(q_time)
summary(data_curr)
plot(data_curr$X, data_curr$weib2)
#summary(surv_fit1$time)

####################### plot the two curves#################
 plot(surv_fit1,lty = 1, col="red",conf.int = FALSE)
 lines(data_curr[order(data_curr$X),][,c("X","weib2")],type="l", col ='blue')
 abline(h = true_est(q_time), col = "green",lty = 2 )
 #abline(v = c(q_time), col=c("yellow"), lty=c(2), lwd=c(2))
 #plot(data_curr[order(data_curr$X),][,c("X","weib1")],type="s", main = "Joint Gumbel Survival Probability",ylab = "Probability")

