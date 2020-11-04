#t : screening round
#h: censoring time
#v: number of prior false positives
#i=0: true negative
#i=1: false positive
#i=2: competing event

library("plyr")
library("dplyr")
library("nnet")
library("ggplot2")

#expit function
expit = function(x){
  return(exp(x)/(1+exp(x)))
}

# function for simulating one state from a transition matrix P with initial state x1
MC.sim = function(P,x1) {
  sim = as.numeric(1)
  m = ncol(P)
  sim = sample(0:(m-1),1, prob = P[(x1+1),])
  sim 
}

#add the variable previous result (s1) created from result
addCol = function(df){
  n = nrow(df)
  df = df[order(df$StudyID_c,df$round),]
  df$s1 = c(NA,df$result[-n])
  df$s1[match(unique(df$StudyID_c),df$StudyID_c )] = NA 
  return (df)
}



#this matrix gives p_ij^(h,t) where 2<=t and t=h 
#m40 is for row 1 ( when s1=0) and m41 is for row 2 (s1=1) (both multinomial regression with covariates h and t)
P_cb = function(m40, m41){
  function(h, t){
    matrix(c(1/(1+exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)+exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)), 1/(1+exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)+exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)), 0,
             exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)/(1+exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)+exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)),exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)/(1+exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)+exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)), 0,
             exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)/(1+exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)+exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)), exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)/(1+exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)+exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)), 1), ncol=3)
  }
}

#this matrix gives p_ij^(t) for all t
##m40_ds is multinomial regression coeffiecients for s1=0 (row 1) and m41_ds is multinomial 
#regression coeffiecient for s1=1 (row2) (regression on round)
P_ds = function(m40_ds, m41_ds){
  function(t){
    matrix(c(1/(1+exp(m40_ds[1]+m40_ds[3]*t)+exp(m40_ds[2]+m40_ds[4]*t)), 1/(1+exp(m41_ds[1]+m41_ds[3]*t)+exp(m41_ds[2]+m41_ds[4]*t)), 0,
             exp(m40_ds[1]+m40_ds[3]*t)/(1+exp(m40_ds[1]+m40_ds[3]*t)+exp(m40_ds[2]+m40_ds[4]*t)),exp(m41_ds[1]+m41_ds[3]*t)/(1+exp(m41_ds[1]+m41_ds[3]*t)+exp(m41_ds[2]+m41_ds[4]*t)), 0,
             exp(m40_ds[2]+m40_ds[4]*t)/(1+exp(m40_ds[1]+m40_ds[3]*t)+exp(m40_ds[2]+m40_ds[4]*t)), exp(m41_ds[2]+m41_ds[4]*t)/(1+exp(m41_ds[1]+m41_ds[3]*t)+exp(m41_ds[2]+m41_ds[4]*t)), 1), ncol=3)
  }
}



#this matrix gives p_ij^(h) for all h
#m40_pa is multinomial regression coeffiecients for s1=0 (row 1) and m41_pa is multinomial 
#regression coeffiecients for s1=1 (row2) (regression on censoring time h)
P_pa = function(m40_pa, m41_pa){
  function(h){
    matrix(c(1/(1+exp(m40_pa[1]+m40_pa[3]*h)+exp(m40_pa[2]+m40_pa[4]*h)), 1/(1+exp(m41_pa[1]+m41_pa[3]*h)+exp(m41_pa[2]+m41_pa[4]*h)), 0,
             exp(m40_pa[1]+m40_pa[3]*h)/(1+exp(m40_pa[1]+m40_pa[3]*h)+exp(m40_pa[2]+m40_pa[4]*h)),exp(m41_pa[1]+m41_pa[3]*h)/(1+exp(m41_pa[1]+m41_pa[3]*h)+exp(m41_pa[2]+m41_pa[4]*h)), 0,
             exp(m40_pa[2]+m40_pa[4]*h)/(1+exp(m40_pa[1]+m40_pa[3]*h)+exp(m40_pa[2]+m40_pa[4]*h)), exp(m41_pa[2]+m41_pa[4]*h)/(1+exp(m41_pa[1]+m41_pa[3]*h)+exp(m41_pa[2]+m41_pa[4]*h)), 1), ncol=3)
  }
}

n_ds = function(data, s1, s2, t){
  nrow(subset(data, s1==s1 & result==s2 & round==t))
}

n_pa  = function(data, s1, s2, h){
  nrow(subset(data, s1==s1 & result==s2 & S==h))
}

n_cb = function(data, s1, s2, t, h){
  nrow(subset(data, s1==s1 & result==s2 & round==t & S==h))
}

################################################################
#                         scenario 1:                          #
#   false positive risk is independent of censoring time,      #
#   number of prior false positive results and is constant     #
#                  across screening rounds                     #
################################################################

#using data generated from the first scenario to calculate AIC and compare the goodness of fits
#P is the transition matrix when 2<=t<h 
P = matrix(c(1-expit(-2.34), 1-expit(-1.86),0,
             expit(-2.34), expit(-1.86), 0,
             0, 0, 1), ncol=3)


#P_prime is the transition matrix when 2<=t and t=h or t>h  
P_prime = matrix(c(1/(1+exp(-2.18)+exp(-4.73)), 1/(1+exp(-1.65)+exp(-4.13)), 0,
                   exp(-2.18)/(1+exp(-2.18)+exp(-4.73)),exp(-1.65)/(1+exp(-1.65)+exp(-4.13)), 0,
                   exp(-4.73)/(1+exp(-2.18)+exp(-4.73)), exp(-4.13)/(1+exp(-1.65)+exp(-4.13)),1), ncol=3)

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==1)
      return(0.23)
    if(i==0)
      return(0.77)
    if(i==2)
      return(0)
  }
}

################################################

################################################################
#                         scenario 2:                          #
#   false positive risk depends on screening round but is      # 
#   independent of censoring time and number of prior          #
#             false positive results                           #
#            Case I : Strong dependency                        #
################################################################



#P is the transition matrix when 2<=t<h 
P = function(t){
  matrix(c(1-expit(-2.12-0.12*t), 1-expit(-1.75-0.08*t),0,
           expit(-2.12-0.12*t), expit(-1.75-0.08*t), 0,
           0, 0, 1), ncol=3)
}

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(t){
  matrix(c(1/(1+exp(-1.97-0.12*t)+exp(-4.78+0.02*t)), 1/(1+exp(-1.55-0.06*t)+exp(-4.33+0.12*t)), 0,
           exp(-1.97-0.12*t)/(1+exp(-1.97-0.12*t)+exp(-4.78+0.02*t)),exp(-1.55-0.06*t)/(1+exp(-1.55-0.06*t)+exp(-4.33+0.12*t)), 0,
           exp(-4.78+0.02*t)/(1+exp(-1.97-0.12*t)+exp(-4.78+0.02*t)), exp(-4.33+0.12*t)/(1+exp(-1.55-0.06*t)+exp(-4.33+0.12*t)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==1)
      return(0.23)
    if(i==0)
      return(0.77)
    if(i==2)
      return(0)
  }
}


################################################################
#                         scenario 2:                          #
#   false positive risk depends on screening round but is      # 
#   independent of censoring time and number of prior          #
#             false positive results                           #
#            Case II : Moderate dependency                      #
################################################################

#P is the transition matrix when 2<=t<h 
P = function(t){
  matrix(c(1-expit(-2.12-0.06*t), 1-expit(-1.75-0.04*t),0,
           expit(-2.12-0.06*t), expit(-1.75-0.04*t), 0,
           0, 0, 1), ncol=3)
}

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(t){
  matrix(c(1/(1+exp(-1.97-0.06*t)+exp(-4.78+0.01*t)), 1/(1+exp(-1.55-0.03*t)+exp(-4.33+0.06*t)), 0,
           exp(-1.97-0.06*t)/(1+exp(-1.97-0.06*t)+exp(-4.78+0.01*t)),exp(-1.55-0.03*t)/(1+exp(-1.55-0.03*t)+exp(-4.33+0.06*t)), 0,
           exp(-4.78+0.01*t)/(1+exp(-1.97-0.06*t)+exp(-4.78+0.01*t)), exp(-4.33+0.06*t)/(1+exp(-1.55-0.03*t)+exp(-4.33+0.06*t)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==1)
      return(0.23)
    if(i==0)
      return(0.77)
    if(i==2)
      return(0)
  }
}
##################################################################
################################################################
#                         scenario 3:                          #
#   false positive risk depends on the number of prior false   # 
#   positive results and is independent of screening round     # 
#              and censoring time.                             #
#            Case I : Strong dependency                        #
################################################################

#P is the transition matrix when 2<=t<h 
P = function(v){
  matrix(c(1-expit(-2.37+0.2*v), 1-expit(-2.07+0.02*v),0,
           expit(-2.37+0.2*v), expit(-2.07+0.02*v), 0,
           0, 0, 1), ncol=3)
}

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(v){
  matrix(c(1/(1+exp(-2.21+0.06*v)+exp(-4.8+0.16*v)), 1/(1+exp(-1.87+0.34*v)+exp(-4.55+0.68*v)), 0,
           exp(-2.21+0.06*v)/(1+exp(-2.21+0.06*v)+exp(-4.8+0.16*v)),exp(-1.87+0.34*v)/(1+exp(-1.87+0.34*v)+exp(-4.55+0.68*v)), 0,
           exp(-4.8+0.16*v)/(1+exp(-2.21+0.06*v)+exp(-4.8+0.16*v)), exp(-4.55+0.68*v)/(1+exp(-1.87+0.34*v)+exp(-4.55+0.68*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==1)
      return(0.23)
    if(i==0)
      return(0.77)
    if(i==2)
      return(0)
  }
}


################################################################
#                         scenario 3:                          #
#   false positive risk depends on the number of prior false   # 
#   positive results and is independent of screening round     # 
#              and censoring time.                             #
#            Case II : Moderate dependency                     #
################################################################

#P is the transition matrix when 2<=t<h 
P = function(v){
  matrix(c(1-expit(-2.37+0.1*v), 1-expit(-2.07+0.01*v),0,
           expit(-2.37+0.1*v), expit(-2.07+0.01*v), 0,
           0, 0, 1), ncol=3)
}

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(v){
  matrix(c(1/(1+exp(-2.21+0.03*v)+exp(-4.8+0.08*v)), 1/(1+exp(-1.87+0.17*v)+exp(-4.55+0.34*v)), 0,
           exp(-2.21+0.03*v)/(1+exp(-2.21+0.03*v)+exp(-4.8+0.08*v)),exp(-1.87+0.17*v)/(1+exp(-1.87+0.17*v)+exp(-4.55+0.34*v)), 0,
           exp(-4.8+0.08*v)/(1+exp(-2.21+0.03*v)+exp(-4.8+0.08*v)), exp(-4.55+0.34*v)/(1+exp(-1.87+0.17*v)+exp(-4.55+0.34*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==1)
      return(0.23)
    if(i==0)
      return(0.77)
    if(i==2)
      return(0)
  }
}

##################################################################

################################################################
#                         scenario 4:                          #
#    false positive risk depends on the number of prior        #
#    false positive results and screening rounds but is        #
#         independent of censoring time                        #   
#           Case I: Strong dependency                          #
################################################################

#P is the transition matrix when 2<=t<h 
P = function(t, v){
  matrix(c(1-expit(-2.07-0.09*t+0.27*v), 1-expit(-1.97-0.1*t+0.32*v),0,
           expit(-2.07-0.09*t+0.27*v), expit(-1.97-0.1*t+0.32*v), 0,
           0, 0, 1), ncol=3)
}

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(t, v){
  matrix(c(1/(1+exp(-1.92-0.09*t+0.27*v)+exp(-4.72-0.02*t+0.32*v)), 1/(1+exp(-1.78-0.09*t+0.33*v)+exp(-4.56+0.01*t+0.33*v)), 0,
           exp(-1.92-0.09*t+0.27*v)/(1+exp(-1.92-0.09*t+0.27*v)+exp(-4.72-0.02*t+0.32*v)),exp(-1.78-0.09*t+0.33*v)/(1+exp(-1.78-0.09*t+0.33*v)+exp(-4.56+0.01*t+0.33*v)), 0,
           exp(-4.72-0.02*t+0.32*v)/(1+exp(-1.92-0.09*t+0.27*v)+exp(-4.72-0.02*t+0.32*v)), exp(-4.56+0.01*t+0.33*v)/(1+exp(-1.78-0.09*t+0.33*v)+exp(-4.56+0.01*t+0.33*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==0)
      return(0.77)
    if(i==1)
      return(0.23)
    if(i==2)
      return(0)
  }
}


################################################################
#                         scenario 4:                          #
#    false positive risk depends on the number of prior        #
#    false positive results and screening rounds but is        #
#         independent of censoring time                        #   
#           Case II: Moderate dependency                       #
################################################################

#P is the transition matrix when 2<=t<h 
P = function(t, v){
  matrix(c(1-expit(-2.07-0.045*t+0.135*v), 1-expit(-1.97-0.05*t+0.16*v),0,
           expit(-2.07-0.045*t+0.135*v), expit(-1.97-0.05*t+0.16*v), 0,
           0, 0, 1), ncol=3)
}

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(t, v){
  matrix(c(1/(1+exp(-1.92-0.045*t+0.135*v)+exp(-4.72-0.01*t+0.16*v)), 1/(1+exp(-1.78-0.045*t+0.165*v)+exp(-4.56+0.005*t+0.165*v)), 0,
           exp(-1.92-0.045*t+0.135*v)/(1+exp(-1.92-0.045*t+0.135*v)+exp(-4.72-0.01*t+0.16*v)),exp(-1.78-0.045*t+0.165*v)/(1+exp(-1.78-0.045*t+0.165*v)+exp(-4.56+0.005*t+0.165*v)), 0,
           exp(-4.72-0.01*t+0.16*v)/(1+exp(-1.92-0.045*t+0.135*v)+exp(-4.72-0.01*t+0.16*v)), exp(-4.56+0.005*t+0.165*v)/(1+exp(-1.78-0.045*t+0.165*v)+exp(-4.56+0.005*t+0.165*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i, h){
  if (h==1){
    if (i==0)
      return(0.74)
    if(i==1)
      return(0.25)
    if (i==2)
      return(0.01)
  } else{
    if(i==0)
      return(0.77)
    if(i==1)
      return(0.23)
    if(i==2)
      return(0)
  }
}

###################

#M is the maximum number of screening rounds
M = 10
# censoring time is simulated using geometric distribution with parameter 0.2
p=0.2
P_s = p*((1-p)^seq(0,(M-1)))
#sum of geometric distribution for censoring time should be 1
Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))



#this function is used for simulation in scenario 1,2,3 and 4
Goodness_of_fit = function(Nsim, M, p, Nsubj, P, P_prime, scen){
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  AIC = matrix (ncol=3, nrow= Nsim) 
  for(jj in 1:Nsim){
    #simulating censoring time from Pr ( one vector of length M for each subject)
    data.s = t(rmultinom(Nsubj,1,Pr))
    #find the round of censoring
    sup = data.s %*% seq(1,M)
    Su = ifelse(sup>M, M, sup)
    #Now Su is the censoring time so we repeat it for M times 
    S = rep(Su, each = M)  # Now S is the censoring time
    # creating round variable 
    round = rep(seq(1,M), Nsubj)
    #creating StudyID_c varisble
    StudyID_c = rep(seq(1,Nsubj), each = M)
    #Now create initial data frame
    data  = data.frame (S = S, round = round, StudyID_c = StudyID_c)
    #add the variable result (0,1,2) to the data using the transition matrix and probabilities
    if (scen==1){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(0.74, 0.25, 0.01)),
                                     sample(0:1, size=1, prob= c(0.77,0.23)))
          } else{
            
            df[i, 'result'] = ifelse(i< df$S[i], MC.sim(P, df$result[i-1]),
                                     MC.sim(P_prime, df$result[i-1]))
          }
        }
        return (df)
      }
      )  
    }
    
    if (scen==2){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(0.74, 0.25, 0.01)),
                                     sample(0:1, size=1, prob= c(0.77,0.23)))
          } else{
            df[i, 'result'] = ifelse(i< df$S[i], MC.sim(P(i), df$result[i-1]),
                                     MC.sim(P_prime(i), df$result[i-1]))
          }
        }
        return (df)
      }
      )
    }
    
    if (scen==3){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(0.74, 0.25, 0.01)),
                                     sample(0:1, size=1, prob= c(0.77,0.23)))
          } else{
            v = sum(df[1:(i-1),'result']==1)
            df[i, 'result'] = ifelse(i< df$S[i], MC.sim(P(v), df$result[i-1]),
                                     MC.sim(P_prime(v), df$result[i-1]))
          }
        }
        return (df)
      }
      )
    }
    
    if(scen==4){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(0.74, 0.25, 0.01)),
                                     sample(0:1, size=1, prob= c(0.77, 0.23)))
          } else{
            v = sum(df[1:(i-1),'result']==1)
            df[i, 'result'] = ifelse(i< df$S[i], MC.sim(P(i,v), df$result[i-1]),
                                     MC.sim(P_prime(i,v), df$result[i-1]))
          }
        }
        return (df)
      }
      )
    }
    
    data2 = addCol(data1)
    p1 = summary(factor(data2$S[!duplicated(data2$StudyID_c)],levels=c(1:M)))/length(data2$S[!duplicated(data2$StudyID_c)]) 
    
    set0 = subset(data2, s1==0)
    set1 = subset(data2, s1==1)
    
    mod_ds0 = multinom(result~round , data = set0)
    m_ds0 = summary(mod_ds0)$coefficients 
    m_ds0[is.na(m_ds0)]=0
    mod_ds1 = multinom(result~ round , data = set1)
    m_ds1 = summary(mod_ds1)$coefficients 
    m_ds1[is.na(m_ds1)]=0
    
    mod_pa0 = multinom(result~S , data = set0)
    m_pa0 = summary(mod_pa0)$coefficients 
    m_pa0[is.na(m_pa0)]=0
    mod_pa1 = glm(result~ S , data = set1)
    m_pa1 = summary(mod_pa1)$coefficients 
    m_pa1[is.na(m_pa1)]=0
    
    
    mod_cb0 = multinom(result~S+round+S*round , data = set0)
    m_cb0 = summary(mod_cb0)$coefficients 
    m_cb0[is.na(m_cb0)]=0
    mod_cb1 = multinom(result~S+round+S*round , data = set1)
    m_cb1 = summary(mod_cb1)$coefficients 
    m_cb1[is.na(m_cb1)]=0
    
    
    
    f1 = P_ds(m_ds0, m_ds1)#function of t
    f2 = P_pa(m_pa0, m_pa1)# function of S
    f3 = P_cb(m_cb0, m_cb1)# function of t and S
    
    A1 = expand.grid(0:2, 0:2, 1:10, 1:10)
    colnames(A1) = c("s1", "s2", "h", "j")
    #A1 = A1[A1$j <= A1$h,]
    L1 = apply(A1, 1, FUN = function(x){(f1(x[4])[x[1]+1, x[2]+1])^n_ds(data2, x[1], x[2], x[4])})
    
    
    A2 = expand.grid(0:2, 0:2, 1:10)
    colnames(A2) = c("s1", "s2", "h")
    L2 =  apply(A2, 1, FUN = function(x){((f2(x[3])[x[1]+1, x[2]+1])*p1[x[3]])^n_pa(data2, x[1], x[2], x[3])})
    
    A3 = expand.grid(0:2, 0:2, 1:10,1:10,1:10)
    colnames(A3) = c("s1", "s2", "h", "j", "l")
    A3 = A3[A3$j <= A3$h,]
    L3 = apply(A3, 1, FUN = function(x){((p1[x[3]])^n_cb(data2,x[1], x[2],x[4], x[3]))*(f3(x[3], x[5])[x[1]+1, x[2]+1]^n_cb(data2,x[1], x[2],x[5], x[3]))})
    
    Lh1 = subset(L1, L1!=0)     
    Lh2 = subset(L2, L2!=0) 
    Lh3 = subset(L3, L3!=0) 
    
    AIC1 = 2*2-2*sum(log(Lh1))
    AIC2 = 2*2-2*sum(log(Lh2))
    AIC3 = 2*4-2*sum(log(Lh3))
    
    AIC[jj, 1] = AIC1
    AIC[jj, 2] = AIC2
    AIC[jj, 3] = AIC3
  }
  
  return(list(AIC = apply(AIC, 2, mean),  estr = apply(AIC, 2, sd)))
}



Goodness_of_fit (10, 10, 0.2, 50000, P, P_prime )


######################################################


################################################################
#                         scenario 5:                          #
#    false positive risk depends on the number of prior        #
#  false positive results, censoring time and screening round  #   
#           Case I: Strong dependency                          #
################################################################

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(t,v){
  matrix(c(1/(1+exp(-2.02-0.09*t+0.27*v)+exp(-5.85+0.01*t+0.32*v)), 1/(1+exp(-1.89-0.09*t+0.32*v)+exp(-5.53+0.02*t+0.3*v)), 0,
           exp(-2.02-0.09*t+0.27*v)/(1+exp(-2.02-0.09*t+0.27*v)+exp(-5.85+0.01*t+0.32*v)),exp(-1.89-0.09*t+0.32*v)/(1+exp(-1.89-0.09*t+0.32*v)+exp(-5.53+0.02*t+0.3*v)), 0,
           exp(-5.85+0.01*t+0.32*v)/(1+exp(-2.02-0.09*t+0.27*v)+exp(-5.85+0.01*t+0.32*v)), exp(-5.53+0.02*t+0.3*v)/(1+exp(-1.89-0.09*t+0.32*v)+exp(-5.53+0.02*t+0.3*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
}


################################################################
#                         scenario 5:                          #
#    false positive risk depends on the number of prior        #
#  false positive results, censoring time and screening round  #   
#           Case II: Moderate dependency                       #
################################################################


#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(t,v){
  matrix(c(1/(1+exp(-2.02-0.045*t+0.135*v)+exp(-5.85+0.005*t+0.16*v)), 1/(1+exp(-1.89-0.045*t+0.16*v)+exp(-5.53+0.01*t+0.15*v)), 0,
           exp(-2.02-0.045*t+0.135*v)/(1+exp(-2.02-0.045*t+0.135*v)+exp(-5.85+0.005*t+0.16*v)),exp(-1.89-0.045*t+0.16*v)/(1+exp(-1.89-0.045*t+0.16*v)+exp(-5.53+0.01*t+0.15*v)), 0,
           exp(-5.85+0.005*t+0.16*v)/(1+exp(-2.02-0.045*t+0.135*v)+exp(-5.85+0.005*t+0.16*v)), exp(-5.53+0.01*t+0.15*v)/(1+exp(-1.89-0.045*t+0.16*v)+exp(-5.53+0.01*t+0.15*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
  
}

################
################################################################
#                         scenario 6:                          #
#    false positive risk depends on the number of prior        # 
#      false positive results and censoring time but is        #
#        constant across screening                             #   
#           Case I: Strong dependency                          #
################################################################


#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(v){
  matrix(c(1/(1+exp(-2.31+0.14*v)+exp(-5.81+0.33*v)), 1/(1+exp(-1.98+0.17*v)+exp(-5.50+0.34*v)), 0,
           exp(-2.31+0.14*v)/(1+exp(-2.31+0.14*v)+exp(-5.81+0.33*v)),exp(-1.98+0.17*v)/(1+exp(-1.98+0.17*v)+exp(-5.50+0.34*v)), 0,
           exp(-5.81+0.33*v)/(1+exp(-2.31+0.14*v)+exp(-5.81+0.33*v)), exp(-5.50+0.34*v)/(1+exp(-1.98+0.17*v)+exp(-5.50+0.34*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
}


################################################################
#                         scenario 6:                          #
#    false positive risk depends on the number of prior        # 
#      false positive results and censoring time but is        #
#        constant across screening                             #   
#           Case II: Moderate dependency                       #
################################################################


#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function(v){
  matrix(c(1/(1+exp(-2.31+0.07*v)+exp(-5.81+0.165*v)), 1/(1+exp(-1.98+0.085*v)+exp(-5.50+0.17*v)), 0,
           exp(-2.31+0.07*v)/(1+exp(-2.31+0.07*v)+exp(-5.81+0.165*v)),exp(-1.98+0.085*v)/(1+exp(-1.98+0.085*v)+exp(-5.50+0.17*v)), 0,
           exp(-5.81+0.165*v)/(1+exp(-2.31+0.07*v)+exp(-5.81+0.165*v)), exp(-5.50+0.17*v)/(1+exp(-1.98+0.085*v)+exp(-5.50+0.17*v)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
}


####################


################################################################
#                         scenario 7:                          #
#   false positive risk depends on censoring time              #
#   and is independent of screening round and number           #
#          of prior false positive results                     #   
#         Strong and moderate dependency (depends on alpha)    #
################################################################

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = matrix(c(1/(1+exp(-2.28)+exp(-5.72)), 1/(1+exp(-1.77)+exp(-5.08)), 0,
                   exp(-2.28)/(1+exp(-2.28)+exp(-5.72)),exp(-1.77)/(1+exp(-1.77)+exp(-5.08)), 0,
                   exp(-5.72)/(1+exp(-2.28)+exp(-5.72)), exp(-5.08)/(1+exp(-1.77)+exp(-5.08)), 1), ncol=3)


#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
  
}

######################################


################################################################
#                         scenario 8:                          #
#         false positive risk depends on censoring             #
#     time and screening round but is independent              #
#       of number of prior false positive results              #   
#           Case I: Strong dependency                          #
################################################################

#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function( t){
  matrix(c(1/(1+exp(-2.07-0.06*t)+exp(-5.9+0.05*t)), 1/(1+exp(-1.67-0.03*t)+exp(-5.32+0.08*t)), 0,
           exp(-2.07-0.06*t)/(1+exp(-2.07-0.06*t)+exp(-5.9+0.05*t)),exp(-1.67-0.03*t)/(1+exp(-1.67-0.03*t)+exp(-5.32+0.08*t)), 0,
           exp(-5.9+0.05*t)/(1+exp(-2.07-0.06*t)+exp(-5.9+0.05*t)), exp(-5.32+0.08*t)/(1+exp(-1.67-0.03*t)+exp(-5.32+0.08*t)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
  
}



#P_prime is the transition matrix when 2<=t and t=h or t>h 
P_prime = function( t){
  matrix(c(1/(1+exp(-2.07-0.03*t)+exp(-5.9+0.025*t)), 1/(1+exp(-1.67-0.015*t)+exp(-5.32+0.04*t)), 0,
           exp(-2.07-0.03*t)/(1+exp(-2.07-0.03*t)+exp(-5.9+0.025*t)),exp(-1.67-0.015*t)/(1+exp(-1.67-0.015*t)+exp(-5.32+0.04*t)), 0,
           exp(-5.9+0.025*t)/(1+exp(-2.07-0.03*t)+exp(-5.9+0.025*t)), exp(-5.32+0.04*t)/(1+exp(-1.67-0.015*t)+exp(-5.32+0.04*t)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) for the round 1 
round_1 = function(i){
  if (i==0)
    return(0.75)
  if(i==1)
    return(0.24)
  if (i==2)
    return(0.01)
  
}

#this function is used for simulation in scenario 1,2,3 and 4
Goodness_of_fit1 = function(Nsim, M, p, Nsubj, P, P_prime, scen, d="strong"){
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  AIC = matrix (ncol=3, nrow= Nsim) 
  for(jj in 1:Nsim){
    #Now create initial data frame
    round = rep(seq(1,M), Nsubj)
    #creating StudyID_c varisble
    StudyID_c = rep(seq(1,Nsubj), each = M)
    #Now create initial data frame
    data  = data.frame (round = round, StudyID_c = StudyID_c)
    #add the variable result (0,1,2) to the data using the transition matrix and probabilities
    if (scen==5){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = sample(0:2, size=1, prob=c(0.75, 0.24, 0.01))
            
          } else{
            v = sum(df[1:(i-1),'result']==1)
            df[i, 'result'] =  MC.sim(P_prime(i,v), df$result[i-1])
            
          }
        }
        return (df)
      }
      )
    }
    if (scen==6){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = sample(0:2, size=1, prob=c(0.75, 0.24, 0.01))
            
          } else{
            v = sum(df[1:(i-1),'result']==1)
            df[i, 'result'] =  MC.sim(P_prime(v), df$result[i-1])
            
          }
        }
        return (df)
      }
      ) 
    }
    if (scen==7){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = sample(0:2, size=1, prob=c(0.75, 0.24, 0.01))
            
          } else{
            
            df[i, 'result'] =  MC.sim(P_prime, df$result[i-1])
            
          }
        }
        return (df)
      }
      )
    }
    
    if (scen==8){
      data1 = ddply(.data=data, .variables='StudyID_c', .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, 'result'] = sample(0:2, size=1, prob=c(0.75, 0.24, 0.01))
            
          } else{
            
            df[i, 'result'] =  MC.sim(P_prime(i), df$result[i-1])
            
          }
        }
        return (df)
      }
      )
      
    }
    
    data2 = addCol(data1)
    #add the variable for the round of the first fp
    data2 = ddply(.data = data2, .variable ='StudyID_c', .fun = function(df){
      df$fp1 = ifelse(is.na(which(df$result==1)[1]), M+1, which(df$result==1)[1])
      return(df)
    })
    
    #add the variable for the round of second fp (this is the event time for at least 2 FP)
    data2 = ddply(.data = data2, .variable ='StudyID_c', .fun = function(df){
      df$fp2 = ifelse(is.na(which(df$result==1)[2]), M+1, which(df$result==1)[2])
      return(df)
    })
    
    fp2 = ddply(.data=data2, .variable ='StudyID_c', .fun = function(df){ 
      fp2 = unique(df$fp2)
      return(fp2)
    })
    
    #define the round of competing event (first 2)
    data2 = ddply(.data = data2, .variable ='StudyID_c', .fun = function(df){
      df$comp = ifelse(is.na(which(df$result==2)[1]), M+1, which(df$result==2)[1])
      return(df)
    })
    
    if(scen==5|scen==8){
      if(d=="strong"){
        alp=-0.1
      }else{
        alp=-0.075
      }
    }
    if(scen==6|scen==7){
      if(d=="strong"){
        alp=-0.075
      }else{
        alp=-0.05
      }
    }
    p1= do.call(rbind,lapply(fp2[,2],function(x){expit(alp*(x-(seq(1,M+1)+1)))}))
    cump = t(apply(1-cbind(0,p1[,-(M+1)]),1,cumprod))
    PS = t(apply(cbind(p1,cump),1,function(x){x[1:(M+1)]*x[(M+2):(2*(M+1))]}))
    Su = apply(PS,1,function(x){rmultinom(1,1,x)})
    Sup    = t(Su)%*%seq(1,M+1)
    S = ifelse(Sup > M, M, Sup)
    data2$S = rep(S,each = M)
    data2$S = ifelse(data2$S<=data2$comp, data2$S, data2$comp)
    p1 = summary(factor(data2$S[!duplicated(data2$StudyID_c)],levels=c(1:M)))/length(data2$S[!duplicated(data2$StudyID_c)]) 
    
    set0 = subset(data2, s1==0)
    set1 = subset(data2, s1==1)
    
    mod_ds0 = multinom(result~round , data = set0)
    m_ds0 = summary(mod_ds0)$coefficients 
    m_ds0[is.na(m_ds0)]=0
    mod_ds1 = multinom(result~ round , data = set1)
    m_ds1 = summary(mod_ds1)$coefficients 
    m_ds1[is.na(m_ds1)]=0
    
    mod_pa0 = multinom(result~S , data = set0)
    m_pa0 = summary(mod_pa0)$coefficients 
    m_pa0[is.na(m_pa0)]=0
    mod_pa1 = glm(result~ S , data = set1)
    m_pa1 = summary(mod_pa1)$coefficients 
    m_pa1[is.na(m_pa1)]=0
    
    
    mod_cb0 = multinom(result~S+round+S*round , data = set0)
    m_cb0 = summary(mod_cb0)$coefficients 
    m_cb0[is.na(m_cb0)]=0
    mod_cb1 = multinom(result~S+round+S*round , data = set1)
    m_cb1 = summary(mod_cb1)$coefficients 
    m_cb1[is.na(m_cb1)]=0
    
    
    
    f1 = P_ds(m_ds0, m_ds1)#function of t
    f2 = P_pa(m_pa0, m_pa1)# function of S
    f3 = P_cb(m_cb0, m_cb1)# function of t and S
    
    A1 = expand.grid(0:2, 0:2, 1:10, 1:10)
    colnames(A1) = c("s1", "s2", "h", "j")
    #A1 = A1[A1$j <= A1$h,]
    L1 = apply(A1, 1, FUN = function(x){(f1(x[4])[x[1]+1, x[2]+1])^n_ds(data2, x[1], x[2], x[4])})
    
    
    A2 = expand.grid(0:2, 0:2, 1:10)
    colnames(A2) = c("s1", "s2", "h")
    L2 =  apply(A2, 1, FUN = function(x){((f2(x[3])[x[1]+1, x[2]+1])*p1[x[3]])^n_pa(data2, x[1], x[2], x[3])})
    
    A3 = expand.grid(0:2, 0:2, 1:10,1:10,1:10)
    colnames(A3) = c("s1", "s2", "h", "j", "l")
    A3 = A3[A3$j <= A3$h,]
    L3 = apply(A3, 1, FUN = function(x){((p1[x[3]])^n_cb(data2,x[1], x[2],x[4], x[3]))*(f3(x[3], x[5])[x[1]+1, x[2]+1]^n_cb(data2,x[1], x[2],x[5], x[3]))})
    
    Lh1 = subset(L1, L1!=0)     
    Lh2 = subset(L2, L2!=0) 
    Lh3 = subset(L3, L3!=0) 
    
    AIC1 = 2*2-2*sum(log(Lh1))
    AIC2 = 2*2-2*sum(log(Lh2))
    AIC3 = 2*4-2*sum(log(Lh3))
    
    AIC[jj, 1] = AIC1
    AIC[jj, 2] = AIC2
    AIC[jj, 3] = AIC3
  }
  
  return(list(AIC = apply(AIC, 2, mean),  estr = apply(AIC, 2, sd)))
}



Goodness_of_fit (5000, 10, 0.2, 50000, P, P_prime )


