#probability of two consecutive FPs
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

#i=0: True negative
#i=1: False positive
#i=2: Competing event
#t: screening round
#h: censoring time (total number of rounds a woman participates in)
#v: number of prior false positives
#fm: baseline family history of breast cancer
#fd: baseline breast density
#rc: race\ethnicity
#int: interval (annual or biennial)
######################################################################
#probability of two consecutive FP based on censoring bias model

# P gives the p_ij ^(h,t,X) where 2<=t<h (X is a set of covariates)
P = function(t,h,v,age,fm,fd,rc,int){ 
  matrix(c(1-expit(-0.81-0.04*h-0.1*t+0.25*v-0.02*age-0.01*(fm==1)-0.44*(fd==1)+0.13*(fd==3)-0.01*(fd==4)+0.05*(rc==2)-0.32*(rc==3)-0.36*(rc==5)-0.21*(rc==8)-0.05*(rc==9)+0.06*(int==2)+0.005*h*t), 1-expit(-1.27+0.01*h-0.01*t+0.34*v-0.02*age-0.05*(fm==1)-0.33*(fd==1)+0.16*(fd==3)-0.05*(fd==4)+0.11*(rc==2)-0.04*(rc==3)+0.39*(rc==5)-0.15*(rc==8)+0.09*(rc==9)+0.07*(int==2)-0.01*h*t),0,
           expit(-0.81-0.04*h-0.1*t+0.25*v-0.02*age-0.01*(fm==1)-0.44*(fd==1)+0.13*(fd==3)-0.01*(fd==4)+0.05*(rc==2)-0.32*(rc==3)-0.36*(rc==5)-0.21*(rc==8)-0.05*(rc==9)+0.06*(int==2)+0.005*h*t), expit(-1.27+0.01*h-0.01*t+0.34*v-0.02*age-0.05*(fm==1)-0.33*(fd==1)+0.16*(fd==3)-0.05*(fd==4)+0.11*(rc==2)-0.04*(rc==3)+0.39*(rc==5)-0.15*(rc==8)+0.09*(rc==9)+0.07*(int==2)-0.01*h*t), 0,
           0, 0, 1), ncol=3)
}

#P_prime gives p_ij^(h,t,X) where 2<=t and t=h  (X is a set of covariates)
P_prime = function(t,h,v,age,fm,fd,rc,int){
  matrix(c(1/(1+exp(-1.22-0.04*h-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2)-0.001*h*t)+exp(-5.63-0.02*h-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2)-0.01*h*t)), 1/(1+exp(-1.49+0.03*h+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2)-0.01*h*t)+exp(-6.24-0.02*h-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2)-0.01*h*t)), 0,
           exp(-1.22-0.04*h-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2)-0.001*h*t)/(1+exp(-1.22-0.04*h-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2)-0.001*h*t)+exp(-5.63-0.02*h-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2)-0.01*h*t)),exp(-1.49+0.03*h+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2)-0.01*h*t)/(1+exp(-1.49+0.03*h+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2)-0.01*h*t)+exp(-6.24-0.02*h-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2)-0.01*h*t)), 0,
           exp(-5.63-0.02*h-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2)-0.01*h*t)/(1+exp(-1.22-0.04*h-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2)-0.001*h*t)+exp(-5.63-0.02*h-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2)-0.01*h*t)), exp(-6.24-0.02*h-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2)-0.01*h*t)/(1+exp(-1.49+0.03*h+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2)-0.01*h*t)+exp(-6.24-0.02*h-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2)-0.01*h*t)), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i,h,age,fm,fd,rc){
  if (h==1){
    if (i==0)
      return(1/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
    if(i==1)
      return(exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
    if (i==2)
      return(exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
  } else{
    if(i==1)
      return(expit(-1.97-0.03*h+0.02*age+0.09*(fm==1)-0.49*(fd==1)+0.17*(fd==3)-0.37*(fd==4)+0.01*(rc==2)-0.51*(rc==3)+0.11*(rc==5)-0.24*(rc==8)-0.18*(rc==9)))
    if(i==0)
      return(1-expit(-1.97-0.03*h+0.02*age+0.09*(fm==1)-0.49*(fd==1)+0.17*(fd==3)-0.37*(fd==4)+0.01*(rc==2)-0.51*(rc==3)+0.11*(rc==5)-0.24*(rc==8)-0.18*(rc==9)))
    if(i==2)
      return(0)
  }
}
#---------------------------------------------------------------------------------

MfpRisk_sim_real = function(Nsim, M, p, Nsubj,age,fm,fd,rc,int){
  if (int==2){M=M/2}
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  #simulating censoring time from Pr 
  out1 = vector(length=Nsim)
  for (jj in 1:Nsim){
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
    
    if(int==1){
      data1 = ddply(.data=data, .variables="StudyID_c", .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, "result"] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(round_1(0,1,age,fm,fd,rc),round_1(1,1,age,fm,fd,rc), round_1(2,1,age,fm,fd,rc))),
                                     sample(0:1, size=1, prob= c(round_1(0,df$S[1],age,fm,fd,rc),round_1(1,df$S[1],age,fm,fd,rc))))
          } else{
            v = sum(df[1:(i-1),"result"]==1)
            df[i, "result"] = ifelse(i< df$S[i], MC.sim(P(i,df$S[i],v,(age+(i-1)),fm,fd,rc,int), df$result[i-1]),
                                     MC.sim(P_prime(i,df$S[i],v,(age+(i-1)),fm,fd,rc,int), df$result[i-1]))
          }
        }
        return (df)
      }
      )
    }else{
      
      data1 = ddply(.data=data, .variables="StudyID_c", .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, "result"] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(round_1(0,1,age,fm,fd,rc),round_1(1,1,age,fm,fd,rc), round_1(2,1,age,fm,fd,rc))),
                                     sample(0:1, size=1, prob= c(round_1(0,df$S[1],age,fm,fd,rc),round_1(1,df$S[1],age,fm,fd,rc))))
          } else{
            v = sum(df[1:(i-1),"result"]==1)
            df[i, "result"] = ifelse(i< df$S[i], MC.sim(P(i,df$S[i],v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1]),
                                     MC.sim(P_prime(i,df$S[i],v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1]))
          }
        }
        return (df)
      }
      )
      
    }
    
    data2 = addCol(data1)
    
    data_new = ddply(.data = data2, .variable ="StudyID_c", .fun = function(df){
      x=which(df[,"result"]==1)
      y=x[-1] - x[-length(x)]
      cons=any(y==1)
      return(cons)
    })
    
    out1[jj]= table(data_new$V1)[2]/sum(table(data_new$V1))
  } 
  return(list(cumr = mean(out1), qtls = quantile(out1, probs = c(0,0.25, 0.75, 0.975, 1),  na.rm = TRUE)))
}

a1 = MfpRisk_sim_real(5000,10, 0.2, 1000, 40,1,3,2,1)#high risk, annual screener

a2 = MfpRisk_sim_real(5000,10, 0.2, 1000, 40,1,3,2,2)#high risk, biennial screener

a3 = MfpRisk_sim_real(5000,10, 0.2, 1000, 50,0,1,3,1) #low risk, annual screener

a4 =MfpRisk_sim_real(5000,10, 0.2, 1000, 50,0,1,3,2) #low risk, biennial screeners

####################################################################################

#probability of two consecutive FP based on population average model

# P gives the p_ij ^(h,t,X) where 2<=t<h (X is a set of covariates)
P = function(h,v,age,fm,fd,rc,int){ 
  matrix(c(1-expit(-0.81-0.04*h+0.25*v-0.02*age-0.01*(fm==1)-0.44*(fd==1)+0.13*(fd==3)-0.01*(fd==4)+0.05*(rc==2)-0.32*(rc==3)-0.36*(rc==5)-0.21*(rc==8)-0.05*(rc==9)+0.06*(int==2)), 1-expit(-1.27+0.01*h+0.34*v-0.02*age-0.05*(fm==1)-0.33*(fd==1)+0.16*(fd==3)-0.05*(fd==4)+0.11*(rc==2)-0.04*(rc==3)+0.39*(rc==5)-0.15*(rc==8)+0.09*(rc==9)+0.07*(int==2)),0,
           expit(-0.81-0.04*h+0.25*v-0.02*age-0.01*(fm==1)-0.44*(fd==1)+0.13*(fd==3)-0.01*(fd==4)+0.05*(rc==2)-0.32*(rc==3)-0.36*(rc==5)-0.21*(rc==8)-0.05*(rc==9)+0.06*(int==2)), expit(-1.27+0.01*h+0.34*v-0.02*age-0.05*(fm==1)-0.33*(fd==1)+0.16*(fd==3)-0.05*(fd==4)+0.11*(rc==2)-0.04*(rc==3)+0.39*(rc==5)-0.15*(rc==8)+0.09*(rc==9)+0.07*(int==2)), 0,
           0, 0, 1), ncol=3)
}

#P_prime gives p_ij^(h,t,X) where 2<=t and t=h  (X is a set of covariates) 
P_prime = function(h,v,age,fm,fd,rc,int){
  matrix(c(1/(1+exp(-1.22-0.04*h+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))+exp(-5.63-0.02*h+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))), 1/(1+exp(-1.49+0.03*h+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))+exp(-6.24-0.02*h+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))), 0,
           exp(-1.22-0.04*h+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))/(1+exp(-1.22-0.04*h+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))+exp(-5.63-0.02*h+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))),exp(-1.49+0.03*h+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))/(1+exp(-1.49+0.03*h+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))+exp(-6.24-0.02*h+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))), 0,
           exp(-5.63-0.02*h+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))/(1+exp(-1.22-0.04*h+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))+exp(-5.63-0.02*h+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))), exp(-6.24-0.02*h+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))/(1+exp(-1.49+0.03*h+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))+exp(-6.24-0.02*h+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i,h,age,fm,fd,rc){
  if (h==1){
    if (i==0)
      return(1/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
    if(i==1)
      return(exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
    if (i==2)
      return(exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
  } else{
    if(i==1)
      return(expit(-1.97-0.03*h+0.02*age+0.09*(fm==1)-0.49*(fd==1)+0.17*(fd==3)-0.37*(fd==4)+0.01*(rc==2)-0.51*(rc==3)+0.11*(rc==5)-0.24*(rc==8)-0.18*(rc==9)))
    if(i==0)
      return(1-expit(-1.97-0.03*h+0.02*age+0.09*(fm==1)-0.49*(fd==1)+0.17*(fd==3)-0.37*(fd==4)+0.01*(rc==2)-0.51*(rc==3)+0.11*(rc==5)-0.24*(rc==8)-0.18*(rc==9)))
    if(i==2)
      return(0)
  }
}
#-------------------------------------------------------------------------------------
MfpRisk_sim_real = function(Nsim, M, p, Nsubj,age,fm,fd,rc,int){
  if (int==2){M=M/2}
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  #simulating censoring time from Pr 
  out1 = vector(length=Nsim)
  for (jj in 1:Nsim){
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
    
    if(int==1){
      data1 = ddply(.data=data, .variables="StudyID_c", .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, "result"] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(round_1(0,1,age,fm,fd,rc),round_1(1,1,age,fm,fd,rc), round_1(2,1,age,fm,fd,rc))),
                                     sample(0:1, size=1, prob= c(round_1(0,df$S[1],age,fm,fd,rc),round_1(1,df$S[1],age,fm,fd,rc))))
          } else{
            v = sum(df[1:(i-1),"result"]==1)
            df[i, "result"] = ifelse(i< df$S[i], MC.sim(P(df$S[i],v,(age+(i-1)),fm,fd,rc,int), df$result[i-1]),
                                     MC.sim(P_prime(df$S[i],v,(age+(i-1)),fm,fd,rc,int), df$result[i-1]))
          }
        }
        return (df)
      }
      )
    }else{
      
      data1 = ddply(.data=data, .variables="StudyID_c", .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, "result"] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(round_1(0,1,age,fm,fd,rc),round_1(1,1,age,fm,fd,rc), round_1(2,1,age,fm,fd,rc))),
                                     sample(0:1, size=1, prob= c(round_1(0,df$S[1],age,fm,fd,rc),round_1(1,df$S[1],age,fm,fd,rc))))
          } else{
            v = sum(df[1:(i-1),"result"]==1)
            df[i, "result"] = ifelse(i< df$S[i], MC.sim(P(df$S[i],v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1]),
                                     MC.sim(P_prime(df$S[i],v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1]))
          }
        }
        return (df)
      }
      )
      
    }
    
    data2 = addCol(data1)
    
    data_new = ddply(.data = data2, .variable ="StudyID_c", .fun = function(df){
      x=which(df[,"result"]==1)
      y=x[-1] - x[-length(x)]
      cons=any(y==1)
      return(cons)
    })
    
    out1[jj]= table(data_new$V1)[2]/sum(table(data_new$V1))
  } 
  return(list(cumr = mean(out1),qtls = quantile(out1, probs = c(0,0.25, 0.75, 0.975, 1),  na.rm = TRUE)))
}

b1 = MfpRisk_sim_real(5000,10, 0.2, 50000, 40,1,3,2,1)#high risk, annual screener

b2 = MfpRisk_sim_real(500,10, 0.2, 50000, 40,1,3,2,2)#high risk, biennial screener

b3 = MfpRisk_sim_real(500,10, 0.2, 50000, 50,0,1,3,1) #low risk, annual screener

b4 =MfpRisk_sim_real(500,10, 0.2, 50000, 50,0,1,3,2) #low risk, biennial screeners

##############################################################################################
#probability of two consecutive FP based on discrete survival model

# P gives the p_ij ^(h,t,X) where 2<=t<=h (X is a set of covariates)
P = function(t,v,age,fm,fd,rc,int){ 
  matrix(c(1-expit(-0.81-0.1*t+0.25*v-0.02*age-0.01*(fm==1)-0.44*(fd==1)+0.13*(fd==3)-0.01*(fd==4)+0.05*(rc==2)-0.32*(rc==3)-0.36*(rc==5)-0.21*(rc==8)-0.05*(rc==9)+0.06*(int==2)), 1-expit(-1.27-0.01*t+0.34*v-0.02*age-0.05*(fm==1)-0.33*(fd==1)+0.16*(fd==3)-0.05*(fd==4)+0.11*(rc==2)-0.04*(rc==3)+0.39*(rc==5)-0.15*(rc==8)+0.09*(rc==9)+0.07*(int==2)),0,
           expit(-0.81-0.1*t+0.25*v-0.02*age-0.01*(fm==1)-0.44*(fd==1)+0.13*(fd==3)-0.01*(fd==4)+0.05*(rc==2)-0.32*(rc==3)-0.36*(rc==5)-0.21*(rc==8)-0.05*(rc==9)+0.06*(int==2)), expit(-1.27-0.01*t+0.34*v-0.02*age-0.05*(fm==1)-0.33*(fd==1)+0.16*(fd==3)-0.05*(fd==4)+0.11*(rc==2)-0.04*(rc==3)+0.39*(rc==5)-0.15*(rc==8)+0.09*(rc==9)+0.07*(int==2)), 0,
           0, 0, 1), ncol=3)
}

#P_prime gives p_ij^(h,t,X) where 2<=t and t=h  (X is a set of covariates)#not using this for discrete survival model
P_prime = function(t,v,age,fm,fd,rc,int){
  matrix(c(1/(1+exp(-1.22-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))+exp(-5.63-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))), 1/(1+exp(-1.49+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))+exp(-6.24-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))), 0,
           exp(-1.22-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))/(1+exp(-1.22-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))+exp(-5.63-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))),exp(-1.49+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))/(1+exp(-1.49+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))+exp(-6.24-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))), 0,
           exp(-5.63-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))/(1+exp(-1.22-0.04*t+0.23*v-0.01*age+0.08*(fm==1)-0.41*(fd==1)+0.14*(fd==3)+0.07*(fd==4)+0.12*(rc==2)-0.23*(rc==3)-0.61*(rc==5)-0.07*(rc==8)-0.11*(rc==9)-0.004*(int==2))+exp(-5.63-0.02*t+0.19*v+0.04*age+0.49*(fm==1)-0.15*(fd==1)-0.11*(fd==3)-0.09*(fd==4)+0.10*(rc==2)-0.30*(rc==3)+0.18*(rc==5)-0.51*(rc==8)-0.43*(rc==9)+0.13*(int==2))), exp(-6.24-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))/(1+exp(-1.49+0.03*t+0.24*v-0.01*age+0.06*(fm==1)-0.31*(fd==1)+0.15*(fd==3)+0.16*(fd==4)+0.15*(rc==2)-0.30*(rc==3)+0.04*(rc==5)-0.05*(rc==8)+0.19*(rc==9)+0.04*(int==2))+exp(-6.24-0.02*t+0.2*v+0.05*age-0.06*(fm==1)-0.44*(fd==1)+0.10*(fd==3)+0.41*(fd==4)-0.58*(rc==2)-0.60*(rc==3)+1.36*(rc==5)-1.27*(rc==8)-0.28*(rc==9)+0.33*(int==2))), 1), ncol=3)
}

#This function gives the probability of true negative (i=0),
#false positive (i=1) and competing event (i=2) when t=1 (first round)
#with censoring time (h) that is either 1 or greater than 1
round_1 = function(i,h,age,fm,fd,rc){
  if (h==1){
    if (i==0)
      return(1/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
    if(i==1)
      return(exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
    if (i==2)
      return(exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))/(1+exp(-2.18+0.02*age+0.20*(fm==1)-0.67*(fd==1)+0.18*(fd==3)-0.37*(fd==4)+0.26*(rc==2)-0.41*(rc==3)-0.01*(rc==5)-0.04*(rc==8)-0.17*(rc==9))+exp(-7.44+0.09*age+0.19*(fm==1)+0.25*(fd==1)+0.02*(fd==3)+0.002*(fd==4)-0.03*(rc==2)-0.63*(rc==3)+0.70*(rc==5)-0.76*(rc==8)-0.63*(rc==9))))
  } else{
    if(i==1)
      return(expit(-1.97+0.02*age+0.09*(fm==1)-0.49*(fd==1)+0.17*(fd==3)-0.37*(fd==4)+0.01*(rc==2)-0.51*(rc==3)+0.11*(rc==5)-0.24*(rc==8)-0.18*(rc==9)))
    if(i==0)
      return(1-expit(-1.97+0.02*age+0.09*(fm==1)-0.49*(fd==1)+0.17*(fd==3)-0.37*(fd==4)+0.01*(rc==2)-0.51*(rc==3)+0.11*(rc==5)-0.24*(rc==8)-0.18*(rc==9)))
    if(i==2)
      return(0)
  }
}
#---------------------------------------------------------------------------------

MfpRisk_sim_real = function(Nsim, M, p, Nsubj,age,fm,fd,rc,int){
  if (int==2){M=M/2}
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  #simulating censoring time from Pr 
  out1 = vector(length=Nsim)
  for (jj in 1:Nsim){
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
    
    if(int==1){
      data1 = ddply(.data=data, .variables="StudyID_c", .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, "result"] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(round_1(0,1,age,fm,fd,rc),round_1(1,1,age,fm,fd,rc), round_1(2,1,age,fm,fd,rc))),
                                     sample(0:1, size=1, prob= c(round_1(0,df$S[1],age,fm,fd,rc),round_1(1,df$S[1],age,fm,fd,rc))))
          } else{
            v = sum(df[1:(i-1),"result"]==1)
            #df[i, "result"] = ifelse(i< df$S[i], MC.sim(P(i,v,(age+(i-1)),fm,fd,rc,int), df$result[i-1]),
                                     #MC.sim(P(i,v,(age+(i-1)),fm,fd,rc,int), df$result[i-1]))
            
            df[i, "result"] = MC.sim(P(i,v,(age+(i-1)),fm,fd,rc,int), df$result[i-1])
                                      
          }
        }
        return (df)
      }
      )
    }else{
      
      data1 = ddply(.data=data, .variables="StudyID_c", .fun= function(df){
        df = df[order(df$round),]
        df$result = NA
        for (i in 1:nrow(df)){
          if (i == 1){
            df[1, "result"] = ifelse(df$S[1] == 1, sample(0:2, size=1, prob=c(round_1(0,1,age,fm,fd,rc),round_1(1,1,age,fm,fd,rc), round_1(2,1,age,fm,fd,rc))),
                                     sample(0:1, size=1, prob= c(round_1(0,df$S[1],age,fm,fd,rc),round_1(1,df$S[1],age,fm,fd,rc))))
          } else{
            v = sum(df[1:(i-1),"result"]==1)
            #df[i, "result"] = ifelse(i< df$S[i], MC.sim(P(i,v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1]),
                                     #MC.sim(P(i,v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1]))
            
            df[i, "result"] =  MC.sim(P(i,v,(2*i+(age-2)),fm,fd,rc,int), df$result[i-1])
                                     
          }
        }
        return (df)
      }
      )
      
    }
    
    data2 = addCol(data1)
  
    data_new = ddply(.data = data2, .variable ="StudyID_c", .fun = function(df){
      x=which(df[,"result"]==1)
      y=x[-1] - x[-length(x)]
      cons=any(y==1)
      return(cons)
    })
    
    out1[jj]= table(data_new$V1)[2]/sum(table(data_new$V1))
  } 
  return(list(cumr = mean(out1), qtls = quantile(out1, probs = c(0,0.25, 0.75, 0.975, 1),  na.rm = TRUE)))
}

c1 = MfpRisk_sim_real(5000,10, 0.2, 50000, 40,1,3,2,1)#high risk, annual screener

c2 = MfpRisk_sim_real(5000,10, 0.2, 50000, 40,1,3,2,2)#high risk, biennial screener

c3 = MfpRisk_sim_real(5000,10, 0.2, 50000, 50,0,1,3,1) #low risk, annual screener

c4 =MfpRisk_sim_real(5000,10, 0.2, 50000, 50,0,1,3,2) #low risk, biennial screeners

