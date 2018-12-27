                                 #R code used for simulation (scenario 1-8)
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

#matrices needed for censoring bias model (lower triangle) and also for uncensored subjects
#in population average model, p_i^(h,1) in the manuscript
round_11 = function(m1, m2){
  function(i, h){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1[1])+exp(m1[2])))
      if(i==1)
        return(exp(m1[1])/(1+exp(m1[1])+exp(m1[2])))
      if (i==2)
        return(exp(m1[2])/(1+exp(m1[1])+exp(m1[2])))
    } else{
      if(i==1)
        return(expit(m2[1]+m2[2]*h))
      if(i==0)
        return(1 -expit(m2[1]+m2[2]*h))
      if(i==2)
        return(0)
    }
  }
}

# this matrix gives the p_ij ^(h,t) where 2<=t<h
#m30 is for row 1 (when s1=0), m31 is for row 2 (when s1=1) (both logistic regression with covarites h and t)
P2 = function(m30 , m31){
  function(h, t){
    matrix(c(1-expit(m30[1]+m30[2]*h+m30[3]*t+m30[4]*h*t), 1-expit(m31[1]+m31[2]*h+m31[3]*t+m31[4]*h*t),0,
             expit(m30[1]+m30[2]*h+m30[3]*t+m30[4]*h*t), expit(m31[1]+m31[2]*h+m31[3]*t+m31[4]*h*t), 0,
             0, 0, 1), ncol=3)}
}

#this matrix gives p_ij^(h,t) where 2<=t and t=h 
#m40 is for row 1 ( when s1=0) and m41 is for row 2 (s1=1) (both multinomial regression with covariates h and t)
P_prime2 = function(m40, m41){
  function(h, t){
    matrix(c(1/(1+exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)+exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)), 1/(1+exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)+exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)), 0,
             exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)/(1+exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)+exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)),exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)/(1+exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)+exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)), 0,
             exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)/(1+exp(m40[1]+m40[3]*h+m40[5]*t+m40[7]*h*t)+exp(m40[2]+m40[4]*h+m40[6]*t+m40[8]*h*t)), exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)/(1+exp(m41[1]+m41[3]*h+m41[5]*t+m41[7]*h*t)+exp(m41[2]+m41[4]*h+m41[6]*t+m41[8]*h*t)), 1), ncol=3)
  }
}


#matrices needed for censoring bias model (upper triangle)

#receiving result i (i can be 0,1,2) in the first round when 
# censoring time is h (h can be 1 or greater than 1) and n_g=x 
#where n_g is number of false positives by and including time g
round_1_nh = function(m1_nh, m2_nh){
  function(i, h, g, x){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)+exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)))
      if(i==1)
        return(exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)/(1+exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)+exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)))
      if (i==2)
        return(exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)/(1+exp(m1_nh[2*(g-1)+1,1]+m1_nh[2*(g-1)+1,2]*x)+exp(m1_nh[2*g,1]+m1_nh[2*g,2]*x)))
    } else{
      if(i==1)
        return(expit(m2_nh[g,1]+m2_nh[g,2]*h+m2_nh[g,3]*x))
      if(i==0)
        return(1-expit(m2_nh[g,1]+m2_nh[g,2]*h+m2_nh[g,3]*x))
      if(i==2)
        return(0)
    }
  }
}

# this matrix gives the p_ij ^(h,t,n_g=x) where 2<=t<h 
P_nh = function(m30_nh, m31_nh){
  function(h,t,g,x){
    matrix(c(1-expit(m30_nh[g,1]+m30_nh[g,2]*x+m30_nh[g,3]*h+m30_nh[g,4]*t+m30_nh[g,5]*h*t), 1-expit(m31_nh[g,1]+m31_nh[g,2]*x+m31_nh[g,3]*h+m31_nh[g,4]*t+m31_nh[g,5]*h*t),0,
             expit(m30_nh[g,1]+m30_nh[g,2]*x+m30_nh[g,3]*h+m30_nh[g,4]*t+m30_nh[g,5]*h*t), expit(m31_nh[g,1]+m31_nh[g,2]*x+m31_nh[g,3]*h+m31_nh[g,4]*t+m31_nh[g,5]*h*t), 0,
             0, 0, 1), ncol=3)}
}

#this matrix gives p_ij^(h,t,n_g=x) where 2<=t and t=h, for all i and j
P_prime_nh = function(m40_nh, m41_nh){
  function(h, t, g, x){
    matrix(c(1/(1+exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)+exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)), 1/(1+exp(m41_nh[2*(g-1)+1,1]+m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)+exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)), 0,
             exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)/(1+exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)+exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)),exp(m41_nh[2*(g-1)+1,1]+m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)/(1+exp(m41_nh[2*(g-1)+1,1]+m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)+exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)), 0,
             exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)/(1+exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)+exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)), exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)/(1+exp(m41_nh[2*(g-1)+1,1]+m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)+exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)), 1), ncol=3)
  }}


#analytical variance estimation for censoring bias model
cum.se = function(theta, M, G, N, N2, p){
  G = rep(G,M)
  for (j in (M-1):1){
    for (l in (j+1):M){
      for (s in 1:M){ 
        N[j,l] = N[j,l] + exp(G[j]*(s-(j+1)))*N2[l,s]
      }
    }		
  }
  
  var.theta = theta*(1-theta)
  cov.theta = vector(length = M)
  temp = NULL
  for (i in 1:M){
    temp = -matrix(theta[i,],ncol =1)%*%matrix(theta[i,],nrow =1)
    diag(temp) = 0
    cov.theta[i] = sum(c(temp))
  }	
  var.theta = apply(sweep(var.theta,1,p[1:M]^2/apply(N,1,sum),"*"),2,sum)-cov.theta/apply(N,1,sum)
  var.theta = sum(var.theta)
}

cens.bias_Markov = function(data, M, alpha, l, functions){
  xround_11 = functions$round_11
  xP_prime2 = functions$P_prime2
  xP2 = functions$P2
  xround_1_nh = functions$round_1_nh
  xP_nh = functions$P_nh
  xP_prime_nh = functions$P_prime_nh
  #data = subset(data, round <= cens2)
  theta = matrix(data = 0, nrow=M, ncol=M)
  # compute theta[h, j] for lower triangle
  cartesianProd = function(a){
    if(a==1){
      if(l==1){
        return(data.frame('i1'=1))
      } else{ stop('not possible')}
    } else{
      df = expand.grid(rep(list(c(0,1)), a-1))
      names(df) = paste0('i', 1:(a-1))
      df = subset(df, rowSums(df)==l-1)
      df[, paste0('i',a )] = 1
      return(df)}}
  ########################## j = 1, h= 1:M creates the first column
  for (j in l:M){
    df = cartesianProd(j)
    for (h in j:M){
      df$term = apply(df,1 , FUN = function(x){
        x = as.numeric(x)
        C = xround_11(x[1], h) #pi1
        if (j>1){
          C = C*prod(sapply(1:(j-1), function(k){
            if (k+1 == h){
              return (xP_prime2(h, k+1)[x[k]+1, x[k+1]+1])
            } else{
              return (xP2(h, k+1)[x[k]+1, x[k+1]+1])
            }
          }))}
        return (C)
      })
      theta[h, j] = sum(df$term, na.rm = T)
    }}
  #####################################################################
  B = list() # B[[h]][b,c] = P(S=b, nh = c-1)  
  A = list() # A[[h]][a,b,c] = P(w_l=a | S=b, nh = c-1)
  for(h in 1:(M-1)){
    B[[h]] = matrix(data=0, nrow=M, ncol=min(l-1, h)+1)  
    for (b in h:M)
      for (c in 0:min(l-1, h))
        B[[h]][b, c+1] = nrow(subset(data, S==b & get(paste0('n',h))==c))/nrow(data)
      
      A[[h]] = array(data=0, dim=c(M+1, M, 1+min(l-1, h)))
      # case I: a<=b
      for (b in max(l, h+1):M)
        for (a in l:b)
          for (c in 0:min(l-1, h)){
            df = cartesianProd(a)
            df$term = apply(df,1 , FUN = function(x){
              x = as.numeric(x)
              C = xround_1_nh(x[1], b, h, c)
              if (a>1){
                C = C*prod(sapply(1:(a-1), function(k){
                  if (k+1 == b){
                    return (xP_prime_nh(b, k+1, h, c)[x[k]+1, x[k+1]+1])
                  } else{
                    return (xP_nh(b, k+1, h, c)[x[k]+1, x[k+1]+1])
                  }  
                }))}
              return (ifelse(is.nan(C), 0, C))
            })
            A[[h]][a, b, c+1] = sum(df$term, na.rm = T)
          }
      # case II: a > b = M
       for (c in 0:min(l-1, h))
    A[[h]][M+1, M, c+1] = 1 - sum(sapply(l:M, function(jj) {A[[h]][jj, M, c+1]}))
      
      # case III: a >b & b < M. Also in the final round we calculate A[[h]][a,h,c+1] 
      for (b in (M-1):h)
        for (a in max(l,(b+1)):ifelse(b>h, M+1, M))
          for (c in 0:min(l-1, h)){
            if(a< M+1){
              E = diag(exp(alpha*(l-c)*seq(M+1)))
              A[[h]][a, b, c+1] = sum(A[[h]][a, , c+1]*B[[h]][, c+1])*exp(alpha*(l-c)*a)/sum(E%*%as.matrix(A[[h]][, , c+1]%*%B[[h]][, c+1, drop=F]))
            } else{
              A[[h]][a, b, c+1] = 1 - sum(sapply(l:M, function(jj) {A[[h]][jj, b, c+1]}))
            }
            
          }
  }
  
  # compute theta[h, j] for upper triangle
  for (h in 1:(M-1))
    for (j in max(h+1, l):M)
      theta[h, j] = sum(A[[h]][j,h, ]*B[[h]][h, ])*nrow(data)/nrow(subset(data, S==h))
  
  p1 = summary(factor(data$S[!duplicated(data$StudyID_c)],levels=c(1:M)))/length(data$S[!duplicated(data$StudyID_c)]) 
  D = apply(theta,1,sum)
  D = as.matrix(D, drop=F)
  risk = p1%*%D
  risk = as.numeric(risk)
  data$RESINIT_C = ifelse(data$fp2==data$round, 1,0)
  N1 = table(data$S[data$RESINIT_C==1], data$round[data$RESINIT_C==1])
  N2 = matrix(0, ncol=M, nrow=M)
  N2[1:nrow(N1), 1:ncol(N1)]=N1
  N <- table(data$S,data$round)
  varp = cum.se(theta, M, alpha, N, N2, p1)
  se = sqrt(varp)
  return (list(risk=risk, se = se))
}

#matrices needed for discrete survival model
#matrix for first round, m1_ds is when 1=t=h and m2_ds is when 1=t <h, this gives p_i^(1) 
#m1_ds is multinomial regression on 1 when t=h=1 and m2_ds is logistic regression on 1 when 1=t<h
round_11_ds = function(m1_ds, m2_ds){ 
  function(i, h){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1_ds[1])+exp(m1_ds[2])))
      if(i==1)
        return(exp(m1_ds[1])/(1+exp(m1_ds[1])+exp(m1_ds[2])))
      if (i==2)
        return(exp(m1_ds[2])/(1+exp(m1_ds[1])+exp(m1_ds[2])))
    } else{
      if(i==1)
        return(expit(m2_ds[1]))
      if(i==0)
        return(1-expit(m2_ds[1]))
      if(i==2)
        return(0)
    }
  }
}

# this matrix gives the p_ij ^(t) where 2<=t<h
#m30_ds is logistic regression coefficients for s1=0 (row 1) and m31_ds is 
#logistic regression coeffiecients for s1=1 (row 2) (regression on t)
P2_ds = function(m30_ds , m31_ds){ 
  function(t){
    matrix(c(1-expit(m30_ds[1]+m30_ds[2]*t), 1-expit(m31_ds[1]+m31_ds[2]*t),0,
             expit(m30_ds[1]+m30_ds[2]*t), expit(m31_ds[1]+m31_ds[2]*t), 0,
             0, 0, 1), ncol=3)}
}

#this matrix gives p_ij^(t) where 2<=t and t=h
##m40_ds is multinomial regression coeffiecients for s1=0 (row 1) and m41_ds is multinomial 
#regression coeffiecient for s1=1 (row2) (regression on round)
P_prime2_ds = function(m40_ds, m41_ds){
  function(t){
    matrix(c(1/(1+exp(m40_ds[1]+m40_ds[3]*t)+exp(m40_ds[2]+m40_ds[4]*t)), 1/(1+exp(m41_ds[1]+m41_ds[3]*t)+exp(m41_ds[2]+m41_ds[4]*t)), 0,
             exp(m40_ds[1]+m40_ds[3]*t)/(1+exp(m40_ds[1]+m40_ds[3]*t)+exp(m40_ds[2]+m40_ds[4]*t)),exp(m41_ds[1]+m41_ds[3]*t)/(1+exp(m41_ds[1]+m41_ds[3]*t)+exp(m41_ds[2]+m41_ds[4]*t)), 0,
             exp(m40_ds[2]+m40_ds[4]*t)/(1+exp(m40_ds[1]+m40_ds[3]*t)+exp(m40_ds[2]+m40_ds[4]*t)), exp(m41_ds[2]+m41_ds[4]*t)/(1+exp(m41_ds[1]+m41_ds[3]*t)+exp(m41_ds[2]+m41_ds[4]*t)), 1), ncol=3)
  }
}

#analytical variance estimation for discrete survival model
disc.surv.var = function(phat.gw, S, fp2.cens, comp, M){ #phat.gw is the estimated risk by discrete survival model
  r = vector (length = M)
  s = vector (length =  M)
  for (k in 1:M){ 
    r[k] =  length(fp2.cens[S>=k & (fp2.cens>=k | fp2.cens==0) ]) 
    s[k]= sum(fp2.cens[S>=k & (fp2.cens>=k | fp2.cens==0) ]>k|fp2.cens[S>=k & (fp2.cens>=k | fp2.cens==0) ]==0)
  }
  inner = r-s/r ; inner = inner/s
  varp = ((1-phat.gw)^2) * sum(inner)
  return(varp)
}

#function needed for defining discrete survival model in Markov setting
createCartesian = function(l, h){
  L = list()
  for (i in 1:h)
    L[[i]] = c(0,1,2)
  if (h ==1){
    df = data.frame('i1' = 1)
  } else{
    df = expand.grid(L)
    for (i in 1:h)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=l)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==h-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
  }
  return (df)
}

#discrete survival model in Markov setting
f = function(x, p1, functions0){
  x = as.numeric(x)
  round_11_ds = functions0$round_11_ds
  P2_ds = functions0$P2_ds
  P_prime2_ds = functions0$P_prime2_ds
  h = length(x)
  if(length(x)==1)
    return (round_11_ds(1,h)*p1[h])
  c = round_11_ds(x[1], h)*p1[h]
  for (k in 1:(length(x)-1)){
    if (k+1 == h){
      c = c*P_prime2_ds(k+1)[x[k]+1, x[k+1]+1]
    } else{
      c = c*P2_ds(k+1)[x[k]+1, x[k+1]+1]
    }
  }
  return (c)
}


disc.surv_Markov = function(data, functions0, l, M){ 
  p1 = summary(factor(data$S[!duplicated(data$StudyID_c)],levels=c(1:M)))/length(data$S[!duplicated(data$StudyID_c)]) 
  P = 0 # initial probability
  for (j in l:M){
    if(j==1){
      df = data.frame('i1'=1)
    } else{
      df = createCartesian(l, j)
    }
    df$term = apply(df,1, f, p1, functions0=functions0)
    P = P + sum(df$term)
  }
  varp = disc.surv.var(P, data$S, data$fp2.cens, data$comp, M)
  se = sqrt(varp)
  return (list(risk_ds=P, se_ds = se))
}

#discrete survival model 
disc.surv = function(data,M,l){
  theta_hat = vector(length = M)
  r = vector (length = M)
  s = vector (length =  M)
  varp = vector (length =M)
  for (k in 1:M){
    theta_hat[k] = mean(data$fp2.cens[data$S>=k & (data$fp2.cens>=k | data$fp2.cens==0)] == k)
    r[k] =  length(data$fp2.cens[data$S>=k & (data$fp2.cens>=k | data$fp2.cens==0) ]) 
    s[k]= sum(data$fp2.cens[data$S>=k & (data$fp2.cens>=k | data$fp2.cens==0) ]>k|data$fp2.cens[data$S>=k & (data$fp2.cens>=k | data$fp2.cens==0) ]==0)
  }
  theta_hat = ifelse(is.na(theta_hat),0, theta_hat)
  s = ifelse(is.na(s),0, s)
  r = ifelse(is.na(r),0, r)
  cum = sum(theta_hat*cumprod(1-c(0,theta_hat[1:M-1])))
  inner = (r-s)/r; inner = inner/s
  varp = ((1-cum)^2) * sum(inner)
  se = sqrt(varp)
  return (list(risk_ds=cum, se_ds = se))
}  

#matrices needed for population average model
#for uncensored subjects we use round_11, P_2 and P_prime2 that I already defined for censoring bias model
##m1_pa is multinomial regression on 1 when t=h=1 and m2_pa is logistic regression on h when 1=t <h
round_11_pa = function(m1_pa, m2_pa){ 
  function(i, h){
    if (h==1){
      if (i==0)
        return(1/(1+exp(m1_pa[1])+exp(m1_pa[2])))
      if(i==1)
        return(exp(m1_pa[1])/(1+exp(m1_pa[1])+exp(m1_pa[2])))
      if (i==2)
        return(exp(m1_pa[2])/(1+exp(m1_pa[1])+exp(m1_pa[2])))
    } else{
      if(i==1)
        return(expit(m2_pa[1]+m2_pa[2]*h))
      if(i==0)
        return(1-expit(m2_pa[1]+m2_pa[2]*h))
      if(i==2)
        return(0)
    }
  }
}

# this matrix gives the p_ij ^(h) where 2<=t<h
#m30_pa is logistic regression coefficients for s1=0 (row 1) and m31 is
#logistic regression coeffiecients for s1=1 (row 2) (regression on censoring time h)
P2_pa = function(m30_pa , m31_pa){ 
  function(h){
    matrix(c(1-expit(m30_pa[1]+m30_pa[2]*h), 1-expit(m31_pa[1]+m31_pa[2]*h),0,
             expit(m30_pa[1]+m30_pa[2]*h), expit(m31_pa[1]+m31_pa[2]*h), 0,
             0, 0, 1), ncol=3)}
}


#this matrix gives p_ij^(h) where 2<=t and t=h
#m40_pa is multinomial regression coeffiecients for s1=0 (row 1) and m41_pa is multinomial 
#regression coeffiecients for s1=1 (row2) (regression on censoring time h)
P_prime2_pa = function(m40_pa, m41_pa){
  function(h){
    matrix(c(1/(1+exp(m40_pa[1]+m40_pa[3]*h)+exp(m40_pa[2]+m40_pa[4]*h)), 1/(1+exp(m41_pa[1]+m41_pa[3]*h)+exp(m41_pa[2]+m41_pa[4]*h)), 0,
             exp(m40_pa[1]+m40_pa[3]*h)/(1+exp(m40_pa[1]+m40_pa[3]*h)+exp(m40_pa[2]+m40_pa[4]*h)),exp(m41_pa[1]+m41_pa[3]*h)/(1+exp(m41_pa[1]+m41_pa[3]*h)+exp(m41_pa[2]+m41_pa[4]*h)), 0,
             exp(m40_pa[2]+m40_pa[4]*h)/(1+exp(m40_pa[1]+m40_pa[3]*h)+exp(m40_pa[2]+m40_pa[4]*h)), exp(m41_pa[2]+m41_pa[4]*h)/(1+exp(m41_pa[1]+m41_pa[3]*h)+exp(m41_pa[2]+m41_pa[4]*h)), 1), ncol=3)
  }
}


#analytical variance estimation for population average model
pop.ave.var = function(S, fp2.cens, comp, delta, M){
  sig = vector(length = M)
  N = max(S)*(max(S)+3)/2-1
  # sigma_1
  eta = vector(length = N)
  k = 1
  for (i in 1:max(S)){
    for (j in i:max(S)){
      eta[k] = mean(S == j & fp2.cens == i & delta == 1 & comp>=i)
      k = k+1
    }	
  }
  for (i in 1:(max(S)-1)){
    eta[k+i-1] = mean(S == i & delta == 0) 
  }
  
  sigma = matrix(ncol = N, nrow = N)
  sigma = -outer(eta,eta)
  diag(sigma) = eta*(1-eta)
  
  sig[1] = t(c(rep(1,max(S)),rep(0,N-max(S))))%*%sigma%*%c(rep(1,max(S)),rep(0,N-max(S)))/length(S)
  
  # sigma_2 to sigma_K
  mu0 = vector(length = M)
  for (i in 1:M){
    mu0[i] = mean(S == i & delta == 0)
  }
  mu = vector("list", max(S))
  for (i in 1:max(S)){
    for (j in i:max(S)){
      mu[[i]][j-i+1] = mean(S == j & fp2.cens == i & delta == 1 & comp>=i)
    }	
  }
  
  for (i in 2:max(S)){
    k = 1
    a = rep(0,N)
    for (j in 1:i){
      a[k:(k+max(S)-j)]  = (1 + (1-i/j)*mu0[j]^(i/j)*(mu0[j]+rep(1,(max(S)-j+1))%*%mu[[j]])^(-i/j))%*%rep(1,(max(S)-j+1))
      k = k+max(S)-j+1
    }
    a[(length(a)-(max(S)-2)):length(a)] = c(rep(1,(i-1)),rep(0,(max(S)-i)))
    for (j in 1:(i-1)){
      a[(length(a)-(max(S)-2)+j-1)] = a[(length(a)-(max(S)-2)+j-1)] +
        (mu0[j]^(i/j-1)*(mu0[j]+sum(mu[[j]]))^(-i/j)*(mu0[j]+i/j*sum(mu[[j]])))
    }
    sig[i] = a%*%sigma%*%a/length(S)				
  }
  
  return(varp = sig[M])
}

#functions needed for defining population average model in Markov setting
createCartesian1 = function(l, j){
  L = list()
  for (i in 1:(j-1))
    L[[i]] = c(0,1)
  if (j ==1){
    df = data.frame('i1' = 1)
  } else{
    df = expand.grid(L)
    for (i in 1:(j-1))
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var==l-1)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==j-ind)})
    df = subset(df, valid==1)
    df$valid = NULL
    df[, paste0('i',j)] = 1}
  return (df)
}


f2 = function(x, h, functions2){
  round_11 = functions2$round_11
  P2 = functions2$P2
  P_prime2 = functions2$P_prime2
  x = as.numeric(x)
  c = round_11(x[1], h)*x[length(x)]
  for (k in 1:(length(x)-2)){
    if (k+1 == h){
      c = c*P_prime2(h,k+1)[x[k]+1, x[k+1]+1]
    } else{
      c = c*P2(h,k+1)[x[k]+1, x[k+1]+1]
    }
  }
  return (c)
}

f3 = function(x, h, functions1){
  round_11_pa = functions1$round_11_pa
  P2_pa = functions1$P2_pa
  P_prime2_pa = functions1$P_prime2_pa
  x = as.numeric(x)
  c = round_11_pa(x[1], h)*x[length(x)]
  for (k in 1:(length(x)-2)){
    if (k+1 >= h){
      c = c*P_prime2_pa(h)[x[k]+1, x[k+1]+1]
    } else{
      c = c*P2_pa(h)[x[k]+1, x[k+1]+1]
    }
  }
  return (c)
}

#population average model in Markov setting
pop.ave_Markov = function(data, functions1, functions2,l,M){
  p1 = summary(factor(data$S[!duplicated(data$StudyID_c)],levels=c(1:M)))/length(data$S[!duplicated(data$StudyID_c)]) 
  P = 0 # initial probability
  for (j in l:M)
    for (h in j:M){
      df = createCartesian1(l,j)
      df[, 'p1'] = p1[h]
      df$term = apply(df,1, f2, h=h, functions2=functions2)
      P = P + sum(df$term)
    }
  for (j in max(2,l):M)
    for (h in 1:(j-1)){
      df = createCartesian1(l,j)
      df[, 'p1'] = p1[h]
      df$term = apply(df,1, f3, h=h, functions1=functions1)
      P = P + sum(df$term)
    }
  varp = pop.ave.var(data$S, data$fp2.cens, data$comp, data$delta, M)
  se = sqrt(varp)
  return(list(risk_pa = P, se_pa = se))
}

#M is the maximum number of screening rounds
M = 10
# censoring time is simulated using geometric distribution with parameter 0.2
p=0.2
P_s = p*((1-p)^seq(0,(M-1)))
#sum of geometric distribution for censoring time should be 1
Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))

#this function is used for simulation in scenario 1,2,3 and 4
MfpRisk_sim1 = function(Nsim, M, p, l, Nsubj, P, P_prime, scen){
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  out1 = matrix (ncol=3, nrow= Nsim) 
  se1 = matrix (ncol=3, nrow= Nsim)
  for (jj in 1:Nsim){
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
    
    #add the variable for the round of the first fp
    data3 = ddply(.data = data2, .variable ='StudyID_c', .fun = function(df){
      df$fp1 = ifelse(is.na(which(df$result==1)[1]), M+1, which(df$result==1)[1])
      return(df)
    })
    
    #add the variable for the round of second fp (this is the event time for at least 2 FP)
    data4 = ddply(.data = data3, .variable ='StudyID_c', .fun = function(df){
      df$fp2 = ifelse(is.na(which(df$result==1)[2]), M+1, which(df$result==1)[2])
      return(df)
    })
    
    #add M variables ni when ni is the number of false positive by and including round i for i=1,...,M
    data5 = ddply(.data= data4, .variables = 'StudyID_c', .fun = function(df) {
      for (i in 1:M){
        df[, paste0("n",i)] =  sum(df$result[1:i]==1)
      }
      return(df)
    })
    
    #add the variable for the time of cancer (when result is 2)
    data6 = ddply(.data = data5, .variable ='StudyID_c', .fun = function(df){
      df$comp = ifelse(is.na(which(df$result==2)[1]), M+1, which(df$result==2)[1])
      return(df)
    })
    
    #add variable for the censoring if first fp is the event (min(S,fp1))
    data6$cens1 = apply(cbind(data6$S,data6$fp1),1,min)  
    
    #add variable for the censoring if second fp is the event (min(S,fp2))
    data6$cens2 = apply(cbind(data6$S,data6$fp2),1,min)  
    
    #define event time for the other two models
    data6$fp2.cens = ifelse(data6$fp2>data6$S,0,data6$fp2) 
    
    #define censorin indicator for the other two models
    data6$delta = ifelse(data6$fp2 <= data6$S, 1, 0)
    
    set1 = subset(data6,  1==round & round ==S)
    set2 = subset(data6,  1==round & round <= (S-1))
    set3 = subset(data6,  2<=round & round <= (S-1))
    set4 = subset(data6,  2<=round & round ==S)
    set30 = subset(set3, s1==0)
    set31 = subset(set3, s1==1)
    set40 = subset(set4, s1==0)
    set41 = subset(set4, s1==1)
    mod1 = multinom(result~1 , data = set1)
    m1 = summary(mod1)$coefficients 
    m1[is.na(m1)]=0
    mod2 = glm(result~ S , data = set2, family=binomial)
    m2 = mod2$coefficients 
    m2[is.na(m2)]=0
    mod30 = glm(result~ S+round+S*round , data = set30, family=binomial)
    m30 = mod30$coefficients 
    m30[is.na(m30)]=0
    mod31 = glm(result~ S+round+S*round , data = set31, family=binomial)
    m31 = mod31$coefficients 
    m31[is.na(m31)]=0
    mod40 = multinom(result~S+round+S*round , data = set40)
    m40 = summary(mod40)$coefficients 
    m40[is.na(m40)]=0
    mod41 = multinom(result~S+round+S*round , data = set41)
    m41 = summary(mod41)$coefficients 
    m41[is.na(m41)]=0
    m1_nh = matrix (0, nrow = 2*M, ncol=2) 
    for (i in 1:M){
      mod1 = multinom(result~ set1[,paste0("n",i)] , data = set1)
      m1_nh[2*(i-1)+1,] = summary(mod1)$coefficients[c(1,3)]
      m1_nh[2*i,] = summary(mod1)$coefficients[c(2,4)]
    }
    
    m1_nh[is.na(m1_nh)]=0
    
    m2_nh = matrix (0, nrow = M, ncol=3) 
    for (i in 1:M){
      mod2 = glm(result~ S+set2[,paste0("n",i)], data = set2, family=binomial)#gives the coefficient estimates of (9) 
      m2_nh[i,] = mod2$coefficients
    }
    
    m2_nh[is.na(m2_nh)]=0
    
    m30_nh = matrix (0, nrow = M, ncol=5)
    for (i in 1:M){
      mod30 = glm(result~set30[,paste0("n",i)]+S+round+S*round ,data = set30, family=binomial)#row 1 of matrix (7)
      m30_nh[i,] = mod30$coefficients
    }
    m30_nh[is.na(m30_nh)]=0
    
    m31_nh = matrix(0,nrow = M, ncol=5) 
    for (i in 1:M){
      mod31 = glm(result~set31[,paste0("n",i)]+S+round+S*round , data =set31 , family=binomial)#row 2 of matrix (7)
      m31_nh[i,] = mod31$coefficients
    }
    
    m31_nh[is.na(m31_nh)]=0
    
    m40_nh = matrix(0, nrow = 2*M, ncol = 5)
    for (i in 1:M){
      mod40 = multinom(result~set40[,paste0("n",i)]+S+round+S*round, data =set40 )   
      m40_nh[2*(i-1)+1,] = summary(mod40)$coefficients[c(1,3,5,7,9)]
      m40_nh[2*i,] = summary(mod40)$coefficients[c(2,4,6,8,10)]
    }
    
    m40_nh[is.na(m40_nh)]=0
    
    m41_nh = matrix(0,nrow = 2*M, ncol = 5) 
    for (i in 1:M){
      mod41 = multinom(result~set41[,paste0("n",i)]+S+round+S*round , data = set41 )   
      m41_nh[2*(i-1)+1,] = summary(mod41)$coefficients[c(1,3,5,7,9)]
      m41_nh[2*i,] = summary(mod41)$coefficients[c(2,4,6,8,10)]
    } 
    
    m41_nh[is.na(m41_nh)]=0
    
    mod30_pa = glm(result~ S, data = set30, family = binomial)
    m30_pa = mod30_pa$coefficients
    
    mod31_pa = glm(result~ S, data = set31, family = binomial)
    m31_pa = mod31_pa$coefficients
    
    mod40_pa = multinom(result~S , data = set40)
    m40_pa = summary(mod40_pa)$coefficients 
    
    mod41_pa = multinom(result~S , data = set41)
    m41_pa = summary(mod41_pa)$coefficients
    
    functions2_sc1 = list (round_11= round_11(m1, m2), P_prime2= P_prime2(m40, m41), P2=  P2(m30 , m31))
    functions1_sc1 = list(round_11_pa = round_11_pa(m1,m2), P2_pa = P2_pa(m30_pa, m31_pa), P_prime2_pa= P_prime2_pa(m40_pa, m41_pa))
    functions_sc1 = list (round_11= round_11(m1, m2),
                          P_prime2= P_prime2(m40, m41),
                          P2=  P2(m30 , m31),
                          round_1_nh= round_1_nh(m1_nh, m2_nh),
                          P_nh=  P_nh(m30_nh, m31_nh),
                          P_prime_nh= P_prime_nh(m40_nh, m41_nh))
    
    g = data6[data6$fp2 <= data6$S,] 
    g.dat = NULL
    for (j in 2:(M-1)) 
      g.dat = tryCatch(rbind(g.dat,cbind( g[g$round <= j & g$S >=j  ,],j)), error = function(e){g.dat})
    
    g.dat = tryCatch(g.dat[g.dat$round <= g.dat$j ,],error=function(e){g.dat})
    
    g.dat$alpha =(g.dat$round - (g.dat$j+1))
    
    G = glm((S == j) ~ alpha, family = "binomial", data = g.dat)$coef[2]
    
    cum_ds = disc.surv(data6,M,l)
    cum_pa = pop.ave_Markov(data6, functions1_sc1, functions2_sc1, l, M) 
    cum_cb = cens.bias_Markov(data6, M, G , l, functions_sc1)
    
    out1[jj,1] = cum_ds$risk_ds
    out1[jj,2] = cum_pa$risk_pa
    out1[jj,3] = cum_cb$risk
    se1[jj,1] = cum_ds$se_ds
    se1[jj,2] = cum_pa$se_pa
    se1[jj,3] = cum_cb$se
  }
  return(list(cumr = apply(out1, 2, mean), tstr = apply(se1, 2, mean), estr = apply(out1, 2, sd), qtls = apply(out1, 2, quantile, probs = c(0,0.25, 0.75, 0.975, 1),  na.rm = TRUE)))
  
}

#this function is used for simulation in scenario 5,6,7 and 8
MfpRisk_sim2 = function(Nsim, M, p, l, Nsubj, scen, P_prime, d="strong"){
  
  out1 = matrix (ncol=3, nrow= Nsim)
  se1 = matrix (ncol=3, nrow= Nsim)
  for (jj in 1:Nsim){
    #creating round variable 
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
    data3 = ddply(.data = data2, .variable ='StudyID_c', .fun = function(df){
      df$fp1 = ifelse(is.na(which(df$result==1)[1]), M+1, which(df$result==1)[1])
      return(df)
    })
    
    #add the variable for the round of second fp (this is the event time for at least 2 FP)
    data4 = ddply(.data = data3, .variable ='StudyID_c', .fun = function(df){
      df$fp2 = ifelse(is.na(which(df$result==1)[2]), M+1, which(df$result==1)[2])
      return(df)
    })
    
    fp2 = ddply(.data=data4, .variable ='StudyID_c', .fun = function(df){ 
      fp2 = unique(df$fp2)
      return(fp2)
    })
    
    #define the round of competing event (first 2)
    data5 = ddply(.data = data4, .variable ='StudyID_c', .fun = function(df){
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
    data5$S = rep(S,each = M)
    data5$S = ifelse(data5$S<=data5$comp, data5$S, data5$comp)
    
    #add M variables ni when ni is the number of false positive by and including round i for i=1,...,M
    data6 = ddply(.data= data5, .variables = 'StudyID_c', .fun = function(df) {
      for (i in 1:M){
        df[, paste0("n",i)] = sum(df$result[1:i]==1)
      }
      return(df)
    })
    
    #add variable for the censoring if first fp is the event (min(S,fp1))
    data6$cens1 = apply(cbind(data6$S,data6$fp1),1,min)  
    
    #add variable for the censoring if second fp is the event (min(S,fp2))
    data6$cens2 = apply(cbind(data6$S,data6$fp2),1,min)  
    
    #define event time 
    data6$fp2.cens = ifelse(data6$fp2>data6$S,0,data6$fp2) 
    
    #define censoring indicator 
    data6$delta = ifelse(data6$fp2 <= data6$S, 1, 0)
    
    set1 = subset(data6,  1==round & round ==S)
    set2 = subset(data6,  1==round & round <= (S-1))
    set3 = subset(data6,  2<=round & round <= (S-1))
    set4 = subset(data6,  2<=round & round ==S)
    set30 = subset(set3, s1==0)
    set31 = subset(set3, s1==1)
    set40 = subset(set4, s1==0)
    set41 = subset(set4, s1==1)
    mod1 = multinom(result~1 , data = set1)
    m1 = summary(mod1)$coefficients 
    m1[is.na(m1)]=0
    mod2 = glm(result~ S , data = set2, family=binomial)
    m2 = mod2$coefficients 
    m2[is.na(m2)]=0
    mod30 = glm(result~ S+round+S*round , data = set30, family=binomial)
    m30 = mod30$coefficients 
    m30[is.na(m30)]=0
    mod31 = glm(result~ S+round+S*round , data = set31, family=binomial)
    m31 = mod31$coefficients 
    m31[is.na(m31)]=0
    mod40 = multinom(result~S+round+S*round , data = set40)
    m40 = summary(mod40)$coefficients 
    m40[is.na(m40)]=0
    mod41 = multinom(result~S+round+S*round , data = set41)
    m41 = summary(mod41)$coefficients 
    m41[is.na(m41)]=0
    m1_nh = matrix (0, nrow = 2*M, ncol=2) 
    for (i in 1:M){
      mod1 = multinom(result~ set1[,paste0("n",i)] , data = set1)
      m1_nh[2*(i-1)+1,] = summary(mod1)$coefficients[c(1,3)]
      m1_nh[2*i,] = summary(mod1)$coefficients[c(2,4)]
    }
    
    m1_nh[is.na(m1_nh)]=0
    
    m2_nh = matrix (0, nrow = M, ncol=3) 
    for (i in 1:M){
      mod2 = glm(result~ S+set2[,paste0("n",i)] , data = set2, family=binomial) 
      m2_nh[i,] = mod2$coefficients
    }
    
    m2_nh[is.na(m2_nh)]=0
    
    m30_nh = matrix (0, nrow = M, ncol=5)
    for (i in 1:M){
      mod30 = glm(result~set30[,paste0("n",i)]+S+round+S*round , data = set30, family=binomial)#row 1 of matrix (7)
      m30_nh[i,] = mod30$coefficients
    }
    m30_nh[is.na(m30_nh)]=0
    
    m31_nh = matrix(0,nrow = M, ncol=5) 
    for (i in 1:M){
      mod31 = glm(result~set31[,paste0("n",i)]+ S+round+S*round , data =set31 , family=binomial)#row 2 of matrix (7)
      m31_nh[i,] = mod31$coefficients
    }
    
    m31_nh[is.na(m31_nh)]=0
    
    m40_nh = matrix(0, nrow = 2*M, ncol = 5)
    for (i in 1:M){
      mod40 = multinom(result~set40[,paste0("n",i)]+S+round+S*round , data =set40 )  
      m40_nh[2*(i-1)+1,] = summary(mod40)$coefficients[c(1,3,5,7,9)]
      m40_nh[2*i,] = summary(mod40)$coefficients[c(2,4,6,8,10)]
    }
    
    m40_nh[is.na(m40_nh)]=0
    
    m41_nh = matrix(0,nrow = 2*M, ncol = 5) 
    for (i in 1:M){
      mod41 = multinom(result~set41[,paste0("n",i)]+S+round+S*round , data = set41 ) #row 2 of matrix (10)  
      m41_nh[2*(i-1)+1,] = summary(mod41)$coefficients[c(1,3,5,7,9)]
      m41_nh[2*i,] = summary(mod41)$coefficients[c(2,4,6,8,10)]
    } 
    
    m41_nh[is.na(m41_nh)]=0
    
    mod30_pa = glm(result~ S, data = set30, family = binomial)
    m30_pa = mod30_pa$coefficients
    
    mod31_pa = glm(result~ S, data = set31, family = binomial)
    m31_pa = mod31_pa$coefficients
    
    mod40_pa = multinom(result~S , data = set40)
    m40_pa = summary(mod40_pa)$coefficients 
    
    mod41_pa = multinom(result~S , data = set41)
    m41_pa = summary(mod41_pa)$coefficients
    
    functions2_sc1 = list (round_11= round_11(m1, m2), P_prime2= P_prime2(m40, m41), P2=  P2(m30 , m31))
    functions1_sc1 = list(round_11_pa = round_11_pa(m1,m2), P2_pa = P2_pa(m30_pa, m31_pa), P_prime2_pa= P_prime2_pa(m40_pa, m41_pa))
    functions_sc1 = list (round_11= round_11(m1, m2),
                          P_prime2= P_prime2(m40, m41),
                          P2=  P2(m30 , m31),
                          round_1_nh= round_1_nh(m1_nh, m2_nh),
                          P_nh=  P_nh(m30_nh, m31_nh),
                          P_prime_nh= P_prime_nh(m40_nh, m41_nh))
    
    g = data6[data6$fp2 <= data6$S,] 
    g.dat = NULL
    for (j in 2:(M-1)) 
      g.dat = tryCatch(rbind(g.dat,cbind( g[g$round <= j & g$S >=j  ,],j)), error = function(e){g.dat})
    
    g.dat = tryCatch(g.dat[g.dat$round <= g.dat$j ,],error=function(e){g.dat})
    
    g.dat$alpha =(g.dat$round - (g.dat$j+1))
    
    G = glm((S == j) ~ alpha, family = "binomial", data = g.dat)$coef[2]
    
    cum_ds = disc.surv(data6,M,l)
    cum_pa = pop.ave_Markov(data6, functions1_sc1, functions2_sc1, l, M) 
    cum_cb = cens.bias_Markov(data6, M, G , l, functions_sc1)
    
    out1[jj,1] = cum_ds$risk_ds
    out1[jj,2] = cum_pa$risk_pa
    out1[jj,3] = cum_cb$risk
    se1[jj,1] = cum_ds$se_ds
    se1[jj,2] = cum_pa$se_pa
    se1[jj,3] = cum_cb$se
  }
  return(list(cumr = apply(out1, 2, mean), tstr = apply(se1, 2, mean), estr = apply(out1, 2, sd), qtls = apply(out1, 2, quantile, probs = c(0,0.25,0.50, 0.75, 0.975, 1),  na.rm = TRUE)))
  
}


                  ################################################################
                  #                         scenario 1:                          #
                  #   false positive risk is independent of censoring time,      #
                  #   number of prior false positive results and is constant     #
                  #                  across screening rounds                     #
                  ################################################################

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

#this function returns the p_ij ^(t,h) in the manuscript 
trnMt = function(i,j,h,t){
  
  if (2<=t & t< h)
    return(P[i+1, j+1])
  
  if (t>=h & t*h!=1)
    return(P_prime[i+1, j+1])
}


#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2)){
    c = c*trnMt(x[k],x[k+1],h,k+1)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least two false positive after 10 rounds
T1= TrueFp(2,10)
#0.3187748

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 1)

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

#this function returns the p_ij ^(t,h) in the manuscript 
trnMt = function(i,j,h,t){
  
  if (2<=t & t< h)
    return(P(t)[i+1, j+1])
  
  if (t>=h & t*h!=1)
    return(P_prime(t)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2))
    c = c*trnMt(x[k],x[k+1],h,k+1)
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T2_s=TrueFp(2,10)
#0.2012562

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 2)

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

trnMt = function(i,j,h,t){
  if (2<=t & t< h)
    return(P(t)[i+1, j+1])
  
  if (t>=h & t*h!=1)
    return(P_prime(t)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2))
    c = c*trnMt(x[k],x[k+1],h,k+1)
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T2_m=TrueFp(2,10)
#0.2790235

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 2)

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

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum((x[1:k] == 1))
    c = c*trnMt(x[k],x[k+1],h,k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T3_s= TrueFp(2,10)
#0.3242431

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 3)


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

trnMt = function(i,j,h,t,v){
  
  if (2<=t & t< h)
    return(P(v)[i+1, j+1])
  
  if (t>=h & t*h!=1)
    return(P_prime(v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum((x[1:k] == 1))
    c = c*trnMt(x[k],x[k+1],h,k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)


TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T3_m= TrueFp(2,10)
#0.3129076

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 3)


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

trnMt = function(i,j,h,t,v){
  if (2<=t & t< h)
    return(P(t, v)[i+1, j+1])
  if (t>=h & t*h!=1)
    return(P_prime(t, v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],h,k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

createCartesian0 = function(l, j){
  L = list()
  for (i in 1:(j-1))
    L[[i]] = c(0,1)
  if (j ==1){
    df = data.frame('i1' = 1)
  } else{
    df = expand.grid(L)
    for (i in 1:(j-1))
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var==l-1)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==j-ind)})
    df = subset(df, valid==1)
    df$valid = NULL
    df[, paste0('i',j)] = 1}
  return (df)
}

TrueFp = function(r,M){
  init.P = 0
  for (ii in r:M){
    df = createCartesian0(r,ii)
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    init.P = init.P + sum(df$term)
  }
  return (init.P)
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T4_s=TrueFp(2,10)
#0.2782955

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 4)

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

trnMt = function(i,j,h,t,v){
  
  if (2<=t & t< h)
    return(P(t, v)[i+1, j+1])
  
  if (t>=h & t*h!=1)
    return(P_prime(t, v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1], h)*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],h,k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

createCartesian0 = function(l, j){
  L = list()
  for (i in 1:(j-1))
    L[[i]] = c(0,1)
  if (j ==1){
    df = data.frame('i1' = 1)
  } else{
    df = expand.grid(L)
    for (i in 1:(j-1))
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var==l-1)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==j-ind)})
    df = subset(df, valid==1)
    df$valid = NULL
    df[, paste0('i',j)] = 1}
  return (df)
}

TrueFp = function(r,M){
  init.P = 0
  for (ii in r:M){
    df = createCartesian0(r,ii)
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    init.P = init.P + sum(df$term)
  }
  return (init.P)
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T4_m=TrueFp(2,10)
#0.328788

MfpRisk_sim1 (5000, 10, 0.2, 2, 50000, P, P_prime, 4)

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

trnMt = function(i,j,t,v){
  if (t>1)
    return(P_prime(t, v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)


TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T5_s=TrueFp(2,10)
# 0.2686583

MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 5, P_prime, d="strong")

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

trnMt = function(i,j,t,v){
  if (t>1)
    return(P_prime(t, v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T5_m=TrueFp(2,10)
# 0.3162108

MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 5, P_prime, d="moderate")
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


trnMt = function(i,j,t,v){
  if (t>1)
    return(P_prime(v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T6_s=TrueFp(2,10)
# 0.3163313

MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 6, P_prime, d="strong")

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


trnMt = function(i,j,t,v){
  if (t>1)
    return(P_prime(v)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],k+1, v)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T6_m=TrueFp(2,10)

MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 6, P_prime, d="moderate")

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

trnMt = function(i,j){
  return(P_prime[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    c = c*trnMt(x[k],x[k+1])}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)


TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T7_ms =TrueFp(2,10)
# 0.3057751

MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 7, P_prime, d="strong")
MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 7, P_prime, d="moderate")


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

trnMt = function(i,j,t){
  if (t>1)
    return(P_prime(t)[i+1, j+1])
}

#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    v = sum(x[1:k] == 1)
    c = c*trnMt(x[k],x[k+1],k+1)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T8_s =TrueFp(2,10)
# 0.2678492
MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 8, P_prime, d="strong")

              ################################################################
              #                         scenario 8:                          #
              #         false positive risk depends on censoring             #
              #     time and screening round but is independent              #
              #       of number of prior false positive results              #   
              #           Case II: Moderate dependency                       #
              ################################################################


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

#this function returns the p_ij ^(h,t) 
trnMt = function(i,j,t){
  if (t>1)
    return(P_prime(t)[i+1, j+1])
}
#Calculating the true cumulative risk of at least two false positive after 10 rounds
f0 = function(x){
  x = as.numeric(x)
  h = x[length(x)]
  c = round_1(x[1])*Pr[h]
  for (k in 1:(length(x)-2)){
    c = c*trnMt(x[k],x[k+1],k+1)}
  return (c)
}

L = list()
for (i in 1:M)
  L[[i]] = c(0,1,2)

TrueFp = function(r,M){
  for (ii in r:M){
    df = expand.grid(L)
    for (i in 1:M)
      names(df)[i] = paste0('i',i)
    df$var = apply(df, 1, function(x) length(which(x==1)))
    df = subset(df, var>=ii)
    df$var = NULL
    df$valid = apply(df, 1, function(x){
      x = as.numeric(x)
      ind = which(x==2)[1]
      ifelse(is.na(ind), 1, sum(x==2)==M-ind+1)})
    df = subset(df, valid==1)
    df$valid = NULL
    df = merge(df, data.frame('h'=1:M))
    df$term = apply(df,1, f0)
    return(sum(df$term))
  }
}

# True value of cumulative risk of at least 2 FP after 10 rounds
T8_m =TrueFp(2,10)
# 0.3153073

MfpRisk_sim2 (5000, 10, 0.2, 2, 50000, 8, P_prime, d="moderate")







