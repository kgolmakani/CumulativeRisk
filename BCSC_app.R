#Application to BCSC

library("plyr")
library("dplyr")
library("nnet")
library("ggplot2")

MFP = read.csv("C:\\Users\\golmmk1\\Desktop\\MFP_dat_new.csv", sep=",", header=TRUE)
#examdate is the date of current screening mammogram
MFP$examdate = as.Date(strptime(x = as.character(MFP$examdate), format= "%d%b%Y"))
#nxtscrdt_c is the date of next screening mammogram
MFP$nxtscrdt_c = as.Date(strptime(x = as.character(MFP$nxtscrdt_c), format= "%d%b%y"))
#nxtscrdt_c is the date of next-next screening mammogram
MFP$nxtscrdt_c_lead = as.Date(strptime(x = as.character(MFP$nxtscrdt_c_lead), format= "%d%b%y"))
#newdxdt is the date of woman's earliest breast cancer diagnosis
MFP$newdxdt= as.Date(strptime(x = as.character(MFP$newdxdt), format= "%d%b%Y"))
#eoccadt is the date of complete cancer capture
MFP$eoccadt = as.Date(strptime(x = as.character(MFP$eoccadt), format= "%d%b%y"))
#creating variable for time until next screen (in months)
interval_n = floor(difftime(MFP$nxtscrdt_c, MFP$examdate, units="days")/30.5)
interval_n =as.numeric(interval_n)
nterval_ne = ifelse (is.na(interval_n), 300,interval_n)
MFP$interval_ne = interval_ne
#ordering the data set
MFP1 = MFP[with(MFP, order(StudyID_c, examdate)),]
#0=true negative (TN), 1=false positive (FP), 2=Cancer (Ca) 
#(we use initial result and cancer within a year of index exam with 
#cut-off at next screen for creating TN, FP and ca)
MFP1$result = ifelse(MFP1$resinit_c==0 & MFP1$cancscrfu1yr_c==0,0,(ifelse(MFP1$resinit_c==1 & MFP1$cancscrfu1yr_c==0,1,(ifelse(MFP1$cancscrfu1yr_c==1,2,NA)))))
MFP1$id1 = factor(MFP1$StudyID_c, levels=unique(MFP1$StudyID_c))
#creating censoring time (S=censoring time: total number of rounds attended) and screening round
MFP1_c =  cbind(MFP1, do.call(rbind,tapply(MFP1$result, MFP1$id1, function(v){
  indx = which(v==1) 
  if (length(indx)>0)
    v = rep(indx[1],length(v))
  S1 = rep(length(v),length(v))
  rd = seq(length(v))
  as.data.frame(list(S1=S1,rd=rd))
})))

#considering first round
round1 = subset(MFP1_c, rd==1) 
nrow(round1)
#we just consider the subjects who their first screen is their first screen ever 
first_ever = subset(round1, round1$prvmam_c== 0)

st_first_ever = unique(first_ever$StudyID_c)

dat_MFP1 = subset(MFP1_c, MFP$StudyID_c %in% st_first_ever)

dat_MFP = dat_MFP1[with(dat_MFP1, order(StudyID_c, examdate)),]
#prvmam_c captures whether or not there is any evidence of a 
#previous mammogram based on all available sources
dat_MFP$prvmam_c = as.factor(dat_MFP$prvmam_c)

levels(dat_MFP$prvmam_c) #0 : No previous mammogram , 1: evidence of prior mamogram, 9 : unknown

dat_MFP1 = subset(dat_MFP, (resinit_c!="8" & resinit_c!="9"))

# for time since last screen we define first, annual and biennial and triennial
#0: first , 1 : annual (9-18 months), 2 : biennial (19-30 months), 3: triennial(31-42 months), 
#4: longer (>42 months), 5: missing time since last mammogram
interval_l = ifelse(dat_MFP1$prvmam_c=="0", 0, ifelse(is.na(dat_MFP1$prvmos_c),5,(ifelse((9<=dat_MFP1$prvmos_c & dat_MFP1$prvmos_c<=18), 1, (ifelse((19<=dat_MFP1$prvmos_c & dat_MFP1$prvmos_c<=30),2,(ifelse((31<=dat_MFP1$prvmos_c & dat_MFP1$prvmos_c<=42),3,4))))))))

interval_l = as.factor(interval_l)

dat_MFP1$interval_l = interval_l
# We just consider first, annual and biennial screeners
dat_MFP2 = subset(dat_MFP1, (interval_l==0 | interval_l==1 | interval_l==2|(rd==1 & prvmam_c==0)))

dat_MFP2$interval_l=factor(dat_MFP2$interval_l)
#first screeners ever 
dat_MFP2_first = subset(dat_MFP2, (interval_l==0|(rd==1 & prvmam_c==0)))
#just considering biennial screening
dat_MFP2_bi = subset(dat_MFP2, interval_l==2)
#just considering annual screening
dat_MFP2_ann = subset(dat_MFP2, interval_l==1)

options(max.print=1000000)

#biennial screens with the complete cancer capture
d=dat_MFP2_bi$eoccadt

d=as.POSIXlt(d)

d$year=d$year-2

dat_MFP2_bi$comp_cap = as.Date(d)

dat_MFP2_bi_c= subset(dat_MFP2_bi, examdate <= comp_cap)

#annual screens with complete cancer capture
d_ann=dat_MFP2_ann$eoccadt

d_ann=as.POSIXlt(d_ann)

d_ann$year=d_ann$year-1

dat_MFP2_ann$comp_cap = as.Date(d_ann)

dat_MFP2_ann_c= subset(dat_MFP2_ann, examdate <= comp_cap)

#first screen with complete cancer capture
d_first=dat_MFP2_first$eoccadt

d_first=as.POSIXlt(d_first)

d_first$year=d_first$year

dat_MFP2_first$comp_cap = as.Date(d_first)

dat_MFP2_first_c= subset(dat_MFP2_first, examdate <= comp_cap)

subset(dat_MFP2_first, examdate > comp_cap)

dat_MFP2_c= rbind(dat_MFP2_bi_c, dat_MFP2_ann_c, dat_MFP2_first_c)

dat_MFP2_c1 = dat_MFP2_c[with(dat_MFP2_c, order(StudyID_c, examdate)),]

#creating the response variable : categorical with 3 values: 0,1,2
#0=TN: A negative screen with no cancer diagnosis within 1 year of the index screen and before the next screen
#1=FP: A positive screen with no cancer diagnosis within 1 year of the index screen and before the next screen
#2=Ca: Cancer diagnosis within 1 year of the index screen and before the next screen (regardless of the result of the index screen)
dat_MFP2_c1$result = ifelse(dat_MFP2_c1$resinit_c==0 & dat_MFP2_c1$cancscrfu1yr_c==0,0,(ifelse(dat_MFP2_c1$resinit_c==1 & dat_MFP2_c1$cancscrfu1yr_c==0,1,(ifelse(dat_MFP2_c1$cancscrfu1yr_c==1,2,NA)))))

#creating censoring time and screening rounds
dat_MFP2_c1$id = factor(dat_MFP2_c1$StudyID_c, levels=unique(dat_MFP2_c1$StudyID_c))

dat_MFP_c =  cbind(dat_MFP2_c1, do.call(rbind,tapply(dat_MFP2_c1$result, dat_MFP2_c1$id, function(v){
  indx = which(v==1) 
  if (length(indx)>0)
    v = rep(indx[1],length(v))
  S = rep(length(v),length(v))
  round = seq(length(v))
  as.data.frame(list(S=S,round=round))
})))

# creating s1 variable (previous result for first order Markov chain)
# for n-th order markov chain we need n variable for n previous results
addCol = function(df){
  n = nrow(df)
  df = df[order(df$StudyID_c,df$round),]
  df$s1 = c(NA,df$result[-n])
  df$s1[match(unique(df$StudyID_c),df$StudyID_c )] = NA 
  return (df)
}

dat_MFP_c_f = addCol(dat_MFP_c)

dat_MFP_c_f$result = as.factor(dat_MFP_c_f$result)

dat_MFP_c_f$s1 = as.factor(dat_MFP_c_f$s1)

dat_MFP_c_f = subset(dat_MFP_c_f, round<=M)

dat_MFP_c_f$S = ifelse(dat_MFP_c_f$S >M, M , dat_MFP_c_f$S)

dat0 = dat_MFP_c_f[, c("StudyID_c", "result","S", "round")]

dat0$ID = seq(1:nrow(dat0))

dat0$StudyID_c = as.character(dat0$StudyID_c)
#add the variable for the round of the first fp
dat1 = ddply(.data =dat0 , .variable ='StudyID_c', .fun = function(df){
  df$fp1 = ifelse(is.na(which(df$result==1)[1]), M+1, which(df$result==1)[1])
  return(df)
})
#add the variable for the round of second fp (this is the event time for at least 2 FP)
dat2 = ddply(.data = dat1, .variable ='StudyID_c', .fun = function(df){
  df$fp2 = ifelse(is.na(which(df$result==1)[2]), M+1, which(df$result==1)[2])
  return(df)
})
#add the variable for the round of the third fp
dat3 = ddply(.data =dat2 , .variable ='StudyID_c', .fun = function(df){
  df$fp3 = ifelse(is.na(which(df$result==1)[3]), M+1, which(df$result==1)[3])
  return(df)
})
#add the variable for the round of the forth fp
dat4 = ddply(.data =dat3 , .variable ='StudyID_c', .fun = function(df){
  df$fp4 = ifelse(is.na(which(df$result==1)[4]), M+1, which(df$result==1)[4])
  return(df)
})
#add the variable for the round of the fifth fp
dat5 = ddply(.data =dat4 , .variable ='StudyID_c', .fun = function(df){
  df$fp5 = ifelse(is.na(which(df$result==1)[5]), M+1, which(df$result==1)[5])
  return(df)
})

#add variable for the number of prior FPs
dat6 = ddply(.data =dat5, .variables ='StudyID_c', .fun = function(df){
  x= ifelse(df$result==1,1,0)
  df[,'prior_fp'] = c(0,cumsum(x[-length(x)]))
  return(df)
} )


dat6 = dat6[with(dat6, order(ID)),]

dat6 = dat6[,c("fp1", "fp2","fp3","fp4","fp5", "prior_fp")]

dat_MFP_c_f = cbind(dat_MFP_c_f, dat6)

dat_MFP_c_f = dat_MFP_c_f[with(dat_MFP_c_f, order(StudyID_c, examdate)),]

#function to create baseline variable
addFirstVariable = function(df, var){
  df = df[order(df$StudyID_c, df$round),]
  indx = match(unique(df$StudyID_c),df$StudyID_c )
  indx2 = indx[-1] - 1
  x=c()
  for (i in 1:length(indx2))
    x = c(x,indx[i],indx2[i])
  x = c(x, indx[length(indx)], nrow(df))
  for (ii in 1:(length(x)/2))
    df[x[2*ii-1]:x[2*ii], paste0('first',var)] = df[x[2*ii-1], var]
  return (df)
}

dat_MFP_c_f = addFirstVariable(dat_MFP_c_f, "age_c")

dat_MFP_c_f$firstage_cat = cut(dat_MFP_c_f$firstage_c, breaks =seq(40, 75, by=5), include.lowest = TRUE, right = FALSE)

dat_MFP_c_f = addFirstVariable(dat_MFP_c_f, "famhx_c")

dat_MFP_c_f = addFirstVariable(dat_MFP_c_f, "benbiopresult_c")

dat_MFP_c_f = addFirstVariable(dat_MFP_c_f, "density_c")
#categorizing bmi variable
dat_MFP_c_f$bmi_cat = ifelse(dat_MFP_c_f$bmi_c <18.5, 1, (ifelse(18.5 <=dat_MFP_c_f$bmi_c & dat_MFP_c_f$bmi_c <=24.9, 2,(ifelse(25 <=dat_MFP_c_f$bmi_c & dat_MFP_c_f$bmi_c <=29.9, 3, (ifelse(30 <=dat_MFP_c_f$bmi_c & dat_MFP_c_f$bmi_c <=34.9 , 4, 5)))))))

dat_MFP_c_f = addFirstVariable(dat_MFP_c_f, "bmi_cat")
#creating a variable for first live birth
dat_MFP_c_f$age1stb_cat = ifelse(dat_MFP_c_f$age1stb==1|dat_MFP_c_f$age1stb==2|dat_MFP_c_f$age1stb==3|dat_MFP_c_f$age1stb==7| (9<=dat_MFP_c_f$age1stb &dat_MFP_c_f$age1stb<30), 1, dat_MFP_c_f$age1stb)

dat_MFP_c_f = addFirstVariable(dat_MFP_c_f, "age1stb_cat")

dat_MFP_c_f = subset(dat_MFP_c_f, !is.na(firstfamhx_c) & !is.na(firstdensity_c) & !is.na(racenci_c) )

dat_MFP_c_f$racenci_c = ifelse(dat_MFP_c_f$racenci_c==3 |dat_MFP_c_f$racenci_c==4, 3, dat_MFP_c_f$racenci_c)
dat_MFP_c_f$racenci_c = ifelse(dat_MFP_c_f$racenci_c==6|dat_MFP_c_f$racenci_c==7|dat_MFP_c_f$racenci_c==9, 9, dat_MFP_c_f$racenci_c)

dat_MFP_c_f$firstage1stb_cat = ifelse(dat_MFP_c_f$firstage1stb_cat>=30 & dat_MFP_c_f$firstage1stb_cat <88 ,2, dat_MFP_c_f$firstage1stb_cat)
dat_MFP_c_f$firstage1stb_cat =ifelse(dat_MFP_c_f$firstage1stb_cat==4|dat_MFP_c_f$firstage1stb_cat==5|dat_MFP_c_f$firstage1stb_cat==6|dat_MFP_c_f$firstage1stb_cat==8, 2, dat_MFP_c_f$firstage1stb_cat)
dat_MFP_c_f$firstage1stb_cat = ifelse(dat_MFP_c_f$firstage1stb_cat==88|dat_MFP_c_f$firstage1stb_cat==99, NA, dat_MFP_c_f$firstage1stb_cat)

dat_MFP_c_f$firstfamhx_c = ifelse(dat_MFP_c_f$firstfamhx_c==9, NA, dat_MFP_c_f$firstfamhx_c)

dat_MFP_c_f$firstage_cat = as.character(dat_MFP_c_f$firstage_cat)

dat_MFP_c_f = subset(dat_MFP_c_f, firstage_cat!="[60,65)" & firstage_cat!="[65,70)" & firstage_cat!="[70,75]" )

dat_MFP_c_f$firstage_cat = ifelse(dat_MFP_c_f$firstage_cat=="[40,45)"|dat_MFP_c_f$firstage_cat=="[45,50)", "[40,50)",dat_MFP_c_f$firstage_cat )
dat_MFP_c_f$firstage_cat = ifelse(dat_MFP_c_f$firstage_cat=="[50,55)"|dat_MFP_c_f$firstage_cat=="[55,60)", "[50,60)",dat_MFP_c_f$firstage_cat )
dat_MFP_c_f$firstage_cat = as.factor(dat_MFP_c_f$firstage_cat)

women_round1 = subset(dat_MFP_c_f, round==1)

women_subsequent =subset(dat_MFP_c_f, round>1 & round<=10)

FP = subset(women_round1, result==1)
(nrow(FP)/nrow(women_round1))*100

FP_s = subset(women_subsequent, result==1)
(nrow(FP_s)/nrow(women_subsequent))*100
#distribution of baseline age
sum(table(women_round1$firstage_cat))

(table(women_round1$firstage_cat)/nrow(women_round1))*100
#percentage of women with FP at first exam among each category
FP = subset(women_round1, result==1)
(table(FP$firstage_cat)/table(women_round1$firstage_cat))*100

FP_s = subset(women_subsequent, result==1)
(table(FP_s$firstage_cat)/table(women_subsequent$firstage_cat))*100
#distribution of racenci_c
sum(table(women_round1$racenci_c))

(table(women_round1$racenci_c)/nrow(women_round1))*100
#percentage of women with FP at first exam among each category
FP = subset(women_round1, result==1)
(table(FP$racenci_c)/table(women_round1$racenci_c))*100

FP_s = subset(women_subsequent, result==1)
(table(FP_s$racenci_c)/table(women_subsequent$racenci_c))*100
#distribution of family history of breast cancer at baseline
table(women_round1$firstfamhx_c)

(table(women_round1$firstfamhx_c)/nrow(women_round1))*100
#percentage of women with FP at first exam among each category
FP = subset(women_round1, result==1)
(table(FP$firstfamhx_c)/table(women_round1$firstfamhx_c))*100

FP_s = subset(women_subsequent, result==1)
(table(FP_s$firstfamhx_c)/table(women_subsequent$firstfamhx_c))*100
#distribution of breast density at baseline
table(women_round1$firstdensity_c)

(table(women_round1$firstdensity_c)/nrow(women_round1 ))*100
#percentage of women with FP at first exam among each category
FP = subset(women_round1, result==1)
(table(FP$firstdensity_c)/table(women_round1$firstdensity_c))*100

FP_s = subset(women_subsequent, result==1)
(table(FP_s$firstdensity_c)/table(women_subsequent$firstdensity_c))*100

women_round1$total_round = ifelse(women_round1$S>=5, 5, women_round1$S)
#distribution of total rounds attended
table(women_round1$total_round)

(table(women_round1$total_round)/nrow(women_round1))*100

#percentage of women with FP at first exam among each category
FP = subset(women_round1, result==1)
(table(FP$total_round)/table(women_round1$total_round))*100

women_subsequent$total_round = ifelse(women_subsequent$S>=5, 5, women_subsequent$S)
FP_s = subset(women_subsequent, result==1)
(table(FP_s$total_round)/table(women_subsequent$total_round))*100
#distribution of reason for censoring
table(women_round1$cancens_cat_c)
(table(women_round1$cancens_cat_c)/nrow(women_round1))*100
#creating a binary variable for result
dat_MFP_c_f$result2 = ifelse(dat_MFP_c_f$result==1, 1,0)

dat_MFP_c_f$racenci_c = as.factor(dat_MFP_c_f$racenci_c)

dat_MFP_c_f$firstbmi_cat = as.factor(dat_MFP_c_f$firstbmi_cat)

dat_MFP_c_f$firstbmi_cat = relevel (dat_MFP_c_f$firstbmi_cat, ref = "2")

dat_MFP_c_f$firstage1stb_cat = as.factor(dat_MFP_c_f$firstage1stb_cat)

dat_MFP_c_f$firstage1stb_cat = relevel (dat_MFP_c_f$firstage1stb_cat, ref = "1") 

dat_MFP_c_f$firstdensity_c = as.factor(dat_MFP_c_f$firstdensity_c)

dat_MFP_c_f$firstdensity_c = relevel (dat_MFP_c_f$firstdensity_c, ref = "2") 

dat_MFP_c_f$firstage_cat = as.factor(dat_MFP_c_f$firstage_cat)

dat_MFP_c_f$firstfamhx_c  = as.factor(dat_MFP_c_f$firstfamhx_c )

library("rms")
library("contrast")

#Odds ratio when the outcome is false positive
cov = c("racenci_c", "firstdensity_c", "firstfamhx_c", "age_c")
df = data.frame('var' = as.character(0), 'level'=as.character(0),
                'exp.contrast'=as.numeric(0), 'exp.lower'=as.numeric(0),'exp.upper'=as.numeric(0))

mod=glm(result2~racenci_c+firstdensity_c+firstfamhx_c+age_c, data= dat_MFP_c_f, family= binomial(link="logit"))

exp(cbind("Odds ratio" = coef(mod), confint.default(mod, level = 0.95)))
for (var in cov){
  if (is.factor(dat_MFP_c_f[,var])) {
    lv = levels(dat_MFP_c_f[,var])
    for(j in lv[-1]){
      list1 = list()
      list1[[var]] = j
      list2 = list()
      list2[[var]] = lv[1]
      diff = setdiff(cov, var)
      for (vri in diff){
        list1[[vri]]=dat_MFP_c_f[1,vri] 
        list2[[vri]]=dat_MFP_c_f[1,vri]
      }
      OR = contrast(mod, list1, list2)
      df = rbind(df, data.frame(list('var'=as.character(var),'level'=as.character(j),'exp.contrast'=as.numeric(exp(OR$Contrast)),
                                     'exp.lower'=as.numeric(exp(OR$Lower)),'exp.upper'=as.numeric(exp(OR$Upper)))))
    }
    
  }else{
    list1 = list()
    list1[[var]] = dat_MFP_c_f[1,var]+1
    list2 = list()
    list2[[var]] = dat_MFP_c_f[1,var]
    diff = setdiff(cov, var)
    for (vri in diff){
      list1[[vri]] = dat_MFP_c_f[1,vri]
      list2[[vri]] = dat_MFP_c_f[1,vri]
    }
    OR = contrast(mod, list1,list2)
    
    df=rbind(df,data.frame(list('var'=var, 'level'='NA','exp.contrast'=exp(OR$Contrast),
                                'exp.lower'=exp(OR$Lower),'exp.upper'=exp(OR$Upper))))
  }
}
df = df[-1,]
write.csv(df, file="~/Desktop/OR_combined_fp_MFP")

#Odds raio when the outcome is cancer
dat_MFP_c_f$result3 = ifelse(dat_MFP_c_f$result==2, 1,0)

cov = c("racenci_c",  "firstdensity_c", "firstfamhx_c", "age_c"  )
df = data.frame('var' = as.character(0), 'level'=as.character(0),
                'exp.contrast'=as.numeric(0), 'exp.lower'=as.numeric(0),'exp.upper'=as.numeric(0))

mod=glm(result3~racenci_c+firstdensity_c+firstfamhx_c+age_c, data= dat_MFP_c_f, family= binomial(link="logit"))

exp(cbind("Odds ratio" = coef(mod), confint.default(mod, level = 0.95)))
for (var in cov){
  if (is.factor(dat_MFP_c_f[,var])) {
    lv = levels(dat_MFP_c_f[,var])
    for(j in lv[-1]){
      list1 = list()
      list1[[var]] = j
      list2 = list()
      list2[[var]] = lv[1]
      diff = setdiff(cov, var)
      for (vri in diff){
        list1[[vri]]=dat_MFP_c_f[1,vri] 
        list2[[vri]]=dat_MFP_c_f[1,vri]
      }
      OR = contrast(mod, list1, list2)
      df = rbind(df, data.frame(list('var'=as.character(var),'level'=as.character(j),'exp.contrast'=as.numeric(exp(OR$Contrast)),
                                     'exp.lower'=as.numeric(exp(OR$Lower)),'exp.upper'=as.numeric(exp(OR$Upper)))))
    }
    
  }else{
    
    list1 = list()
    list1[[var]] = dat_MFP_c_f[1,var]+1
    list2 = list()
    list2[[var]] = dat_MFP_c_f[1,var]
    diff = setdiff(cov, var)
    for (vri in diff){
      list1[[vri]] = dat_MFP_c_f[1,vri]
      list2[[vri]] = dat_MFP_c_f[1,vri]
    }
    OR = contrast(mod, list1,list2)
    
    df=rbind(df,data.frame(list('var'=var, 'level'='NA','exp.contrast'=exp(OR$Contrast),
                                'exp.lower'=exp(OR$Lower),'exp.upper'=exp(OR$Upper))))
  }
}
df = df[-1,]

#raw FP rate
fp_rate = (nrow(subset(dat_MFP_c_f, result==1))/ nrow(dat_MFP_c_f))*100
#raw FP rate of the first screen
fp_rate_first = (nrow(subset(dat_MFP_c_f, round==1 & result==1))/nrow(subset(dat_MFP_c_f, round==1)))*100
#observed FP rate of each round
fp_rate_all = vector(length=M)
for (i in 1:M){
  fp_rate_all[i]= (nrow(subset(dat_MFP_c_f, round==i & result==1))/nrow(subset(dat_MFP_c_f, round==i)))
}

plot(seq(10), fp_rate_all , xlab = ' Screening round',ylab = 'FP probability',ylim = c(0, 0.26) , type="b", col=" deepskyblue3", lty=1, pch=18, font.main=2, xaxt="n")
axis(1, at = seq(0, 10, by = 1), las=2)

#observed FP rate for each round for subjects who attended seven or more rounds
fp_rate_all = vector(length=M)
for (i in 1:M){
  fp_rate_all[i]= (nrow(subset(dat_MFP_c_f, round==i & result==1 & S >=7))/nrow(subset(dat_MFP_c_f, round==i & S>=7)))
}

plot(seq(10), fp_rate_all , xlab = ' Screening round',ylab = 'FP probability',ylim = c(0, 0.26) , type="b", col=" deepskyblue3", lty=1, pch=18, main="Observed false positive probability by screening round", font.main=2, xaxt="n")
axis(1, at = seq(1, 10, by = 1), las=2)

fp_rate_all = vector(length=M)
for (i in 1:M){
  fp_rate_all[i]= (nrow(subset(dat_MFP_c_f, round==i & result==1 & S >=4 & S<=6))/nrow(subset(dat_MFP_c_f, round==i & S>=4 & S<=6)))
}

par(new=T)

plot(seq(10), fp_rate_all , xlab = ' ',ylab = '',ylim = c(0, 0.26) , type="b", col=" deepskyblue3", lty=2, pch=18, font.main=2, xaxt="n")

fp_rate_all = vector(length=M)
for (i in 1:M){
  fp_rate_all[i]= (nrow(subset(dat_MFP_c_f, round==i & result==1 & S >=1 & S<=3))/nrow(subset(dat_MFP_c_f, round==i & S>=1 & S<=3)))
}

par(new=T)

plot(seq(10), fp_rate_all , xlab = ' ',ylab = '',ylim = c(0, 0.26) , type="b", col=" deepskyblue3", lty=3, pch=18, font.main=2, xaxt="n")

par(xpd=TRUE)

legend("bottomleft", legend=c("Women with 7 or more observations", "Women with 4-6 observations", "Women with 1-3 observations"),
       lty=1:3, cex=0.8,col=c("deepskyblue3", "deepskyblue3", "deepskyblue3"),
       text.font=1, bty = "n")

prob_first_fp = vector(length=(M))
for (i in 1:10){
  prob_first_fp[i]= (nrow(subset(dat_MFP_c_f, fp1==i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_first_fp*100

prob_second_fp = vector(length=(M))
for (i in 1:10){
  prob_second_fp[i]= (nrow(subset(dat_MFP_c_f, fp2==i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_second_fp*100

prob_third_fp = vector(length=(M))
for (i in 1:10){
  prob_third_fp[i]= (nrow(subset(dat_MFP_c_f, fp3==i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_third_fp *100

prob_forth_fp = vector(length=(M))
for (i in 1:10){
  prob_forth_fp[i]= (nrow(subset(dat_MFP_c_f, fp4==i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_forth_fp *100

prob_fifth_fp = vector(length=(M))
for (i in 1:10){
  prob_fifth_fp[i]= (nrow(subset(dat_MFP_c_f, fp5==i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_fifth_fp *100

prob_fifth_fp[4]=NA
prob_fifth_fp[3]=NA
prob_fifth_fp[2]=NA
prob_fifth_fp[1]=NA
prob_forth_fp[3]=NA
prob_forth_fp[2]=NA
prob_forth_fp[1]=NA
prob_third_fp[1]=NA
prob_third_fp[2]=NA
prob_second_fp[1]=NA

par (mfrow=c(1,2))

# plot of observed probability of multiple FP results
plot(seq(1,10, by=1), prob_first_fp , xlab = 'Screening round',ylab = 'probability of FP result',ylim = c(0,  0.24588847) , type="b", col=" orange", lty=1, pch=18, font.main=1, xaxt="n", yaxt="n", cex.axis=0.6, main = "Observed probability of multiple FP results", cex.main=0.85)
axis(1, at = seq(1, 10, by = 1), las=2)
axis(2, at=seq(0, 0.25, by=0.05), cex.axis=0.8)

par(new=T)

plot(seq(1,10, by=1), prob_second_fp , xlab = ' ',ylab = '',ylim = c(0,  0.24588847) , type="b", col="deepskyblue3", lty=1, pch=18, font.main=2, xaxt="n", yaxt="n")

par(new=T)

plot(seq(1,10, by=1), prob_third_fp , xlab = ' ',ylab = '',ylim = c(0,  0.24588847) , type="b", col="darkgreen", lty=1, pch=18, font.main=2, xaxt="n", yaxt="n")

par(new=T)

plot(seq(1,10, by=1), prob_forth_fp , xlab = ' ',ylab = '',ylim = c(0,  0.24588847) , type="b", col="black", lty=1, pch=18, font.main=2, xaxt="n", yaxt="n")

par(new=T)

plot(seq(1,10, by=1), prob_fifth_fp , xlab = ' ',ylab = '',ylim = c(0,  0.24588847) , type="b", col="red", lty=1, pch=18, font.main=2, xaxt="n", yaxt="n")

legend("topright", legend=c(expression(paste("P(W"[1],"=j|S≥j)")),expression(paste("P(W"[2],"=j|S≥j)")), expression(paste("P(W"[3],"=j|S≥j)")), expression(paste("P(W"[4],"=j|S≥j)")), expression(paste("P(W"[5],"=j|S≥j)"))), cex=0.62,col=c("orange", "deepskyblue3", "darkgreen", "black","red"), lty=c(1,1), text.font=1, bty = "n")



prob_first_fp_cum = vector(length=M)
for (i in 1:10){
  prob_first_fp_cum[i]= (nrow(subset(dat_MFP_c_f, fp1<=i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_first_fp_cum*100

prob_second_fp_cum = vector(length=M)
for (i in 1:10){
  prob_second_fp_cum[i]= (nrow(subset(dat_MFP_c_f, fp2<=i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_second_fp_cum*100


prob_third_fp_cum = vector(length=M)
for (i in 1:10){
  prob_third_fp_cum[i]= (nrow(subset(dat_MFP_c_f, fp3<=i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_third_fp_cum*100

prob_forth_fp_cum = vector(length=M)
for (i in 1:10){
  prob_forth_fp_cum[i]= (nrow(subset(dat_MFP_c_f, fp4<=i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_forth_fp_cum*100

prob_fifth_fp_cum = vector(length=M)
for (i in 1:10){
  prob_fifth_fp_cum[i]= (nrow(subset(dat_MFP_c_f, fp5<=i &S>=i))/nrow(subset(dat_MFP_c_f, S>=i)))
}
prob_fifth_fp_cum*100

prob_fifth_fp_cum[4]=NA
prob_fifth_fp_cum[3]=NA
prob_fifth_fp_cum[2]=NA
prob_fifth_fp_cum[1]=NA
prob_forth_fp_cum[3]=NA
prob_forth_fp_cum[2]=NA
prob_forth_fp_cum[1]=NA
prob_third_fp_cum[1]=NA
prob_third_fp_cum[2]=NA
prob_second_fp_cum[1]=NA

#plot of empirical cumulative probability of multiple FP results
plot(seq(1,10, by=1), prob_first_fp_cum , xlab = 'Screening round',ylab = 'Cumulative FP probability',ylim = c(0, 0.61) , type="b", col=" orange", lty=1, pch=18, font.main=1, xaxt="n" , yayt="n", cex.axis=0.7, main="Empirical cumulative probability of multiple FP results", cex.main=0.85)
axis(1, at = seq(1, 10, by = 1), las=2)

par(new=T)

plot(seq(1,10, by=1), prob_second_fp_cum , xlab = '',ylab = '',ylim = c(0, 0.61) , type="b", col="deepskyblue3", lty=1, pch=18, font.main=2, xaxt="n", cex.axis=0.7)

par(new=T)

plot(seq(1,10, by=1), prob_third_fp_cum , xlab = '',ylab = '',ylim = c(0, 0.61) , type="b", col="darkgreen", lty=1, pch=18, font.main=2, xaxt="n", cex.axis=0.7)

par(new=T)

plot(seq(1,10, by=1), prob_forth_fp_cum , xlab = '',ylab = '',ylim = c(0, 0.61) , type="b", col="black", lty=1, pch=18, font.main=2, xaxt="n", cex.axis=0.7)

par(new=T)

plot(seq(1,10, by=1), prob_fifth_fp_cum , xlab = '',ylab = '',ylim = c(0, 0.61) , type="b", col="red", lty=1, pch=18, font.main=2, xaxt="n", cex.axis=0.7)

legend("topleft", legend=c(expression(paste("P(W"[1],"≤j|S≥j)")), expression(paste("P(W"[2],"≤j|S≥j)")),expression(paste("P(W"[3],"≤j|S≥j)")), expression(paste("P(W"[4],"≤j|S≥j)")), expression(paste("P(W"[5],"≤j|S≥j)")) ), cex=0.62,col=c("orange", "deepskyblue3", "darkgreen","black","red"), lty=c(1,1), text.font=1, bty = "n")

#############################################################
# In order to compute cumulative risk of one, two, three, four and five false positives 
#for high and low risk group with some baseline characteristics, we simulated
# a dataset using markov chain setting for each set of characteristics.
#these simulated datasets are close to real data and give similar results as real dataset
#for cumulative risk of one, two, three, four and five false positives

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
round_1 = function(i, h,age,fm,fd,rc){
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

#this matrix gives p_ij^(h,t,n_g=x) where 2<=t and t=h, for all i and 
P_prime_nh = function(m40_nh, m41_nh){
  function(h, t, g, x){
    matrix(c(1/(1+exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)+exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)),
             1/(1+exp(m41_nh[2*(g-1)+1,1]+m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)+exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)), 0,
             exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)/(1+exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)
                                                                                                                                 +exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)),exp(m41_nh[2*(g-1)+1,1]+m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)/(1+exp(m41_nh[2*(g-1)+1,1]+
                                                                                                                                                                                                                                                                                                                                                     m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)+exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)), 0,
             exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)/(1+exp(m40_nh[2*(g-1)+1,1]+m40_nh[2*(g-1)+1,2]*x+m40_nh[2*(g-1)+1,3]*h+m40_nh[2*(g-1)+1,4]*t+m40_nh[2*(g-1)+1,5]*h*t)
                                                                                                   +exp(m40_nh[2*g,1]+m40_nh[2*g,2]*x+m40_nh[2*g,3]*h+m40_nh[2*g,4]*t+m40_nh[2*g,5]*h*t)), exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)/(1+exp(m41_nh[2*(g-1)+1,1]+
                                                                                                                                                                                                                                                                                          m41_nh[2*(g-1)+1,2]*x+m41_nh[2*(g-1)+1,3]*h+m41_nh[2*(g-1)+1,4]*t+ m41_nh[2*(g-1)+1,5]*h*t)+exp(m41_nh[2*g,1]+m41_nh[2*g,2]*x+m41_nh[2*g,3]*h+m41_nh[2*g,4]*t +m41_nh[2*g,5]*h*t)), 1), ncol=3)
  }
}

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
  cov.theta <- vector(length = M)
  temp <- NULL
  for (i in 1:M){
    temp = -matrix(theta[i,],ncol =1)%*%matrix(theta[i,],nrow =1)
    diag(temp) = 0
    cov.theta[i] = sum(c(temp))
  }	
  var.theta = apply(sweep(var.theta,1,p[1:M]^2/apply(N,1,sum),"*"),2,sum)-cov.theta/apply(N,1,sum)
  var.theta = sum(var.theta)
}


cens.bias_Markov = function(data, M,alpha, l, functions){
  xround_11 = functions$round_11
  xP_prime2 = functions$P_prime2
  xP2 = functions$P2
  xround_1_nh = functions$round_1_nh
  xP_nh = functions$P_nh
  xP_prime_nh = functions$P_prime_nh
  theta = matrix(data = 0, nrow=M, ncol=M)
  # compute theta[h, j] for lower triangle
  cartesianProd = function(a){
    if(a==1){
      if(l==1){
        return(data.frame("i1"=1))
      } else{ stop("not possible")}
    } else{
      df = expand.grid(rep(list(c(0,1)), a-1))
      names(df) = paste0("i", 1:(a-1))
      df = subset(df, rowSums(df)==l-1)
      df[, paste0("i",a )] = 1
      return(df)}}
  # j = 1, h= 1:M creates the first column
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
              return (xP_prime2(h, k+1)[x[k]+1, x[k+1]+1])
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
        B[[h]][b, c+1] = nrow(subset(data, S==b & get(paste0("n",h))==c))/nrow(data)
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
                    return (xP_prime_nh(b, k+1, h, c)[x[k]+1, x[k+1]+1])
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
  data$RESINIT_C = ifelse(data[,paste0("fp",l)]==data$round, 1,0)
  N1 = table(data$S[data$RESINIT_C==1], data$round[data$RESINIT_C==1])
  N2 = matrix(0, ncol=M, nrow=M)
  N2[1:nrow(N1), 1:ncol(N1)]=N1
  N = table(data$S,data$round)
  varp = cum.se(theta, M, alpha, N, N2, p1)
  se = sqrt(varp)
  return (list(risk=risk, se = se, theta=theta))
}


#discrete survival model 
disc.surv = function(data,M,l){
  theta_hat = vector(length = M)
  r = vector (length = M)
  s = vector (length =  M)
  varp = vector (length =M)
  for (k in 1:M){
    theta_hat[k] = mean(data[,paste0("fp",l,".cens")][data$S>=k & (data[,paste0("fp",l,".cens")]>=k |data[,paste0("fp",l,".cens")]==0)] == k)
    r[k] =  length(data[,paste0("fp",l,".cens")][data$S>=k & (data[,paste0("fp",l,".cens")]>=k | data[,paste0("fp",l,".cens")]==0) ])
    s[k]= sum(data[,paste0("fp",l,".cens")][data$S>=k & (data[,paste0("fp",l,".cens")]>=k | data[,paste0("fp",l,".cens")]==0) ]>k|data[,paste0("fp",l,".cens")][data$S>=k & (data[,paste0("fp",l,".cens")]>=k | data[,paste0("fp",l,".cens")]==0) ]==0)
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
pop.ave.var = function(S, cens, comp, delta, M){
  sig = vector(length = M)
  N = max(S)*(max(S)+3)/2-1
  # sigma_1
  eta = vector(length = N)
  k = 1
  for (i in 1:max(S)){
    for (j in i:max(S)){
      eta[k] = mean(S == j & cens == i & delta == 1 & comp>=i)
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
      mu[[i]][j-i+1] = mean(S == j & cens == i & delta == 1 & comp>=i)
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
  varp = pop.ave.var(data$S,data[,paste0("fp",l,".cens")], data$comp, data$delta, M)
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

#this function is used to simulate a data similar to real data for each set of women's characteristics 
#function for cumulative risk of >=2 false positives
MfpRisk_sim_real = function(Nsim, M, p, l, Nsubj,age,fm,fd,rc,int){
  if (int==2){M=M/2}
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  #simulating censoring time from Pr 
  out1 = matrix (ncol=3, nrow= Nsim) 
  se1 = matrix (ncol=3, nrow= Nsim)
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
  }
  else{
    
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
  
  #add the variable for the time of the l-th false positive 
  data3 = ddply(.data = data2, .variable ="StudyID_c", .fun = function(df){
  df[,paste0("fp",l)] = ifelse(is.na(which(df$result==1)[l]), M+1, which(df$result==1)[l])
    return(df)
  })
  
  #add M variables ni when ni is the number of false positive by and including round i for i=1,...,M
  data4 = ddply(.data= data3, .variables = "StudyID_c", .fun = function(df) {
    for (i in 1:M){
      df[, paste0("n",i)] = sum(df$result[1:i]==1)
    }
    return(df)
  })
  
  #add the variable for the time of cancer (when result is 2) (competing event)
  data5 = ddply(.data = data4, .variable ="StudyID_c", .fun = function(df){
    df$comp = ifelse(is.na(which(df$result==2)[1]), M+1, which(df$result==2)[1])
    return(df)
  })
  
  #add variable for the censoring if l-th fp 
  data5[,paste0("cens",l)] = apply(cbind(data5$S,data5[,paste0("fp",l)]),1,min)  
  
  #define event time for the other two models
  data5[,paste0("fp",l,".cens")] = ifelse(data5[,paste0("fp",l,".cens")]>data5$S,0,data5[,paste0("fp",l,".cens")]) 
  
  #define censoring indicator for the other two models
  data5$delta = ifelse(data5[,paste0("fp",l)] <= data5$S, 1, 0)
  
  set1 = subset(data5,  1==round & round ==S)
  set2 = subset(data5,  1==round & round <= (S-1))
  set3 = subset(data5,  2<=round & round <= (S-1))
  set4 = subset(data5,  2<=round & round ==S)
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
    mod30 = glm(result~set30[,paste0("n",i)]+S+round+S*round , data = set30, family=binomial)
    m30_nh[i,] = mod30$coefficients
  }
  m30_nh[is.na(m30_nh)]=0
  
  m31_nh = matrix(0,nrow = M, ncol=5) 
  for (i in 1:M){
    mod31 = glm(result~set31[,paste0("n",i)]+ S+round+S*round , data =set31 , family=binomial)
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
  
  g = data6[data6$cens3 <= data6$S,]   
  g.dat <- NULL
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

#two FPs
MfpRisk_sim_real(5000,10, 0.2, 2, 50000, 40,1,3,2,1)#high risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 2, 50000, 40,1,3,2,2)#high risk, biennial screener
MfpRisk_sim_real(5000,10, 0.2, 2, 50000, 50,0,1,3,1) #low risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 2, 50000, 50,0,1,3,2) #low risk, biennial screeners

#Three FPs
MfpRisk_sim_real(5000,10, 0.2, 3, 50000, 40,1,3,2,1)#high risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 3, 50000, 40,1,3,2,2)#high risk, biennial screener
MfpRisk_sim_real(5000,10, 0.2, 3, 50000, 50,0,1,3,1) #low risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 3, 50000, 50,0,1,3,2) #low risk, biennial screeners

#four FPs
MfpRisk_sim_real(5000,10, 0.2, 4, 50000, 40,1,3,2,1)#high risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 4, 50000, 40,1,3,2,2)#high risk, biennial screener
MfpRisk_sim_real(5000,10, 0.2, 4, 50000, 50,0,1,3,1) #low risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 4, 50000, 50,0,1,3,2) #low risk, biennial screeners

#Five FPs
MfpRisk_sim_real(5000,10, 0.2, 5, 50000, 40,1,3,2,1)#high risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 5, 50000, 40,1,3,2,2)#high risk, biennial screener
MfpRisk_sim_real(5000,10, 0.2, 5, 50000, 50,0,1,3,1) #low risk, annual screener
MfpRisk_sim_real(5000,10, 0.2, 5, 50000, 50,0,1,3,2) #low risk, biennial screeners


#For one False positive the following code works faster
f = function( X,  A, alpha){
  M = nrow(A)
  for (j in (M-1):1)
    for(i in (j+1):M){
      Y = 1 - apply(A[,(j+1):M,drop=F],2,sum)
      A[i,j] = exp(alpha*i) * sum(X[(j+1):length(X)] * A[i, (j+1):M]) / exp(alpha*seq(M+1)) %*% rbind(A[,(j+1):M,drop=F], Y) %*% as.matrix(X[(j+1):M])
      
    }
  A
}
cens.bias_c = function(data,M,alpha){
  data$RESINIT_C = ifelse(data$fp1==data$round, 1,0)#event is the first FP
  comp=data$comp
  t=summary(factor(comp,levels=c(1:(M))))
  cr  = vector(length=M)
  cr[1]=1
  for(i in 2:(M)){
    cr[i]= (sum( table(comp))-cumsum(t)[i-1][[1]])/sum(table(comp))
  }
  cr =ifelse(is.na(cr), 0, cr)
  data0 =data[data$round <= data$cens1,]
  data1=data0
  round.main = paste("round==",seq(1,M))
  S.main = paste("S==",seq(1,M))
  rS.int = vector(length = (M-1)*(M-2)/2)
  l <- 1
  for (i in 2:(M-1)){
    for (j in 2:i){
      rS.int[l] <- paste("(S==",i,")*(round==",j,")",sep="")
      l <- l + 1
    }
  }       
  xnam = paste(rS.int,collapse="+")
  
  fm = formula(paste("RESINIT_C ~ (S == ",M,") +","(round == ",M,") +",xnam,sep=""))
  
  # full data estimate upper triangular, rows represent S and cols represent round
  fp.glm = glm(fm, data = data1, family = "binomial")
  theta = matrix(ncol = M, nrow = M)
  for (i in 1:M){
    for (j in 1:i){
      theta[i,j] = fitted(fp.glm)[data1$S ==i & data1$round == j & data1$comp>=j][1]            
    }
  }
  p = summary(factor(data0$S[!duplicated(data0$StudyID_c)],levels=c(1:M)))/length(data0$S[!duplicated(data0$StudyID_c)]) 
  
  thetastar = matrix(0,ncol = M, nrow = M)
  # create thetastar from theta
  thetastar[,1] = theta[,1]
  for (i in 2:M){#S
    for (j in 2:i)#W
      thetastar[i,j] = theta[i,j]*prod(1-theta[i,(1:(j-1))])} #transform from conditional to joint probs
  thetastar = ifelse(is.na(thetastar),0,thetastar)
  v = f(p,t(thetastar),alpha)
  h=cr%*%v
  risk=h%*%p
  return (risk_cb=risk)
}

disc.surv.c = function(data,M){ 
  S=data$S
  W= data$fp1.cens
  comp=data$comp
  t=summary(factor(comp,levels=c(1:M)))
  cr  = vector(length=(M))
  cr[1]=1
  for(i in 2:(M)){
    cr[i]= (sum(t)-cumsum(t)[i-1][[1]])/sum(t)
  }
  cr =ifelse(is.na(cr), 0, cr)
  W.cens = ifelse(W>S,0, W) 
  theta_hat = vector(length = (M))
  for (k in 1:(M)){
    theta_hat[k] = mean(W.cens[S>=k & (W.cens>=k | W.cens==0) & comp>=k] == k)
  }
  theta_hat = ifelse(is.na(theta_hat),0, theta_hat)
  # create joint probs from conditional probs
  thetastar = vector(length = (M))
  thetastar[1] <- theta_hat[1]
  for (j in 2:M)
  thetastar[j] <- theta_hat[j]*prod(1-theta_hat[(1:(j-1))]) 
  thetastar <- thetastar*cr
  risk = cumsum(thetastar)[M]
  return(risk_ds = risk)
}

#population average model
pop_ave_c = function(data,M){
  comp=data$comp
  S=data$S
  W=data$fp1.cens
  t=summary(factor(data$comp,levels=c(1:M)))
  cr  = vector(length=(M))
  cr[1]=1
  for(i in 2:(M)){
    cr[i]= (sum( table(data$comp))-cumsum(t)[i-1][[1]])/sum( table(data$comp))
  }
  cr =ifelse(is.na(cr), 0, cr)
  
  s_hat = vector(length = M)
  p_hat = matrix(nrow = M, ncol = M)
  theta_hat = vector(length = M)
  
  # j (rows) indexes number screens, k (cols) indexes event time
  for (j in 1:M){
    s_hat[j] = mean(S == j)
    theta_hat[j] = 1-mean(W[S==j] == 0)^(1/j) 
    for (k in 1:M){
      p_hat[j,k] = mean(W[S==j & comp>=k] == k) 
    }
  }
  
  theta_hat = ifelse(is.na(theta_hat),0,theta_hat)
  p_hat = ifelse(is.na(p_hat),0,p_hat)
  
  # j indexes number of screens, k indexes event time
  xi_hat <- vector(length = M)
  for (k in 1:M){
    if (k == 1){
      xi_hat[k] = 0	
      for (j in 1:M)
        xi_hat[k] = xi_hat[k] + p_hat[j,k]*s_hat[j]*cr[k] # everyone has at least one screen	
    }
    if (k > 1 & k <= M){
      xi_hat[k] = 0	
      for (j in 1:(k-1)) # contribution for subjects with fewer than k screens (s < k)
        xi_hat[k] = xi_hat[k]+theta_hat[j]*(1-theta_hat[j])^(k-1)*s_hat[j]*cr[k]
      for (j in k:M) # contribution for subjects with at least k screens
        xi_hat[k] = xi_hat[k]+p_hat[j,k]*s_hat[j]*cr[k]
    }
    if (k > M){
      xi_hat[k] = 0
      for (j in 1:M) # no one has more than M screens
        xi_hat[k] = xi_hat[k]+theta_hat[j]*(1-theta_hat[j])^(k-1)*s_hat[j]*cr[k]		
    }			
  }	
  phat_ave = cumsum(xi_hat)
  return(risk_pa= phat_ave[M])
}

#function for cumulative risk of one false positive
MfpRisk_sim_real_one = function(Nsim,M, p, l, Nsubj,age,fm,fd,rc,int){
  if (int==2){M=M/2}
  #geometric ditribution with prob=p for censoring time
  P_s = p*((1-p)^seq(0,(M-1)))
  #sum of geometric distribution for censoring time should be 1
  Pr = c(P_s[1:(M-1)],1-sum(P_s[1:(M-1)]))
  out1 = matrix (ncol=3, nrow= Nsim) 
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
  #change this so it depends on the number of prior FPs
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
  }
  else{
    
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
  
  #add the variable for the round of the first fp
  data3 = ddply(.data = data1, .variable ="StudyID_c", .fun = function(df){
    df$fp1 = ifelse(is.na(which(df$result==1)[1]), M+1, which(df$result==1)[1])
    return(df)
  })
  
  #add the variable for the time of cancer (when result is 2)
  data6 = ddply(.data = data3, .variable ="StudyID_c", .fun = function(df){
    df$comp = ifelse(is.na(which(df$result==2)[1]), M+1, which(df$result==2)[1])
    return(df)
  })
  
  #add variable for the censoring if first fp is the event (min(S,fp1))
  data6$cens1 = apply(cbind(data6$S,data6$fp1),1,min)  
  
  
  #define event time for the other two models
  data6$fp1.cens = ifelse(data6$fp1>data6$S,0,data6$fp1) 
  
  g = data6[data6$fp1 <= data6$S,] 
  g.dat <- NULL
  for (j in 2:(M-1)) 
    g.dat = tryCatch(rbind(g.dat,cbind( g[g$round <= j & g$S >=j  ,],j)), error = function(e){g.dat})
  
  g.dat = tryCatch(g.dat[g.dat$round <= g.dat$j ,],error=function(e){g.dat})
  
  g.dat$alpha =(g.dat$round - (g.dat$j+1))
  
  G = glm((S == j) ~ alpha, family = "binomial", data = g.dat)$coef[2]
  
  cum_ds = disc.surv.c(data6,M)
  cum_pa = pop_ave_c (data6,M)
  cum_cb = cens.bias_c (data6,M,G)
  
  out1[jj,1] = cum_ds$risk_ds
  out1[jj,2] = cum_pa$risk_pa
  out1[jj,3] = cum_cb$risk

  }
  return(list(cumr = apply(out1, 2, mean), estr = apply(out1, 2, sd), qtls = apply(out1, 2, quantile, probs = c(0,0.25, 0.75, 0.975, 1),  na.rm = TRUE)))
}

#one FP
MfpRisk_sim_real_one(5000,10, 0.2, 1, 50000, 40,1,3,2,1)#high risk, annual screener
MfpRisk_sim_real_one(5000,10, 0.2, 1, 50000, 40,1,3,2,2)#high risk, biennial screener
MfpRisk_sim_real_one(5000,10, 0.2, 1, 50000, 50,0,1,3,1) #low risk, annual screener
MfpRisk_sim_real_one(5000,10, 0.2, 1, 50000, 50,0,1,3,2) #low risk, biennial screeners

