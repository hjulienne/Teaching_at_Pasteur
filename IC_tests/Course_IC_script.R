





#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@     Optional exercices
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

############### Q1 ######################

## Extract from the entire table
temp_ART_week1  = data$haem.hb.level[data$pv.lm=="yes" & data$group=="ART" & data$week==1]
temp_ART_week12 = data$haem.hb.level[data$pv.lm=="yes" & data$group=="ART" & data$week==12]

## Get the quantile
u = qnorm(0.95)

## Parameter estimation
mean.week1 = mean(temp_ART_week1);  mean.week12 = mean(temp_ART_week12)
sd.week1 = sd(temp_ART_week1);      sd.week12 = sd(temp_ART_week12)
n.week1 = length(temp_ART_week1);   n.week12 = length(temp_ART_week12)

## Confidence interval (95%)
ICsup.week1 = mean.week1 + u*sd.week1/sqrt(n.week1);  ICsup.week12 = mean.week12 + u*sd.week12/sqrt(n.week12)
ICinf.week1 = mean.week1 - u*sd.week1/sqrt(n.week1);  ICinf.week12 = mean.week12 - u*sd.week12/sqrt(n.week12)

cat(paste("95% confidence interval = [",round(ICinf.week1,2),";",round(ICsup.week1,2),"]"))
cat(paste("95% confidence interval = [",round(ICinf.week12,2),";",round(ICsup.week12,2),"]"))

############### Q2 ######################

## Extract from the entire table
temp_ART_week1  = data$haem.hb.level[data$pv.lm=="yes" & data$group=="ART" & data$week==1 & data$village=="UTAMUP"]
temp_ART_week12 = data$haem.hb.level[data$pv.lm=="yes" & data$group=="ART" & data$week==12  & data$village=="UTAMUP"]

## Get the size
n1 = length(temp_ART_week1)
n12 = length(temp_ART_week12)

## Get the quantile
t1 = qt(0.975,n1-1)
t12 = qt(0.975,n12-1)

## Parameter estimation
mean.week1 = mean(temp_ART_week1);  mean.week12 = mean(temp_ART_week12)
sd.week1 = sd(temp_ART_week1);      sd.week12 = sd(temp_ART_week12)

## Confidence interval (95%)
ICsup.week1 = mean.week1 + t1*sd.week1/sqrt(n1);  ICsup.week12 = mean.week12 + t12*sd.week12/sqrt(n12)
ICinf.week1 = mean.week1 - t1*sd.week1/sqrt(n1);  ICinf.week12 = mean.week12 - t12*sd.week12/sqrt(n12)

cat(paste("95% confidence interval = [",round(ICinf.week1,2),";",round(ICsup.week1,2),"]"))
cat(paste("95% confidence interval = [",round(ICinf.week12,2),";",round(ICsup.week12,2),"]"))



############### Proportions ######################

## Extract from the entire table
tmp  = data$village[data$pv.lm=="yes" & data$group=="ART&PQ" & data$week==12]

## Get the size
n = length(tmp)

## Get the quantile
u = qnorm(0.95)

## Parameter estimation
prop = sum(tmp=="UTAMUP")/n
sd.prop = sqrt(prop*(1-prop)/n)

## Confidence interval (90%)
ICsup = prop + u*sd.prop; ICsupBin = qbinom(0.95,n,prop)/n
ICinf = prop - u*sd.prop; ICinfBin = qbinom(0.05,n,prop)/n

cat(paste("95% confidence interval = [",round(ICinf,2),";",round(ICsup,2),"]"))
cat(paste("95% confidence interval = [",round(ICinfBin,2),";",round(ICsupBin,2),"]"))






#################################################################################################           
#################################           Tests               #################################           
#################################################################################################  

###########################################
###
###       Comparaison of means
###
###########################################

###################### EXERCISES ON LUNG DATASET
## Get the data
data = read.csv("lungA.csv")
xmale = data$tumor_size[data$gender=="male"]; xfemale = data$tumor_size[data$gender=="female"]

## Parameter estimation
mu.male = mean(xmale); mu.female = mean(xfemale)
sigma2.male = var(xmale); sigma2.female = var(xfemale)
n.male = length(xmale); n.female = length(xfemale)

## Test statistic
S = (mu.male - mu.female)/sqrt(sigma2.male/(n.male) + sigma2.female/(n.female))

## Critical value
cv = qnorm(0.025)

## p-value (Normal approx)
pv.norm = 2*pnorm(-S)

## p-value (Welch test)
pv.welch = 2*pt(-S,df=83.76)


## Run t.test
res.test1 = t.test(xmale,xfemale)

## Alternative 
res.test2 = t.test(data$tumor_size~data$gender)

###################### EXERCISES ON MALARIA DATASET
## Get the data
data <- read.csv("malaria_longitudinal_data_simul.csv")

## Compare mean temperature of sick and not sick W1
temp_sick_w1 <- data$temp[data$week == 1 & data$pv.lm == "yes"]
temp_nosick_w1 <- data$temp[data$week == 1 & data$pv.lm == "no"]
t.test(temp_sick_w1, temp_nosick_w1)

## Compare mean temperature of sick and not sick W1
haem_control_sick_w12 <- data$haem.hb.level[data$week == 12 & data$pv.lm == "yes" & data$group == "Control"]
haem_control_nosick_w12 <- data$haem.hb.level[data$week == 12 & data$pv.lm == "no" & data$group == "Control"]
t.test(haem_control_sick_w12, haem_control_nosick_w12)

## Compare proportion of PCR-sick in control group and ART&PQ group at week 12
sick_control_w12 <- data$pv.qpcr[data$week == 12 & data$group == "Control"]
sick_artpq_w12 <- data$pv.qpcr[data$week == 12 & data$group == "ART&PQ"]
n_control_w12 <- length(sick_control_w12)
n_artpq_w12 <- length(sick_artpq_w12)
prop.test(x = c(sum(sick_control_w12 == "yes"), sum(sick_artpq_w12 == "yes", na.rm = TRUE)), n = c(n_control_w12, n_artpq_w12))
