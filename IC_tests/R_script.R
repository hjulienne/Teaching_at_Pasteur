
########################################################
####
####      Confidence intervals (Bootstrapping)
####
########################################################

## data
x = c(43, 55, 45, 46, 57, 52, 61, 47, 52, 57, 48, 57, 52, 62, 50, 53, 72, 58, 58, 45)
hist(x,breaks=3)

## Number of resampling
B = 5000

## Resampling
estim = c()
for(i in 1:B){
  subsamp = sample(x,size = 20,replace = TRUE)
  print(subsamp)
  estim = c(estim, mean(subsamp))    
}
## Representation (histogram)
hist(estim,main = paste("Mean distribution (B=",B,")"))
mean.boot = mean(estim)

## Get the CI
ICinf = quantile(estim,probs = 0.025)
ICsup = quantile(estim,probs = 0.975)
cat(paste("95% confidence interval = [",round(ICinf,2),";",round(ICsup,2),"]"))

## Representation (add IC on the hist)
abline(v=ICinf,col="red",lwd=4)
abline(v=ICsup,col="red",lwd=4)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@     Bootstrap
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Get the data
x = data$haem.hb.level[data$village=="UTAMUP" & data$group=="ART" & data$week==1]

## Number of resampling
B = 15000

## Resampling
estim = c()
for(i in 1:B)
{
  subsamp = sample(x,size = length(x),replace = TRUE)
  estim = c(estim, mean(subsamp))    
}
## Representation (histogram)
hist(estim,main = paste("Mean distribution (B =", B,")"))

## Get the CI
ICinf = quantile(estim,probs = 0.025)
ICsup = quantile(estim,probs = 0.975)
cat(paste("95% confidence interval = [",round(ICinf,2),";",round(ICsup,2),"]"))

## Representation (add IC on the hist)
abline(v=ICinf,col="red",lwd=4)
abline(v=ICsup,col="red",lwd=4)


################################################
####
####    IC for a proportion
####
################################################

#########################################
## 95% IC of the smoker proportion
#########################################

data = read.csv("lungA.csv")

## Get the data
x = data$smoker

## Get the size
n = length(x)

## Get the quantile
u = qnorm(0.975)

## Parameter estimation
prop = sum(x)/n
sd.prop = sqrt(prop*(1-prop)/n)

## Confidence interval (95%)
ICsup = prop + u*sd.prop
ICinf = prop - u*sd.prop
cat(paste("95% confidence interval = [",round(ICinf,2),";",round(ICsup,2),"]"))


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@     Exercices (Ic of the proportion, hugesmall samples)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#########################################
## 95% IC proportion of patients ARTPQ
#########################################
## Load the data
data = read.csv("malaria_longitudinal_data_simul.csv")

## Get the data
x1 = data$pv.lm[data$group=="ART&PQ" & data$week==1]
x12 = data$pv.lm[data$group=="ART&PQ" & data$week==12]

## Get the size
n1 = length(x1)
n12 = length(x12)

## Get the quantile
u = qnorm(0.975)

## Parameter estimation
prop1 = sum(x1=="yes")/n1; 
prop12 = sum(x12=="yes")/n12
sd.prop1 = sqrt(prop1*(1-prop1)/n1); 
sd.prop12 = sqrt(prop12*(1-prop12)/n12)

## Confidence interval (90%)
ICsup1 = prop1 + u*sd.prop1; 
ICsup12 = prop12 + u*sd.prop12
ICinf1 = prop1 - u*sd.prop1; 
ICinf12 = prop12 - u*sd.prop12
cat(paste("95% confidence interval = [",round(ICinf1,2),";",round(ICsup1,2),"]"))
cat(paste("95% confidence interval = [",round(ICinf12,2),";",round(ICsup12,2),"]"))


#######################################################################
## 95% IC proportion of patients ARTPQ (without the gaussian approximation)
#######################################################################

## Get the data
x = data$pv.lm[data$group=="ART&PQ" & data$week==12]

## Get the size
n = length(x)

## Parameter estimation
prop = sum(x=="yes")/n

## Confidence interval (95%)
ICsup = qbinom(0.975,n,prop)/n
ICinf = qbinom(0.025,n,prop)/n

cat(paste("95% confidence interval = [",round(ICinf,2),";",round(ICsup,2),"]"))





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






## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
  
  ```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:
  
  ```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
