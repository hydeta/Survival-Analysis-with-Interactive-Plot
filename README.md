# Survival Analysis with Interactive Plotting



## Background

Survival Analysis is a technique most frequently used in the biological sciences to predict the probability of an individual surviving past time T > t given that they have survived up until time t already.  One of its strengths is the ability to make use of censored data.  In a logistic regression, if an event is not observed within the scope of the study it is considered to not have happened.  Survival analysis is more nuanced such that if an event is not observed within the scope of the study it is considered to not have happened **yet**.  The probability that a subject will survive past time t is given by: $$S(t)=Pr(T>t)=1-F(t)$$ which is the inverse of the cumulative distribution function.

## Simple Survival Regression - Kidney Infection

Let's start with a trivial example using data on catheter placement to determine the probability of developing a kidney infection.  Catheters can be placed either surgically or percutaneously (i.e. through the skin).  The dataset "kidney" from the KMsurv package contains three variables:

1. `time` time until infection, in months
2. `delta` infection indicator (1 = yes, 0 = no)
3. `type` the type of catheter placement (1 = surgical, 2 = percutaneous)


```r
library(dplyr)
library(anytime)
library(survival)
library(plotly)
library(viridis)
library(survrec)
library(KMsurv)
library(rms)
data("kidney")
head(kidney)
```

```
##   time delta type
## 1  1.5     1    1
## 2  3.5     1    1
## 3  4.5     1    1
## 4  4.5     1    1
## 5  5.5     1    1
## 6  8.5     1    1
```



First let's filter to create two separate dataframes based on the type of placement, then determine the survival function using a Kaplan-Meier estimator.


```r
#filter for each user type
surgical <- kidney %>% filter(type==1)
percut <- kidney %>% filter(type==2)

#create Kaplan-Meier estimator.
km.surgical <- survfit(Surv(time,delta)~1, data=surgical, conf.type="log-log")
km.percut <- survfit(Surv(time,delta)~1, data=percut, conf.type="log-log")

summary(km.surgical)
```

```
## Call: survfit(formula = Surv(time, delta) ~ 1, data = surgical, conf.type = "log-log")
## 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##   1.5     43       1    0.977  0.0230       0.8462        0.997
##   3.5     40       1    0.952  0.0329       0.8224        0.988
##   4.5     36       2    0.899  0.0478       0.7532        0.961
##   5.5     33       1    0.872  0.0536       0.7190        0.945
##   8.5     25       2    0.802  0.0683       0.6250        0.902
##   9.5     22       1    0.766  0.0743       0.5803        0.877
##  10.5     20       1    0.728  0.0799       0.5350        0.851
##  11.5     18       1    0.687  0.0851       0.4886        0.822
##  15.5     11       1    0.625  0.0976       0.4058        0.782
##  16.5     10       1    0.562  0.1060       0.3350        0.738
##  18.5      9       1    0.500  0.1111       0.2726        0.691
##  23.5      4       1    0.375  0.1366       0.1311        0.623
##  26.5      2       1    0.187  0.1491       0.0143        0.517
```

```r
summary(km.percut)
```

```
## Call: survfit(formula = Surv(time, delta) ~ 1, data = percut, conf.type = "log-log")
## 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##   0.5     76       6    0.921  0.0309        0.833        0.964
##   2.5     56       2    0.888  0.0376        0.788        0.943
##   3.5     49       1    0.870  0.0409        0.763        0.931
##   6.5     35       1    0.845  0.0467        0.726        0.915
##  15.5     14       1    0.785  0.0726        0.599        0.892
```

The column `survival` above gives the probability of surviving at each time t.  It appears as though surgical catheter placement carries with it a significantly greater risk of infection than percutaneous placement.  It is necessarily a step function because the event of interest is always discrete.  However, it is interesting to note that the majority of infections for percutaneous placements occured within the first month.  Next we can plot the two survival curves:

```r
km.kidney <- npsurv(Surv(time,delta)~type, data=kidney)
```
If we now choose to plot both curves on the same chart using plotly, we can create an interactive plot that allows us to see the probability of survival at each time T > t:
<iframe src="https://plot.ly/~thyde/31.embed" width="800" height="600" id="igraph" scrolling="no" seamless="seamless" frameBorder="0"> </iframe>

Which shows that indeed the risk of infection from percutaneous catheter placement is greater than from surgical placement, but only for the first 8.5 months, after which surgical placements have a much greater risk of infection.

## Survival Analysis in a Business Context

In the context of transactional data, one can use survival analysis to predict when a customer will purchase next.  In this case, we redefine our 'death' event to be a purchase/order:  


```r
transact <- read.table("CDNOW_sample.txt", 
               sep="",
               col.names=c("joinID", "userID", "date", "count", "total"),
               fill=FALSE, 
               strip.white=TRUE)
transact$date <- anydate(transact$date)

transact <- transact %>%
    group_by(joinID,userID) %>%
    arrange(date) %>%
    mutate(delta = date - lag(date, default=first(date)),cumCount=cumsum(count),cumTotal=cumsum(total))

#now we have a delta in days and cumulative sums:
head(transact %>% filter(userID==129))
```

```
## Source: local data frame [5 x 8]
## Groups: joinID, userID [1]
## 
## # A tibble: 5 x 8
##   joinID userID       date count total   delta cumCount cumTotal
##    <int>  <int>     <date> <int> <dbl>  <time>    <int>    <dbl>
## 1   1544    129 1997-01-07     1  6.79  0 days        1     6.79
## 2   1544    129 1997-01-09     2  9.58  2 days        3    16.37
## 3   1544    129 1997-01-24     2 19.16 15 days        5    35.53
## 4   1544    129 1997-02-13     1 13.97 20 days        6    49.50
## 5   1544    129 1997-03-01     1 11.77 16 days        7    61.27
```

Since each user makes multiple purchases over time, we use a variation on survival analysis that accounts for recurrent events.  For recurrent events, there are generally two approaches to time:

1. Calendar Time - all events are measured against an initial starting date.
2. Gap Time - time between events.

Each approach provides different information; the calendar time approach gives access to the distribution of events within a certain time window, whereas gap time 
gives access to the distribution of gaps.  Calendar time is often most useful for incident events where the process itself is not altered.  Conversely, gap time is stronger when predicting time until next event, and the assumption is that a user is 'reset' between events.  Finally, since we are interested in the time between events, the last event for each user is always considered a censored observation because we do not know the time delta going forward.  

Using the package `survrec` to build our reccurent event survival object and regressing:

```r
#add event column to identify censored observations
event<-function(dataframe){
  dataframe$uniqueID<-rownames(dataframe)
  last_item<-by(dataframe,dataframe$userID,tail,n=1)
  last_df<-do.call("rbind",as.list(last_item))
  dataframe$event<-1-as.integer(dataframe$uniqueID %in% last_df$uniqueID)
  return(dataframe)
}

regress<-event(transact)
regress<-subset(regress,select=-joinID)
rec.fit <- survfitr(Survr(userID,delta,event)~1,data=regress)
```

```
## 
## Needs to Determine a Seed Value for Alpha
##  Seed Alpha:  4.626288
##  
##  Alpha estimate= 0.777785
## 
```

```r
plot(rec.fit,ylim=c(0,1),conf.int=TRUE)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

The result is a slightly more nuanced survival function that recognizes when certain users make a repeat purchase and incorporates this information into the model.

**Takeaways:** Survival analysis can be an extremely powerful tool for making use of time-to-event data.   For example, other regression methods ignore censored data.  For example, if I am looking to see whether a customer purchased or not, a linear regression would treat "0" as "this customer did not purchase" where survival analysis will treat the observation as "this customer did not purchase yet," which is a subtle but important distinction.  One could also use a logistic regression to compare proportion of events between groups, but this ignores the time component.  Survival analysis is unique in its ability to make use of both censored observations and time-to-event data.

**Generalization:** At Wayfair, I built a survival regression to predict the probability of user clicking on a display ad given an impression: $$Pr(Click|Impression)$$ This is important because online advertising is sold through a real-time second price auction, so if I am able to determine the likelihood of a user clicking on an impression I can scale that amount I am willing to pay for that auction.  In this situation I employed a Cox proportional hazard model because its main assumption is that covariates are multiplicatively related to the hazard i.e. the covariates only create a monotonic transformation of the underlying survival function, which remains constant.  This is a strong assumption, but by limiting my population to users who had previously clicked on an ad, I was able to create a model that increased savings by 10%.
