library(tseries)
library(lmtest)
# library(TTR)
library(forecast)
library(plyr)
library(dplyr)
library(ggplot2)
library(tseries)

setwd("C:/Users/me/Downloads")

df <- read.csv("NYPD_Motor_Vehicle_Collisions.csv",stringsAsFactors = F)
head(df)
dim(df)

str(df)

df$date <- as.Date(df$DATE,"%m/%d/%Y")

dfs <- df %>% 
  group_by(date) %>%
  summarize(cnt=n())

head(dfs)
dim(dfs)


ggplot(dfs,aes(date,cnt)) + geom_line() + geom_smooth()


write.csv(dfs,file="NYPD_MVC.csv", row.names = F)


