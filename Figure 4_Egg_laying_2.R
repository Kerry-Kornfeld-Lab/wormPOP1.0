# Intro -------------------------------------------------------------------
# clears workspace
rm(list = ls())


#load packages
library(ggplot2)
library(multcomp)

#load file
egg_total_data=read.csv()
egg_total_data$culling = as.factor(egg_total_data$culling) #changedata in column culling into factor
attach(egg_total_data) #attach makes datafram accessible



#boxplot
ggplot(data=egg_total_data, aes(x=culling, y=total_egg)) + 
  geom_boxplot(lwd=1)+
  theme_classic()+
  scale_y_continuous(breaks = pretty(c(0,120), n = 5),
                     limits = c(0,120))+
  geom_jitter(position=position_jitter(0.3))

aggregate(total_egg ~ culling, 
  data = egg_total_data,
  function(x) round(c(mean (x), sd(x)), 2)
  )

#do anova on the data
aov.culling = aov(total_egg~as.factor(culling))
summary(aov.culling)
ls(aov.culling)

#Tukey HSD test:
#summary(glht(aov.culling, linfct = mcp(group = "Tukey")))

#summary(post_test)
TukeyHSD(aov.culling)
