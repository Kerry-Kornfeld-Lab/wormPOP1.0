# Intro -------------------------------------------------------------------
# clears workspace
rm(list = ls())
#save the start time so that total run time can be calculated at the end, just for fun :)
start.time <- Sys.time() 
#install and load all packages and libraries
#the install lines only need to be run the first time you do this
# install.packages("ggplot2")
library(psych)
library(dplyr)
library(tidyr) #needed to arrange dataframes with gather() and spread()
library(FSA) #needed for dunnTest()
library(gridExtra)
library(survival)
library(ggplot2)
library(ggfortify)

#this sets the font size to be bigger by default
# theme_set(theme_gray(base_size = 24))
theme_set(theme_classic(base_size = 14))
WD=5    #size for smallest panels
HT=3  #size for smallest panels


#a more succint function than summary()
meansdN <- function(x){c(mean=mean(x), sd=sd(x), N=length(x))} 



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Lifespan ----------------------------------------------------------------
#just enter the name of the lifespan file, without the.csv
file.name<-"Lifespan"
input.file.name<-paste(file.name,".csv", sep="")

#EXAMPLE OF FORMAT:
# Genotype	Letter	Day	Alive	Died	Bagged	Censored
# N2	      N	      1	  34	  0	    0	      0
# N2	      N	      2	  33	  0	    0	      1
# N2	      N	      3	  33	  0	    0	      0
# N2	      N	      4	  33	  0	    0	      0
# N2	      N	      5	  31	  0	    1	      1

#type in correct order of graphs
CANDIDATES <- c("2","1","0.5","s1") #used for pairwise graphing and ordering within graphs

#read in the excel spreadsheet, saved as a .csv, and turn it into a dataframe called survivalrawdata
survivalrawdata=read.csv(file=input.file.name,header=TRUE)
levels(survivalrawdata$Letter)
survivalrawdata$Letter <- factor(survivalrawdata$Letter, levels = c("4", CANDIDATES)) #this is the order groups will appear in graphs
levels(survivalrawdata$Letter)


#convert the data to a list of ages etc. for each worm
#first for the worms that actually died
Ages<-rep(survivalrawdata$Day,survivalrawdata$Died)
Treatments<-as.factor(rep(survivalrawdata$Letter,survivalrawdata$Died))
Censoreds<-rep(1,length(Ages))				                     #1 means not censored (died)
IndividualsDied<-cbind.data.frame(Ages,Treatments,Censoreds)
head(IndividualsDied)
str(IndividualsDied)

#output summary stats
IndividualsDiedDataFrame<-cbind.data.frame(Ages,Treatments,Censoreds)
str(IndividualsDiedDataFrame$Treatments)
#to get the core stats:  n mean sd median trimmed mad(median absolute deviation) min max range skew kurtosis se
SummaryStatsLifespan<-describeBy(IndividualsDiedDataFrame$Ages, IndividualsDiedDataFrame$Treatments, mat = TRUE)
output.file.name<-paste("graphs/",file.name,"_","SummaryStatsLifespan1.csv", sep="")
write.csv(SummaryStatsLifespan, file=output.file.name, row.names = FALSE)

#and again for the censored worms
Ages<-rep(survivalrawdata$Day,survivalrawdata$Censored)
Treatments<-as.factor(rep(survivalrawdata$Letter,survivalrawdata$Censored))
Censoreds<-rep(0,length(Ages))				###0 means censored
IndividualsCensored<-cbind.data.frame(Ages,Treatments,Censoreds)
head(IndividualsCensored)
str(IndividualsCensored$Treatments)

#and now put them together into a data frame
#so that you can refer to the appropriate columns later
LifespanIndividuals<-rbind.data.frame(IndividualsDied,IndividualsCensored)
output.file.name<-paste("graphs/",file.name,"_","Individuals1.csv", sep="")
write.csv(LifespanIndividuals, file=output.file.name, row.names = FALSE)
#this dataframe is one individual per line


#here actually the survival analysis begins
#create a survival object for Kaplan-Meier survival curve
#it is looking at columns of the data frame, NOT the vectors directly.
survivalobject=(Surv(LifespanIndividuals[["Ages"]], 
                     LifespanIndividuals[["Censoreds"]]==1)
                ~(LifespanIndividuals[["Treatments"]]))
#fit the standard Kaplan-Meier survival model to the object
zk.fit=survfit(survivalobject)


#make pairwise plots: wild-type and one mutant at a time using multiplot
output.file.name<-paste("graphs/",file.name,"_","survivalstepNEW.pdf", sep="")
pdf(file=output.file.name, width=WD, height=(HT*length(CANDIDATES)))
listofplots <- list()
for (i in CANDIDATES)
{print(paste("4",i, sep=""))
    part.fit=survfit(survivalobject,
                 subset=((LifespanIndividuals[["Treatments"]]) %in% c("4",i)))
    listofplots[[i]] <- 
    autoplot(part.fit,
             facets = F,
             ncol = 1,
             # surv.linetype = 'dashed',
             surv.geom = 'step',
             # surv.geom = 'point',
             # surv.shape = 1,
             surv.size = 1,
             surv.alpha = 1,
             conf.int = FALSE,
             conf.int.alpha = .25,
             censor = TRUE,
             censor.shape = '+',
             censor.size = 2,
             censor.colour = "gray50",
             xlim = c(0,31),
             # ylim = c(0,1),
             strip_swap = T,
             xlab = "Age (days)",
             ylab = "Survival (%)")+
      scale_colour_hue(h.start=180)
}
multiplotpage <- multiplot(plotlist=listofplots, cols=1)
dev.off()

#make a all plots: wild-type and all mutants
#make pairwise plots: wild-type and one mutant at a time using multiplot
allplot <- autoplot(zk.fit,
             facets = F,
             ncol = 1,
             # surv.linetype = 'dashed',
             surv.geom = 'step',
             # surv.geom = 'point',
             # surv.shape = 1,
             surv.size = 1,
             surv.alpha = 1,
             conf.int = FALSE,
             conf.int.alpha = .25,
             censor = TRUE,
             censor.shape = '+',
             censor.size = 2,
             censor.colour = "gray50",
             #next two lines for setting the axis-limits for comparison
             # xlim = c(0,25),
             # ylim = c(0,1),
             strip_swap = T,
             xlab = "adult age [days]",
             ylab = "survival [%]")+
    scale_colour_grey()
ggsave(filename="graphs/allplot.pdf",plot=allplot,width=WD,height=HT)

#ggsurvplot(survfit(survivalobject=(Surv(LifespanIndividuals[["Ages"]], 
#                                        LifespanIndividuals[["Censoreds"]]==1)
#                                  ~(LifespanIndividuals[["Treatments"]]))))


# END: time taken ---------------------------------------------------------
end.time <- Sys.time()  
time.taken <- end.time - start.time
time.taken

