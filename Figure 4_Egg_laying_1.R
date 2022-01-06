# Intro -------------------------------------------------------------------
# clears workspace
rm(list = ls())
#load IndividualWorms.csv file and cut off the header and the first three summary columns
ind0=read.csv(file="file_location",skip=61,header=T)
ind=ind0[,-c(1,2,3)]
state=c("A","+")
nvar=4

loc1=which(ind[1,]==state[1])+nvar #ind[1,] gets row 1 and all egg-laying columns
for(i in 2:length(state)){
  loc1tmp=which(ind[1,]==state[i])+nvar #in welcher column taucht das Zeichen + fuer tot auf
  loc1=c(loc1,loc1tmp) #alle egg-laying columns und die nummer der columns wo der wurm stirbt
}
loc=c(loc1,rep(0,(dim(ind)[2]-length(loc1))))
#find age of mature then plus 4 to get location of 1st egg laying step
# e.g. plus 3 to get location of food eaten and so on
# fill the empty loci with 0
# a sequence about 1st worm indicate "location" (location but not the #eggs yet) 

for(j in 2:dim(ind)[1]){  #anzahl rows in ind)
  locj=which(ind[j,]==state[1])+nvar
  
  for(i in 2:length(state)){
    locjtmp=which(ind[j,]==state[i])+nvar
    locj=c(locj,locjtmp)
  }
  
  locjj=c(locj,rep(0,(dim(ind)[2]-length(locj))))
  loc=rbind(loc,locjj)
}
#loc: a matrix of "location of #eggs" for all worms

locvar=loc+1
indnew=cbind(rep(0,dim(ind)[1]),ind)
locvarnew=locvar[,-which(colSums(locvar)==dim(locvar)[1])]
# if elements of a column is all 0, then it will be all 1 after locvar=loc+1
# and colsum will equal to its length
# remove the extra columns which nodes are all 0 


tablenew=locvarnew
ndeath=rep(0,dim(locvarnew)[2])
for(n in 1:dim(tablenew)[2]){
  for(m in 1:dim(tablenew)[1]){tablenew[m,n]=indnew[m,locvarnew[m,n]]}
  ndeath[n]=sum(which(locvarnew[n]==1))
}
#calculate #death within each column as previous, the dead node would be 1

sum=colSums(tablenew)
nalive=dim(locvarnew)[1]-ndeath
mean=sum/nalive
mean_short=mean[1:120]
mean_data=data.frame(times=c(1:120), mean_eggs=mean_short)
timesteps_day <- 8;
mean_data2=aggregate(mean_data,list(rep(1:(nrow(mean_data)%/%timesteps_day+1),each=timesteps_day,len=nrow(mean_data))),sum)[-1]
mean_data2=transform(mean_data2, days=c(1:15))


plot(mean_data2$days, mean_data2$mean_eggs, type="l", main="egg-laying", xlab="Time(day)", ylab="Number of Eggs")


#average total number of eggs
number_indiv = nrow(ind)#number of inidivudals in the table
total_number_eggs = sum(sum)/number_indiv

number_indiv
total_number_eggs

