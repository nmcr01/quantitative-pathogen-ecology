library(DiagrammeR)
library(popbio)

# SIR

grViz("
 digraph SIR {
 graph [rankdir = LR]
 node [shape=circle]
 S;I;R
 
 #edges
 S->I [label=pIS][fontsize=10]
 I->R [label=pRI][fontsize=10]
 S->S [label='1-pIS'][fontsize=10]
 I->I [label='1-pRI'][fontsize=10]
 R->R [label=1][fontsize=10]
 }      
      
")
# Simple matrix projection model for SIR process
# Assume 3 states as in standard SIR model: Susceptible, Infected, Removed
# There is no immigration or emigration and no mortality.  We assume that 
# it is possible to define a time interval such that meaningful 
# probabilities can be assumed for the state transitions.

# Start by listing the probabiliites

p_is <-0.25
p_ss <-1-p_is
p_ri <-0.3
p_ii <-1-p_ri
p_rr <-1

# Now define a matrix to take the transition probabiliites

A <-matrix(nrow=3, ncol=3, 
   data=c(p_ss,p_is,0,
          0,p_ii,p_ri,
           0,0,p_rr))
#print A to make sure it has been constucted correctly
A

# Now we need a matrix to hold the output.  Declaring it
# in advance is more efficient than trying to dynamically
# add columns during calculation.

projection <-matrix(nrow=3, ncol=21)

# Now we define the initial states and assign them to the first colomun of
# the output matrix

initial <-c(100,0,0)
projection[,1]<-initial

# The numerical projection of the SIR model is now performed using a 
# "for" loop to iterate a matrix multiplication

for (i in seq(2,21,1)) {
  projection[,i]<-A%*%projection[,i-1]
}

# Make a vector to hold the timesteps and generate a plot of
# the numerical output
t<-seq(1,21,1)
graph1<-plot(t,projection[1,], ty="l", lwd=3, col="darkblue",
             xlab="time", ylab="individuals(%)")
lines(t,projection[2,], lwd=3, lty=2, col="red")
lines(t,projection[3,], lwd=3, lty=3, col="black")

legend("right", inset=0.015, legend=c("S","I","R"), lty=c(1,2,3),
       col=c("darkblue","red","black"), lwd=2)
pop_stats_SIR<-eigen.analysis(A)
pop_stats_SIR
barplot(pop_stats_SIR$stable.stage)


# We can find the time step with the number of infectious
# reaches its peak by using the which.max() function. This
# finds the index value of a vector at which the highest
# value of the vector occurs.

peak_I<-which.max(projection[2,])
peak_I
abline(v=peak_I, col="darkgray", lty=3, lwd=3, ylim=c(0,100))

# A little detail to end, we output the spectrum of A and
# see that the dominant eigenvalue =1, which is what we'd
# expect since the matrix is only shuffling individuals 
# among states  - there's no growth as such. The barplot
# of the stable stage distribution confirms what the 
# time plot already showed - that R is an absorbing
# state and everyone ends up in R eventually.

pop_stats_SIR<-eigen.analysis(A)
pop_stats_SIR
barplot(pop_stats_SIR$stable.stage)


# Making an SIRS model
# If we want to have some of the recovered individuals 
# become susceptible again, we need to add a transition
# probability to the matrix to allow for that eventuality.
# Let's preserve the initial example, so we make a copy
# of A and assign it to a new object before editing it.

A2 <-A

# Now we want to give the transition from R to S a 
# probability>0 and put this in A2.  We have to adjust
# the entry for remaining in R because R is no longer
# an absorbing state.
grViz("
 digraph SIRS {
 graph [rankdir = LR]
 node [shape=circle]
 S;I;R
 
 #edges
 S->I [label=pIS][fontsize=10]
 I->R [label=pRI][fontsize=10]
 S->S [label='1-pIS'][fontsize=10]
 I->I [label='1-pRI'][fontsize=10]
 R->R [label='1-pSR'][fontsize=10]
 R->S [label=pSR][fontsize=10]
 }      
      
")
p_sr <-0.05
A2[1,3]<-p_sr
A2[3,3]<-1-p_sr
A2

# Now we repeat the projection with the new matrix
projection2 <-matrix(nrow=3, ncol=21)
projection2[,1]<-projection[,1]

for (i in seq(2,21,1)) {
  projection2[,i]<-A2%*%projection2[,i-1]
}


pop_stats_SIRS<-eigen.analysis(A2)
pop_stats_SIRS
barplot(pop_stats_SIRS$stable.stage,
        col=c("darkblue","red","black"))


graph2<-plot(t,projection2[1,], ty="l", lwd=3, col="darkblue",
             xlab="time", ylab="individuals (%)", ylim=c(0,100))
lines(t,projection2[2,], lwd=3, lty=2, col="red")
lines(t,projection2[3,], lwd=3, lty=3, col="black")
legend("right", inset=0.015, legend=c("S","I","R"), lty=c(1,2,3),
       col=c("darkblue","red","black"), lwd=2)
peak_I2<-which.max(projection2[2,])
peak_I2
abline(v=peak_I2, col="darkgray", lty=3, lwd=3, ylim=c(0,100))