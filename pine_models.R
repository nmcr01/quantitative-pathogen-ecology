library(DiagrammeR)

grViz("
   digraph pitch_canker {
   graph [rankdir=TD]
   node [shape=circle, style=filled, color=grey]
   node [fillcolor=Sienna]
   Cone
   node [fillcolor=Goldenrod]
   Seed
   node [fillcolor=YellowGreen]
   Seedling
   node [fillcolor=Green]
   SSap; MSap; LSap
   node [fillcolor=ForestGreen]
   Mature
   
   #edges
   edge [color=grey]   
   Cone->Seed [arrowhead=vee][label=R]
   Seed->Seedling[arrowhead=vee][label='p(seed-seedling)']
   Seedling->SSap[arrowhead=vee][label='p(seedling-ssapling)']
   SSap->MSap[arrowhead=vee][label='p(ssapling-msapling)']
   MSap->LSap[arrowhead=vee][label='p(msapling-lsapling)']
   LSap->Mature[arrowhead=vee][label='p(LSap-Mature)']
   SSap->Cone[arrowhead=box][label=fSsap]
   MSap->Cone[arrowhead=box][label=fMsap]
   LSap->Cone[arrowhead=box][label=fLsap]
  Mature->Cone[arrowhead=box][label=fMature]
   
   }      
")


# Analysis of the dynamics inherent in the matrix.
tmax<- 20     #At one decade per step, this is equivalent to 200 years
times<- seq(2,tmax,1)    #define a sequence of time steps
stages <-c("Cone","Seed","Seedling","SSap","MSap","LSap","Mature")    # stage labels

# define a matrix to hold the output. Each column of the matrix will be one time step
pine_proj<- matrix(nrow=length(stages), ncol=tmax, dimnames=list(stages,seq(1,tmax,1)))

pine_proj[,]<-0  #put zero in every cell of the matrix

pine_proj[,1]<-c(20,1000,1000,10,2,2,1) #define the initial population

#Make the projection matrix and fill it with interaction coefficients
A<-matrix(nrow=length(stages), ncol=length(stages), dimnames=list(stages,stages),
          byrow=TRUE, data=c(0,0,0,20,25,50,75,
                             40,0,0,0,0,0,0,
                             0,0.5,0,0,0,0,0,
                             0,0,0.001,0,0,0,0,
                             0,0,0,0.05,0.5,0,0,
                             0,0,0,0,0.5,0.6,0,
                             0,0,0,0,0,0.4,0.90))

for (i in times) {  #Here we see the  use of a for() loop to iterate a calculation
  pine_proj[,i]<-A%*%pine_proj[,i-1] #note that we filled the first column of 
  #data outside of the loop, and then started
  #counting the loop iterations at i=2, so we
  #had already established the first i, i-1 pair.
}

#The next block of commands builds a line plot of each row of the projections
#The example illustrates various line options and one way to add a legend

pineplot<-plot(seq(1,20,1), log(pine_proj[1,]), ty="l", lwd=2, col=1,
               ylab="ln(individuals)", xlab="decades", ylim=c(-2,10))
lines(seq(1,20,1), log(pine_proj[2,]), lwd=2, lty=2, col=2)
lines(seq(1,20,1), log(pine_proj[3,]), lwd=2, lty=3, col=3)
lines(seq(1,20,1), log(pine_proj[4,]), lwd=2, lty=4, col=4)
lines(seq(1,20,1), log(pine_proj[5,]), lwd=2, lty=1, col=5)
lines(seq(1,20,1), log(pine_proj[6,]), lwd=2, lty=2, col=6)
lines(seq(1,20,1), log(pine_proj[7,]), lwd=2, lty=3, col=9)

legend("right", inset=0.015, legend=stages, lty=c(1,2,3,4,1,2,3),
       col=c(1,2,3,4,5,6,9), lwd=2)

#Now we want to look at the eigenvalues of A
eigvals_A<-eigen(A)
eigvals_A

# Do the eigen analysis using popbio
# First install popbio from the R repository then load it
# install.packages("popbio")
library(popbio)
pop_stats_Pradiata<-eigen.analysis(A)
barplot(pop_stats_Pradiata$stable.stage)
barplot(log(pop_stats_Pradiata$repro.value))

#Make the new projection matrix and and the effect of disease incidence, p.
pine_projD<- matrix(nrow=length(stages), ncol=tmax, dimnames=list(stages,seq(1,tmax,1)))
pine_projD[,]<-0  #put zero in every cell of the matrix
pine_projD[,1]<-c(20,1000,1000,10,2,2,1) #define the initial population
p<-0.1 #assumed disease incidence across the population
AD<-matrix(nrow=length(stages), ncol=length(stages), dimnames=list(stages,stages),
           byrow=TRUE, data=c(0,0,0,20*(1-p),25*(1-p),50*(1-p),75*(1-p),
                              40,0,0,0,0,0,0,
                              0,0.5*(1-p),0,0,0,0,0,
                              0,0,0.001*(1-p),0,0,0,0,
                              0,0,0,0.05*(1-p),0.5*(1-p),0,0,
                              0,0,0,0,0.5*(1-p),0.6,0,
                              0,0,0,0,0,0.4,0.90))

for (i in times) {  #Here we see the  use of a for() loop to iterate a calculation
  pine_projD[,i]<-AD%*%pine_projD[,i-1] #note that we filled the first column of 
  #data outside of the loop, and then started
  #counting the loop iterations at i=2, so we
  #had already established the first i, i-1 pair.
}

#The next block of commands builds a line plot of each row of the projections
#The example illustrates various line options and one way to add a legend

pineplotD<-plot(seq(1,20,1), log(pine_projD[1,]), ty="l", lwd=2, col=1,
                ylab="ln(individuals)", xlab="decades", ylim=c(-2,10))
lines(seq(1,20,1), log(pine_projD[2,]), lwd=2, lty=2, col=2)
lines(seq(1,20,1), log(pine_projD[3,]), lwd=2, lty=3, col=3)
lines(seq(1,20,1), log(pine_projD[4,]), lwd=2, lty=4, col=4)
lines(seq(1,20,1), log(pine_projD[5,]), lwd=2, lty=1, col=5)
lines(seq(1,20,1), log(pine_projD[6,]), lwd=2, lty=2, col=6)
lines(seq(1,20,1), log(pine_projD[7,]), lwd=2, lty=3, col=9)

legend("right", inset=0.015, legend=stages, lty=c(1,2,3,4,1,2,3),
       col=c(1,2,3,4,5,6,9), lwd=2)
pop_stats_PradiataD<-eigen.analysis(AD)
pop_stats_PradiataD
barplot(pop_stats_PradiataD$stable.stage)
barplot(log(pop_stats_PradiataD$repro.value))
