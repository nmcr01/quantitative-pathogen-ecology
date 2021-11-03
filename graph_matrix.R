# Script for lecture on graphs and matrices


# Draw the graph of the toy ecosystem
# library(DiagrammeR) un-comment and run if DiagrammeR not
# already loaded to the library


grViz("
   digraph toynet {
   graph [layout=circo]
   node [shape=circle, style=filled, color=grey]
   node [fillcolor=red] # use separate node statements to allocate different
                        # attributes to sets of nodes
   FP
   
   node [fillcolor=blue]
   S1;S2;S3
   
   #edges
   edge [color=grey]   
   FP->{S1 S2 S3}[dir=both][arrowhead=dot][arrowtail=dot]
   S1->{S2,S3}[dir=both][arrowhead=dot][arrowtail=dot]
   S2->S3[dir=both][arrowhead=dot][arrowtail=dot]
   }      
")

# Analysis of the dynamics inherent in the matrix, following May's work
tmax<- 20 #it's a good idea to make a constant that defines the length of time series
times<- seq(2,tmax,1) #define a sequence of time steps
species <-c("FP","S1","S2","S3") #make a character vector of species labels
projection<- matrix(nrow=length(species), ncol=tmax, dimnames=list(species,seq(1,20,1))) # define a 
                            #matrix to hold the output. Each column of the matrix will be one time step
projection[,]<-0 #put zero in every cell of the matrix

projection[,1]<-c(10,1,1,1)
A<-matrix(nrow=length(species), ncol=length(species), dimnames=list(species,species),
                      byrow=TRUE, data=c(-0.1,-0.05,-0.05,-0.05,
                                        -0.2,-0.6,-0.1,-0.1,
                                        -0.15,-0.15,-0.5,-0.15,
                                        -0.01,-0.2,-0.2,-0.3))

for (i in times) {  #Here we see the  use of a for() loop to iterate a calculation
   projection[,i]<-A%*%projection[,i-1] #note that we filled the first column of 
                                        #data outside of the loop, and then started
                                        #counting the loop iterations at i=2, so we
                                        #had already established the first i, i-1 pair.
}

#The next block of commands builds a line plot of each row of the projections
#The example illustrates various line options and one way to add a legend
projplot<-plot(seq(1,20,1), projection[1,], ty="l", lwd=2, col=1,
               ylab="perturbation", xlab="time")
  lines(seq(1,20,1), projection[2,], lwd=2, lty=2, col=2)
  lines(seq(1,20,1), projection[3,], lwd=2, lty=3, col=3)
  lines(seq(1,20,1), projection[4,], lwd=2, lty=4, col=4)
  legend("topright", inset=0.015, legend=species, lty=c(1,2,3,4),
         col=c(1,2,3,4), lwd=2)
  
#Now we want to look at the eigenvalues of A
eigvals_A<-eigen(A)
eigvals_A
max_lambda<-max(eigvals_A$values)
max_lambda