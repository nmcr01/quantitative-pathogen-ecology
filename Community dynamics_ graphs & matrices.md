---
title: "Community dynamics: graphs & matrices"
author: Neil McRoberts, Robin Choudhury
output: html_document

---
### Introduction
In this lecture we will be taking a brief look at some ideas introduced by physicist turned theoretical ecologist the late [Sir Robert May](https://en.wikipedia.org/wiki/Robert_May,_Baron_May_of_Oxford).  The material that we'll cover can be found in the first couple of chapters in May's classic monograph *Stability & Complexity in Model Ecosystems* (SCIME) which was first published in 1973.  The book was a synthesis of ideas that May had developed during the preceding three or four years and had published in a series of papers that appeared between 1970 and 1972. For the purposes of this short course there are two main motivations to tackle this subject. First, the ideas May discusses provide a nice example of the homology between graphs and matrices, so following the process that May describes for thinking about disturbance in ecological communities, is a convenient way to learn about that useful path for turning concepts into models; we will revisit it when we think about cognitive maps towards the end of the course.  Second, May's ideas have direct application in community ecology and are useful for anyone thinking about topics such as phytobiomes, biological control, or microbial ecology.

At a more general level, May's work exemplifies a style of science that it is useful to understand and in which to have some ability; _i.e._ the ability to formulate a theoretical conjecture about how some aspect of the world might work and follow through the reasoning to the point of being able to identify features that you would expect to be able to observe if your conjecture were true. It's also important to bear the historical context of the publication of SCIME in mind.  At the time May wrote the book he was in transition from his early career in physics to becoming a theoretical ecologist.  Although theoretical ecology existed as a subject already, May's work moved the field forward by a huge step, as a result of his style of working quantitatively but in a very abstract way.


### Graphs, matrices and communities
Let's start by noting the problem that May was trying to unravel in his work.  Is there any evidence to support the idea that more diverse communities are, as a direct consequence of being more diverse, more stable? To retrace May's analysis we are going to start with a toy community comprising just four species. To link our work to the previous module on microbial ecology, we'll assume that the four species are members of the phylosphere surface microbial ecosystem on the leaves of a species of crop plant.  The analyses themselves are completely unaffected by these particulars, and could apply equally as well to the macro-fauna on a grassland, or the gut inhabitants of an ocean crustacean.

Let's assume that our four species comprise three specialist saprophytes and one facultative plant pathogen that can exist as a surface inhabitant, but which can also adopt a pathogenic mode under conditions that promote a switch in behavior.  That aspect won't be important for our initial analysis.  We assume that all four of the species competes with the other three and so the presence of each species serves to reduce the potential growth of the others.  The first job is draw a graph of this ecological community to represent these interactions.  The set of possible variations on how to draw such a graph is huge, but the useful ones will share a common set of attributes the matter if the graph is to perform as a valid model for the toy ecosystem:
 - each species should be represented by one, and only one, node
 - the nodes for each species should be connected to all three of the other nodes
 - Either by using pairs of directions edges, or by applying appropriate attributes at the ends of edges, the edges should convey the information that there is mutual competition between all pairs of species.

The following code snippet will produce one appropriate graph for this community.
   
+------------------+

grViz("
   digraph toynet {
   graph [layout=circo]
   node [shape=circle, style=filled, color=grey]
   node [fillcolor=red]
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

±-----------------+    
We can "read" the graph as telling us that there are four species which interact and they are all mutually competitive.  Let's suppose that it's possible for such a community to obtain an equilibrium state in which all four species are present.  We can think about community stability by asking what would happen if there was a disturbance, or perturbation, to the equilibrium that caused the population sizes of the four species to be changed from their equilibrium values.  Would the disturbance be temporary, with the population sizes returning to their equilibrium values?  Would the disturbance grow until one or more of the species is eliminated from the community?  Or, perhaps, the disturbance would neither be damped out, or grow out of control, but the populations might settle into a new equilibrium state, with different population sizes than before, or enter into a pattern of endless oscillations in which no species is ever eliminated but the community as a whole never settles to a stable distribution of species densities?  May's analysis was designed to provide a way to decide among these possibilities by looking at what would happen to the perturbation in initial region close to the equilibrium, when the absolute size of the disturbance is still small, and by assuming that in a suitably small interval we can assume that changes in the perturbation are linear.

From here on, we are going to use the same notation that May used so that you can compare these notes directly with his work in SCIME.  We gave our species labels (*e.g.* FP).  May refers to species in a community using italic lowercase letters $i..., j$.  The interactions between species are indicated by a set of coefficients that measure the effect of species $i$ on species $j$ when they are interacting.  May uses the Greek letter alpha,$\alpha$ (lowercase) to refer generically to these interactions and adds subscripts to indicate which pair of species is involved in the interaction: $\alpha_{ij}$, the effect of species $j$ on species $i$.  The key idea is that these are the quantities represented by the ***edges*** in our graph.

With the knowledge that the edges in the graph can be thought of as representing interaction coefficients it's a relatively easy step to see that we can capture all of the information in the graph by writing a square matrix with four rows and four columns in which each row (or column) refers to one species.  If the cells of the matrix are used to store the relevant $\alpha_{ij}$ values, the graph and the matrix will contain exactly the same information; *i.e.* they are homologs.  By convention we construct the matrix so that the columns represent the effect of each species on each of the others.  The resulting matrix is referred to as the *community matrix*.  It is a mathematical object that is a simple model of a set of ecological interactions.  May uses the uppercase Roman letter $A$ to denote the generic community matrix containing a set of $\alpha_{ij}$ coefficients.

If we substitute our labels for the species for the generic $ij$ notation used by May, the community matrix for our toy ecosystem can be written as[^1]:$$
A=\begin{pmatrix}
\alpha_{FPFP} & \alpha_{FPS1} & \alpha_{FPS2} &\alpha_{FPS3} \\
\alpha_{S1FP} & \alpha_{S1S1} & \alpha_{S1S2} &\alpha_{S1S3} \\
\alpha_{S2FP} & \alpha_{S2S1} & \alpha_{S2S2} &\alpha_{S2S3} \\
\alpha_{S3FP} & \alpha_{S3S1} & \alpha_{S3S2} &\alpha_{S3S3} \\
\end{pmatrix} $$

Notice that in completing $A$ we have filled in entries on the main diagonal (i.e. the diagonal running from top left to bottom right).  These are the effects of each species on it's own population size near equilibrium.  In order for the community to be stable these have to be negative values so that the species impose some kind of self-regulation[^2].  With $A$ organized in this way, each of its rows can be thought of as containing the coefficients of a linear difference equation that defines what the next value of the population perturbation will be given the current value.  This allows us to use $A$ as a linear operator (or multiplication device) to propagate a vector containing the initial values of the perturbation over time.  Following May's notation we will use the symbol $\mathbf{x}(t)$ to denote a column vector $\mathbf{x}$ which contains the values for our four species, one value per row in the same order as the rows of $A$; *i.e.* FP, S1,S2, S3.  Successive values for $\mathbf{x}$ are given by the matrix-vector product$$
\mathbf{x}(t+1)=A \cdot \mathbf{x}(t)$$

This gives us a straightforward numerical method for looking at the fate of disturbances in population sizes for ecological communities as long as we have some way of estimating the values for the interaction coefficients, $\alpha_{ij}$.  However, if we know what those value are, or if we can make good guesses, we don't actually need to do the numerical calculations to answer basic questions about the fate of the community if it is disturbed.

#### Analysis of the community matrix
*If* $A$ is a useful model for community dynamics close to equilibrium, May noted that the fate of the perturbation depends on a set of quantities known as the *eigenvalues* of the matrix[^3] $A$, which are traditionally denoted using the lowercase Greek letter lambda, $\lambda$.   For this type of application there will be as many eigenvalues of $A$ as there are rows.  The key result of interest is that the value of the largest eigenvalue determines the growth of the perturbation with time.  If $\lambda_{max}>0$, the perturbation will grow with time and the community is not stable when disturbed from equilibrium[^4].  If $\lambda_{max}<0$, the perturbation will die out with time and equilibrium will return.  In the rare case where $\lambda_{max}=0$ the system displays what is called "neutral stability" and once dislodged from equilibrium it will oscillate endlessly without any of the components being eliminated.

R, like most software for numerical and mathematical analysis, has a set of standard functions for handling matrix eigenvalues so following May's approach to the analysis of community stability is relatively simple.  The following code snippet provides a simple example based on the  toy community we have been considering up till now.  The code is also in the file *graph_matrix.R*.
 

     # Analysis of the dynamics inherent in the matrix, following May's work
    tmax<- 20     #it's a good idea to make a constant that defines the length of time series
    times<- seq(2,tmax,1)    #define a sequence of time steps
    species <-c("FP","S1","S2","S3")    #make a character vector of species labels
    projection<- matrix(nrow=length(species), ncol=tmax, dimnames=list(species,seq(1,20,1))) # define a 
                                #matrix to hold the output. Each column of the matrix will be one time step
    projection[,]<-0            #put zero in every cell of the matrix
    
    projection[,1]<-c(10,1,1,1)
    A<-matrix(nrow=length(species), ncol=length(species), dimnames=list(species,species),
                          byrow=TRUE, data=c(0.1,0.05,0.05,0.05,
                                            0.2,0.6,0.1,0.1,
                                            0.15,0.15,0.6,0.15,
                                            0.25,0.2,0.2,0.3))
    
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

Remember that May established that in order for the system to be stable, the largest eigenvalue has to have a negative (real) value.  The toy model displays this characteristic and the time plot shows the initial perturbation dying away in oscillations of decreasing amplitude.

### Exercises
With some code to work from, there are lots of things you can do to explore this little toy model.  What happens if you change the interaction coefficients?  What would happen if one of the four species wasn't present?  Would the system with the remaining 3 species still be stable?

[^1]: If you're reading this in the pdf version, refer to the Markdown file to see the TeX syntax for laying out the matrix
[^2]:  One widely held tenet of population ecology is that competition within species must be stronger than competition between species.  This line of thought comes from the idea that the strongest competition will come from individuals who exactly share the same ecological niche, and that the individuals who best fit that definition will be other members of the same species.  However, that line of thought seems to contradict the extremely common observation that individuals of the same species often occur in patches, herds, shoals, flocks, families and other dense groupings.
[^3]: We are not going to delve into the properties or derivation of matrix eigenvalues in this class.  They are important in many applications of linear algebra.  For our purposes we will simply note that the eigenvalues of $A$ are related to the interaction coefficients represented by the edges of the graph through the equation $$ 
\lambda x_{i]}(t) = \sum_{j=1} ^{m} \alpha_{ij} x_{j}(t)
$$
[^4]: This is not strictly true.  Eigenvalues can be *[complex* numbers](https://en.wikipedia.org/wiki/Complex_number) with both real and imaginary parts.  The growth of the perturbation depends on only whether the real part is greater than or less than 0; the imaginary parts determine the period of any oscillations in the system.  The special case of neutral stability occurs when one eigenvalue has a real part equal to 0 and all the other eigenvalues have only imaginary parts.

> Written with [StackEdit](https://stackedit.io/).
