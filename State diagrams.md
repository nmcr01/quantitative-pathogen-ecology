---
title: "State Diagrams"
author: Neil McRoberts, Robin Choudhury
output: html_document

---

### A diagram is a good first model

We're going to be making and looking at a lot of mathematical and computer models in this class, but for any biological problem you're likely to work on it's usually a good idea to start with a visual model first and work on translating that into a computable[^1] model later.  Sometimes, just drawing the visual model alone is enough for us to see that we're making conceptual mistakes, or to highlight where we are seriously lacking information.

It's surprising how much biology you can represent with simple diagrams that consist of nothing more than *nodes*[^2] and *edges*.  This type of diagram turns up in lots of contexts.  An obvious example would be any situation where a network is present, which means we can be talking about anything from molecules to species.  One consequence of this is that it's important to make it clear what the nodes are and what the edges mean; edges can simply capture the idea that two nodes are related in someway, or they can convey the idea of causal relationships or flow of something from one node to another.

It should be obvious that if you're going to spend a lot of time drawing diagrams it's a good idea to have software that does a good job of that.  We'll be taking a look at some of the capabilities of R for network analysis later in the class schedule.  Some of the packages for network analysis can handle data sets with thousands of nodes and millions of edges, but for now we'll confine ourselves to thinking about relatively simple diagrams that you could conceivably draw with pencil on paper.

Outside of the R ecosystem there are commonly used office tools such as Powerpoint, Google Slides, LibreOffice Presents and other slide making tools that have some diagramming tools built in.  If nothing else is available, these are a good fallback.  In the world of network analysis, the online tool [Kumu](kumu.io) is popular, powerful and capable of producing really eye-catching diagrams.  The [MentalModeler](https://www.mentalmodeler.com/#home) tool works in web browsers and is specifically aimed at diagrams where the aim is to capture cause and effect relationships - so called cognitive maps.  From experience, of the general diagramming tools available, we think [Lucidchart](https://www.lucidchart.com/pages/) works well and produces clear diagrams in a variety of formats. 

For R there is a package called DiagrammeR that lets you build complex diagrams interactively.  The class is going to use DiagrammeR to help us think about state diagrams as the first step in building models in biology.  As a heads up, where we're heading after that is to use the homology between graphs (in the mathematical sense) and matrices to translate diagrams into computable models.  That will set us up to take a look at one of the classic texts in theoretical ecology as a starting point for thinking about species and community dynamics.

### DiagrammeR
R is a modular system.  Specialist pieces of analysis are usually found in *packages* which can be installed from one of the many online repositories. The base R installation comes with a wide range of common analyses and tools built in, but for something like specialist diagram design, or advanced bioinoformatics, generalized linear modeling, or numerical integration, *etc*, *etc*, you will need to download and install relevant packages.

In R Studio the process for downloading a new package and adding it to your library is partly automated. To install a new package your computer needs to be connected to the internet.  If you click on the **Tools** menu item you should see that the first option on the menu is **install packages**. Selecting that option will open a dialogue box that will let you type in the name of the package you want to install. The box uses predictive text so it will start to offer you choices based on spelling.  Once the name of the package you want appears in the list, select it.  In this case you will (obviously) want to type "diag..." and then select **DiagrammeR** from the list[^3] and click on **install**.

Installing a package is the first step in getting access to the functions and methods it contains.  To actually use the package in a session or running a program, you also have to load it into the *library* of active packages.  There are two ways to do that in R Studio.  The first approach is to add a **library** command to your code before the commands which use the package:

    library(DiagrammeR)
This method works for any package you want to add to the active library.  You can add library() commands at any point in  program, but it's often a good idea to put all "*housekeeping*" code close to the top of the program and then focus on the instructions/calculations that will do the actual analysis. The other way to load a package is to use the **Packages** window.  In a standard R Studio setup the **Packages** window is one of 5 different tabs that share the lower right pane of the application's screen.  If you click on the Packages tab you'll see an alphabetical list of all the packages that are installed on the computer. The ones that are loaded in the library will have a tick mark in box next to their name.  Scroll down the list until you find the package you want and click on the box to activate it.  Note that when you do this, you'll see the corresponding command appear in the **Console** window.  Again, this approach will work for any package you want to load during an interactive session, but it doesn't add the command to your program code.

### Examples
The examples are provided in the file DiagrammeR.R which you can get from the Github repository for the class or from Canvas.  The example file only provides a simple introduction to the range of diagram drawing options available in DiagrammeR.  The first thing to note is that DiagrammeR uses a third-party graphics tool called [GraphViz](https://graphviz.org/) which is available as a separate graphics toolbox for a wide range of programming languages (_e.g._ there are implementations for Python and C).

In R (assuming that DiagrammeR is loaded in the active R session), you start drawing a graph by using the command

    grViz(" ... ")
  where the commands to produce the graph are written where the three dots are.  To make R process the code in the program file you can go to each line one at a time and hit <Ctrl><Enter>, or you can highlight a whole block of code and hit <Ctrl><Enter>, or use the commands from the **Code** menu.  To see the graph you will need to click on the **Viewer** window in the lower right pane in R Studio.  The code for one version of a simple SIR model looks like this:
  
     grViz("
      digraph SIR {
      graph [rankdir = LR]
     
      #nodes
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

Here is another example, motivated by ecological communities.  In this case we are assuming a little ecosystem of 4 species; one primary resource (R1), two primary consumers (C1 and C2) and one predator (P1).

       grViz("
       digraph econet {
       graph [rankdir=TB]
       node [shape=square]
       R1;C1;C2;P1
       
       #edges
     R1->C1 [label=c1][arrowhead=none][fontsize=10]
     R1->C2 [label=c2][arrowhead=none][fontsize=10]
     C1->P1 [label=m1][arrowhead=none][fontsize=10]
     C2->P1 [label=m2][arrowhead=none][fontsize=10]
     C1->C2 [label=b][arrowhead=none][fontsize=10]
       }      
    ")

Notice in the second example the edge statements contain another option [arrowhead=...] with the value set to "none" which results in the edges being simple lines.  This raises the point that the edges (and, indeed the nodes) in a graph can convey meaning depending on the visual syntax used.  In some disciplines the visual syntax is standardized into generally accepted sets of symbols, in others there is standard and each scientist makes their own decisions about how (or even, if) to add meaning to graphs by using options added to nodes and edges.  In the case of DiagrammeR [this page](https://rich-iannone.github.io/DiagrammeR/graphviz_and_mermaid.html) gives a good starting point for seeing the options that are available.

### Exercise
For the little ecological network example add shapes to the ends of the edges to indicate the type of interaction that would occur in a real network of this type.  As an alternative you can start from scratch and make a model of a network that is relevant to your own research.


[^1]: To save typing "mathematical and computer" or "statistical"  or other longhand ways of referring to models which can be calculated or solved in some way, we're going to use the word *computable* as a collective adjective for all of those.  It's a bit of loose terminology but as long as we're clear what we mean it'll save a lot of space.
[^2]: In [graph theory](https://en.wikipedia.org/wiki/Graph_theory), nodes are called vertices (singular, vertex).
[^3]: Note that the dialogue box is not case-sensitive, but R program code **is** case sensitive

> Written with [StackEdit](https://stackedit.io/).
