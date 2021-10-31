# Graphs using DiagrammeR
library(DiagrammeR)

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

# A little ecological community with one primary producer, two consumers,
# and a predator

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
