# Graph layering (Sugiyama Method) 

Step 1. Longest path layering (layer assignment)
Step 2. Barycenter method (edge crossing minimization)
Step 3. Priority layout method (horizontal coordinates assignment)

"With dummy vertices" shows virtual vertices created by the algorithm. It is a dangerous checkbox! It will remove some of the original edges. And on the next auto layout dummy vertices will be dropped.

Many thanks to Ross Kirsling for making a nice [graph editor example]((http://bl.ocks.org/rkirsling/)), which is both simple and functional (and ES6!).
