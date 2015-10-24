
==================shell======
#!/bin/bash

EdgeList= $1
java -jar encode.jar $EdgeList
mace Me_ $EdgeList.en $EdgeList.en.cliques
java -jar decode.jar $EdgeList.map $EdgeList.en.cliques
==============================


For an edge file, we do 3-steps to find the maximal cliques:
eg.
$ cat edge
0	1
0	2
1	2
2	3
3	4
1	3
7	8

1. Encode the nodes in edge file:
===============================================================

$ java -jar bin/encode.jar edge
Map Size: 7
Map file: edge.map
Encoded file: edge.en
NodeVetexCount: edge.nv

$ cat edge.en
0	1
0	2
1	2
2	3
3	4
1	3
5	6


2. Find Connected Component (subgraph) 
===============================================================

$ conngraph edge.nv edge.en > edge.en.cg
$ cat edge.en.cg
0	1	2	3	4
5	6


3. Use mace to find the maximal cliques in the encoded file
===============================================================
genarate Maximal Clique file: edge.en.cliques:

$ .bin/mace Me_ edge.en

$ cat edge.en.cliques
4	3
3 	2	1
2	1	0
6	5

4. Decode the clique file
===============================================================

$ java -jar bin/decode.jar edge.map edge.en.cliques 
Map Size:7
Decoded file: edge.en.cliques.de
Attention: The 1st column is the size of clique

$ cat edge.en.cliques.de
2	4	3
3	3	2	1
3	2	1	0
2	8	7
 
 


