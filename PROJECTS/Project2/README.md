Authors: Mitchell Krystfiak, Olivia Rorke, Lens Desrouleaux
Class: COSC420 - Project2 Progress Report
Date: November 21, 2021


Files Included: matrix.c, matrix.h, BST.h, LinkedList.h, write.c

 matrix.c/.h:

 	- Matrx operations from previous labs
        - Started adjacency matrix function, which reads through the 
          arxivcitations.txt file and builds and adjacency matrix. 
	  This step is needed for many things, including implementing
	  the HITS algorithm, implementing the PageRank algorithm,
	  and implementing the sparse matrix functionality.

 BST.h, LinkedList.h:

 	- Unfinished files to build the indexing functionality.
	- This is going to be used to keep track of keywords in each
	  paper. Each node in the Tree is going to keep a keyword and
	  its associated linked list of the data nodes that share this
	  keyword.
	
 write.c: 

 	- This program calculates the largest lines in the arxiv-metadata.txt
	  file so that we can allocate a spot for these lines. We will probably
	  make a node structure for each paper that holds the id, the title,
	  the authors, and the abstract for quick referencing.
	- This may also be used to split up the large metadata file for more efficient
	  usage and/or test cases.



           
