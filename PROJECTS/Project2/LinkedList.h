//Author: Mitchell Krystofiak
//Class: COSC420 - Project2
//Date: November 21, 2021
//Description: Linked List for connecting associated id's.

#include <stdio.h>
#include <stdlib.h>

typedef struct lnode
{
	//calculated from write.c + 1
	char id[17];
	char title[383];
	char author[36663];
	char abstract[6090];

	lnode *next;
	lnode *prev;
	lnode *head;

} lnode;

