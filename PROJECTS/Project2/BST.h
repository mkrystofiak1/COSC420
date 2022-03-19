//Author: Mitchell Krystofiak
//Class: COSC420 - Project2
//Date: Novemeber 21, 2021
//Description: Binary Search Tree for keywords.

#ifndef _BST_H
#define _BST_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct keynode
{
	char word[100];
	

	struct keynode *left;
	struct keynode *right;
	
	//Need to add list of associated id's
	

} keynode;



#endif
