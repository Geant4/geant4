

//



//
// $Id: dispnodelist.cc,v 1.2 1999-05-21 20:21:08 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/dispnodelist.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

#include <gennode.h>
#include <gennodelist.h>

#include <mgrnode.h>
#include <mgrnodelist.h>
#include <dispnode.h>
#include <dispnodelist.h>

	// deletes node from its previous List & appends
	// actually it puts it at the front of the List.
void DisplayNodeList::Append(GenericNode *node)
{
    InsertBefore(node, head);
}

		// deletes newNode from its previous List & inserts after
		//	existNode
void DisplayNodeList::InsertAfter(GenericNode *newNode, 
					 GenericNode *existNode)
{
    if(newNode->next != 0){	// remove the node from its previous 
	newNode->Remove();	//	display state List
    }
    GenNodeList::InsertAfter(newNode, existNode);
// DON'T DO THIS    ((DisplayNode *)newNode)->displayState = listType;
}

		// deletes newNode from its previous List & inserts before
		//	existNode
void DisplayNodeList::InsertBefore(GenericNode *newNode,
					  GenericNode *existNode)
{
    if(newNode->next != 0){	// remove the node from its previous 
	newNode->Remove();	//	display state List
    }
    GenNodeList::InsertBefore(newNode, existNode);
// DON'T DO THIS!!!   ((DisplayNode *)newNode)->displayState = listType;
}

