

//



//
// $Id: mgrnodelist.cc,v 1.2 1999-05-21 20:21:10 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/mgrnodelist.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

#include <mgrnode.h>
#include <mgrnodelist.h>
#include <dispnode.h>
#include <dispnodelist.h>


MgrNodeList::MgrNodeList(stateEnum type) : GenNodeList(new MgrNode())
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNodeList::MgrNodeList()\n";
    listType = type;
    ((MgrNode *)head)->currState = type;
}

void MgrNodeList::Remove(GenericNode *node)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNodeList::Remove()\n";
    GenNodeList::Remove(node);
// DON'T DO THIS    ((MgrNode *)node)->currState = noStateSE;
}

	// deletes node from its previous List & appends...
	// actually it puts it at the front of the List.
void MgrNodeList::Append(GenericNode *node)
{
    InsertBefore(node, head);
}

		// deletes newNode from its previous List & inserts after
		//	existNode
void MgrNodeList::InsertAfter(GenericNode *newNode, 
				     GenericNode *existNode)
{
    if(newNode->next != 0){	// remove the node from its previous List
	newNode->Remove();
    }
    GenNodeList::InsertAfter(newNode, existNode);
// DON'T DO THIS    ((MgrNode *)newNode)->currState = listType;
}

		// deletes newNode from its previous List & inserts before
		//	existNode
void MgrNodeList::InsertBefore(GenericNode *newNode, 
				      GenericNode *existNode)
{
    if(newNode->next != 0){	// remove the node from its previous 
	newNode->Remove();	//	state List
    }
    GenNodeList::InsertBefore(newNode, existNode);
// DON'T DO THIS!!    ((MgrNode *)newNode)->currState = listType;
}

MgrNode *MgrNodeList::FindFileId(int fileId)
{
    MgrNode *mn = (MgrNode *)head->next;
    while(mn != head)
    {
	if(mn->GetFileId() == fileId)
	    return mn;
	mn = (MgrNode *)mn->next;
    }
    return (MgrNode *)0;
}

