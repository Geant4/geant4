

//



//
// $Id: gennodelist.h,v 1.3 1999-12-15 18:04:17 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef gennodelist_h
#define gennodelist_h

/*
* NIST Utils Class Library
* clutils/gennodelist.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include "G4ios.hh"

//////////////////////////////////////////////////////////////////////////////
// class GenNodeList
// this class implements a doubly linked List by default.
// If you delete this object it does not delete all of its entries,
// only its head.  If you want it to delete all of its entries as well
// as its head, you need to call DeleteEntries().
//////////////////////////////////////////////////////////////////////////////

class GenNodeList {
public:
    GenNodeList(GenericNode *headNode);
    virtual ~GenNodeList() { delete head; }

    GenericNode *GetHead()	{ return head; }

    virtual void ClearEntries(){;}
    virtual void DeleteEntries(){;}
		// deletes node from its previous List & appends
    virtual void Append(GenericNode *node);
		// deletes newNode from its previous List & inserts in
		//	relation to existNode
    virtual void InsertAfter(GenericNode *newNode, GenericNode *existNode);
    virtual void InsertBefore(GenericNode *newNode, GenericNode *existNode);

    virtual void Remove(GenericNode *node){ if (node) node->Remove(); }

protected:
    GenericNode *head;
};

//////////////////////////////////////////////////////////////////////////////
// class GenNodeList inline functions
// these functions don't rely on any inline functions (its own or
//	other classes) that aren't in this file
//////////////////////////////////////////////////////////////////////////////

inline GenNodeList::GenNodeList(GenericNode *headNode)
{
    head = headNode;
    head->next = head;
    head->prev = head;
}

		// inserts after existNode
inline void GenNodeList::InsertAfter(GenericNode *newNode, 
				     GenericNode *existNode)
{
    newNode->next = existNode->next;
    newNode->next->prev = newNode;

    newNode->prev = existNode;
    existNode->next = newNode;
}

		// inserts before existNode
inline void GenNodeList::InsertBefore(GenericNode *newNode, 
				      GenericNode *existNode)
{
    existNode->prev->next = newNode;
    newNode->prev = existNode->prev;

    newNode->next = existNode;
    existNode->prev = newNode;
}

		// inserts before the head node
inline void GenNodeList::Append(GenericNode *node)
{
    InsertBefore(node, head);
}

#endif
