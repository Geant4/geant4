

//



//
// $Id: dispnodelist.h,v 1.2 1999-05-21 20:20:39 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef dispnodelist_h
#define dispnodelist_h

/*
* NIST STEP Editor Class Library
* cleditor/dispnodelist.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <gennode.h>
#include <editordefines.h>
#include <mgrnode.h>
#include <dispnode.h>
#include <gennodelist.h>

///////////////////////////////////////////////////////////////////////////////
// class DisplayNodeList
//	This will be used to represent the display state lists.
//////////////////////////////////////////////////////////////////////////////

class DisplayNodeList : public GenNodeList
{
public:
    DisplayNodeList(displayStateEnum type);
    ~DisplayNodeList() { }

// ADDED functions
    displayStateEnum GetState()	{ return listType; }

// REDEFINED functions
		// deletes node from its previous List & appends
    virtual void Append(GenericNode *node);
		// deletes newNode from its previous List & inserts in
		//	relation to existNode
    virtual void InsertAfter(GenericNode *newNode, GenericNode *existNode);
    virtual void InsertBefore(GenericNode *newNode, GenericNode *existNode);
    virtual void Remove(GenericNode *node);

protected:
private:
    displayStateEnum listType;
};

//////////////////////////////////////////////////////////////////////////////
// class DisplayNodeList inline functions
// these functions don't rely on any inline functions (its own or
//	other classes) that aren't in this file except for Generic functions
//////////////////////////////////////////////////////////////////////////////

inline DisplayNodeList::DisplayNodeList(displayStateEnum type) 
	: GenNodeList(new DisplayNode())
{
    listType = type;
    ((DisplayNode *)head)->displayState = type;
}

inline void DisplayNodeList::Remove(GenericNode *node)
{
    GenNodeList::Remove(node);
// DON'T DO THIS    ((DisplayNode *)node)->displayState = noMapState;
}

#endif
