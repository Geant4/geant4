

//



//
// $Id: mgrnodelist.h,v 1.2 1999-05-21 20:20:42 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef mgrnodelist_h
#define mgrnodelist_h

/*
* NIST STEP Editor Class Library
* cleditor/mgrnodelist.h
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
#include <gennodelist.h>
#include <editordefines.h>

//////////////////////////////////////////////////////////////////////////////
// class MgrNodeList
//	This will be used to represent the state lists.
//////////////////////////////////////////////////////////////////////////////

class MgrNode;

class MgrNodeList : public GenNodeList
{
public:
    MgrNodeList(stateEnum type);
    virtual ~MgrNodeList() { }

// ADDED functions
    virtual MgrNode *FindFileId(int fileId);

// REDEFINED functions
		// deletes node from its previous List & appends
    virtual void Append(GenericNode *node);
		// deletes newNode from its previous List & inserts in
		//	relation to existNode
    virtual void InsertAfter(GenericNode *newNode, GenericNode *existNode);
    virtual void InsertBefore(GenericNode *newNode, GenericNode *existNode);
    virtual void Remove(GenericNode *node);

protected:
    stateEnum listType;
};

//////////////////////////////////////////////////////////////////////////////
// class MgrNodeList inline functions
// these functions don't rely on any inline functions (its own or
//	other classes) that aren't in this file except for Generic functions
//////////////////////////////////////////////////////////////////////////////

#endif
