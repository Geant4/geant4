// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: mgrnode.h,v 1.1 1999-01-07 16:08:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef mgrnode_h
#define mgrnode_h

/*
* NIST STEP Editor Class Library
* cleditor/mgrnode.h
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

class DisplayNode;
class STEPentity;

#include <gennode.h>
#include <gennodelist.h>

#include <editordefines.h>

class InstMgr;

//////////////////////////////////////////////////////////////////////////////
// class MgrNode
// If you delete this object, it deletes the STEPentity it represents,
// the DisplayNode, and removes itself from any List it is in.
//////////////////////////////////////////////////////////////////////////////

class MgrNode : public GenericNode
{
    friend class GenNodeList;
    friend class MgrNodeList;
    friend class InstMgr;
    
protected:
	// currState, next, prev implement several lists
	// based on currState:
	// currState = (completeSE, incompleteSE, newSE, or deleteSE)
	// every node will be on one of the four lists implemented by these:
    stateEnum currState;

	// STEPentity this node is representing info for
    STEPentity *se;
	// this is the index (in the InstMgr master array) of the ptr to
	//   this node.
    int arrayIndex;

	// display info (SEE, etc) for this node
    DisplayNode *di;

public:
	// used for sentinel node on lists of MgrNodes
    MgrNode();
    MgrNode(STEPentity *se);
	// 'listState' ==
	//	completeSE - if reading valid exchange file
	//	incompleteSE or completeSE - if reading working session file
	//	newSE - if instance is created by user using editor (probe)
    MgrNode(STEPentity *se, stateEnum listState);
    MgrNode(STEPentity *se, stateEnum listState, MgrNodeList *List);
    virtual ~MgrNode();

// STATE LIST OPERATIONS
    int MgrNodeListMember(stateEnum s)	{ return (currState == s); }
    stateEnum  CurrState()		{ return  currState;	    }
		// returns next or prev member variables
		// i.e. next or previous node on curr state List
			// searches current List for fileId
    MgrNode *StateFindFileId(int fileId);

			// deletes from previous cmd List, 
			//	& puts on cmd List cmdList
    int ChangeList(MgrNodeList *cmdList);
    int ChangeState(stateEnum s);

	// Removes from current List. 
	// 	Called before adding to diff List or when destructor is called.
    void Remove();

// DISPLAY LIST OPERATIONS
    void *SEE();

    displayStateEnum DisplayState();
    int IsDisplayState(displayStateEnum ds);

		// returns next or prev member variables
		// i.e. next or previous node on display state List
    GenericNode *NextDisplay();
    GenericNode *PrevDisplay();

			// deletes from previous cmd List, 
			//	& puts on cmd List cmdList
    int ChangeList(DisplayNodeList *cmdList);
			// deletes from previous display List, assigns ds to 
			//	displayState & puts on List dsList
    int ChangeState(displayStateEnum ds);

// 	might not want these three? since it won't actually map them?
    void MapModifiable(DisplayNodeList *dnList);
    void MapViewable(DisplayNodeList *dnList);
    void UnMap(DisplayNodeList *dnList);

// ACCESS FUNCTIONS
    int GetFileId();
    STEPentity *GetSTEPentity()	{ return se; }
    DisplayNode *&displayNode() { return di; }
    int ArrayIndex()		{ return arrayIndex; }
    void ArrayIndex(int index)	{ arrayIndex = index; }

protected:

private:
    void Init(STEPentity *s, stateEnum listState, MgrNodeList *List);
};

//////////////////////////////////////////////////////////////////////////////
// class MgrNode inline functions
// these functions don't rely on any inline functions (its own or
//	other classes) that aren't in this file except for Generic functions
//////////////////////////////////////////////////////////////////////////////

#endif
