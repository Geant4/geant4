

//



//
// $Id: instmgr.h,v 1.2 1999-05-21 20:20:41 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef instmgr_h
#define instmgr_h

/*
* NIST STEP Editor Class Library
* cleditor/instmgr.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

///// future? ////////////////
// InstMgr can maintain an undo List for the last operation
//	performed on a node 
// InstMgr can have a startUndo() and endUndo() so it knows when to
//	start a new undo List and delete the old undo List. 
/////////////////////

#ifdef __O3DB__
#include <OpenOODB.h>
#endif



extern char * EntityClassName ( char *);

// IT IS VERY IMPORTANT THAT THE ORDER OF THE FOLLOWING INCLUDE FILES
// BE PRESERVED

#include <gennode.h>
#include <gennodelist.h>
#include <gennodearray.h>

#include <mgrnode.h>
#include <mgrnodelist.h>

#include <dispnode.h>
#include <dispnodelist.h>

#include <mgrnodearray.h>

class InstMgr
{
protected:
    MgrNodeArray *master;	// master array of all MgrNodes made up of
			// complete, incomplete, new, delete MgrNodes lists
			// this corresponds to the display List object by index
    MgrNodeArraySorted *sortedMaster;	// master array sorted by fileId
//    StateList *master; // this will be an sorted array of ptrs to MgrNodes

    int maxFileId;

public:
    InstMgr();
    ~InstMgr() {};

// MASTER LIST OPERATIONS
    int InstanceCount()   		{ return master->Count(); }
    void ClearInstances();
    Severity VerifyInstances(ErrorDescriptor& e);

    // DAS PORT possible BUG two funct's below may create a temp for the cast
    MgrNode *&operator[](int index)   
	{ return (MgrNode* &)((*master)[index]); }
    MgrNode *&GetMgrNode(int index)   
	{ return (MgrNode* &)((*master)[index]); }

    MgrNode *FindFileId(int fileId);
	// get the index into display List given a STEPentity
	//  called by see initiated functions
    int GetIndex(STEPentity *se);
    int GetIndex(MgrNode *mn);

// L. Broglia : new
    int GetIndex(const char* entityKeyword, int starting_index);

    int VerifyEntity(int fileId, const char *expectedType);

//    void Append(MgrNode *node);
    MgrNode *Append(STEPentity *se, stateEnum listState);
		// deletes node from master List structure 
    void Delete(MgrNode *node);
    void Delete(STEPentity *se);

    void ChangeState(MgrNode *node, stateEnum listState);

    int MaxFileId()		{ return maxFileId; }
    int NextFileId()		{ return maxFileId = maxFileId +1; }
    int EntityKeywordCount(const char* name);
    STEPentity *GetSTEPentity(int index);
    STEPentity *GetSTEPentity(const char* entityKeyword, 
			      int starting_index =0);

    STEPentity *GetSTEPentity(MgrNode *node) 
					{ return node->GetSTEPentity(); };
    void *GetSEE(int index);
    void *GetSEE(MgrNode *node) { return node->SEE(); };

    void PrintSortedFileIds();

};

#endif
