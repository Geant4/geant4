

//



//
// $Id: mgrnode.cc,v 1.3 1999-12-15 18:04:21 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/mgrnode.cc
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

#include <instmgr.h>
#include <STEPentity.h>

#include "G4ios.hh"

void *MgrNode::SEE()
{
    return (di ? di->SEE() : 0);
}

int MgrNode::GetFileId()
{
    return (se ? se->GetFileId() : -1);
}

void MgrNode::Remove()
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::Remove()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "MgrNode::this : '" << this << "'\n";
    GenericNode::Remove();
// DON'T DO THIS!!    currState = noStateSE;
}

	// searches current List for fileId
MgrNode *MgrNode::StateFindFileId(int fileId)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::StateFindFileId()\n";
    MgrNode *startNode = this;
    if(startNode->GetFileId() == fileId) return this;
    else
    {
		// mn is really a MgrNode
	MgrNode *mn = (MgrNode *)(startNode->Next());
	while(mn != startNode)
	{
	    if( mn->GetFileId() == fileId)
		return (MgrNode *)mn;
	    mn = ((MgrNode *)mn->Next());
	}
	return (MgrNode *)0;
    }
}

MgrNode::~MgrNode()
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::~MgrNode()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "MgrNode::this : '" << this << "'\n";
    if(se)
	delete se;
    if(di)
	delete di;
//    GenericNode::Remove(); // this is called by default.
}

///////////////////// class MgrNode Display Functions /////////////////////////

displayStateEnum MgrNode::DisplayState() 
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::DisplayState()\n";
    return (di ? di->DisplayState() : noMapState);
}

int MgrNode::IsDisplayState(displayStateEnum ds)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::IsDisplayState()\n";
    return (di ? di->DisplayListMember(ds) : 0);
}

GenericNode *MgrNode::NextDisplay()
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::NextDisplay()\n";
//    return (di ? ((DisplayNode *)di->Next()) : (DisplayNode *)0);
    if(di)
    {
//	GenericNode *dn = di->Next();
//	return (DisplayNode *)dn;
//    	return (DisplayNode *)(di->Next());
    	return di->Next();
    }
    else
	return 0;
}

GenericNode *MgrNode::PrevDisplay()
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::PrevDisplay()\n";
//    return (di ? ((DisplayNode *)di->Prev()) : 0);
    if(di)
	return di->Prev();
    else
	return 0;
}

// STATE LIST OPERATIONS

// deletes from previous cmd List & puts on cmd List cmdList
int MgrNode::ChangeList(DisplayNodeList *cmdList)
{
    if(!di)
	di = new class DisplayNode(this);
    return di->ChangeList(cmdList);
}

// deletes from previous cmd List & puts on cmd List cmdList
int MgrNode::ChangeList(MgrNodeList *cmdList)
{
    Remove();
    cmdList->Append(this);
    return 1;
}

int MgrNode::ChangeState(displayStateEnum s)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::ChangeState()\n";
    if(di)
    {
	return di->ChangeState(s);
    }
    return 0;
}

int MgrNode::ChangeState(stateEnum s)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::ChangeState()\n";
    currState = s;
     // for now, later need to type check somehow and return success or failure
    return 1;
}

void MgrNode::Init(STEPentity *s,
			  stateEnum listState,
			  MgrNodeList *List)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::Init()\n";
    se = s;
    arrayIndex = -1;
    di = 0;
    currState = listState;
    if(List)
    {
	List->Append(this);
    }
}

	// used for sentinel node on lists of MgrNodes
MgrNode::MgrNode()
{ 
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::MgrNode()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "MgrNode::this : '" << this << "'\n";
    Init(0, noStateSE, 0);
}

MgrNode::MgrNode(STEPentity *StepEntPtr)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::MgrNode()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "MgrNode::this : '" << this << "'\n";
    Init(StepEntPtr, noStateSE, 0);
}

	// 'listState' ==
	//	completeSE - if reading valid exchange file
	//	incompleteSE or completeSE - if reading working session file
	//	newSE - if instance is created by user using editor (probe)
MgrNode::MgrNode(STEPentity *StepEntPtr, stateEnum listState)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::MgrNode()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "MgrNode::this : '" << this << "'\n";
    Init(StepEntPtr, listState, 0);
}
	// 'listState' ==
	//	completeSE - if reading valid exchange file
	//	incompleteSE or completeSE - if reading working session file
	//	newSE - if instance is created by user using editor (probe)
MgrNode::MgrNode(STEPentity *StepEntPtr, stateEnum listState, MgrNodeList *List)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNode::MgrNode()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "MgrNode::this : '" << this << "'\n";
    Init(StepEntPtr, listState, List);

}
