

//



//
// $Id: cmdmgr.cc,v 1.2 1999-05-21 20:21:06 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/cmdmgr.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

#include <cmdmgr.h>


ReplicateLinkNode *ReplicateList::FindNode(MgrNode *mn)
{
    ReplicateLinkNode *rln = (ReplicateLinkNode *)GetHead();
    int numEntries = EntryCount();
    while(numEntries--)
    {
	if(rln->ReplicateNode() == mn)
	{
	    return rln;
	}
	rln = (ReplicateLinkNode *)rln->NextNode();
    }
    return 0;
}

BOOLEAN ReplicateList::IsOnList(MgrNode *mn)
{
    return (FindNode(mn) != 0);
/*
    ReplicateLinkNode *rln = (ReplicateLinkNode *)GetHead();
    int numEntries = EntryCount();
    int found = 0;
    while(numEntries--)
    {
	if(rln->ReplicateNode() == mn)
	{
	    found = 1;
	    numEntries = 0;
	}
	rln = (ReplicateLinkNode *)rln->NextNode();
    }
    return found;
*/
}

///////////////////////////////////////////////////////////////////////////////
// returns true if it could delete the node
///////////////////////////////////////////////////////////////////////////////
BOOLEAN ReplicateList::Remove(ReplicateLinkNode *rln)
{
    ReplicateLinkNode *rnFollow = (ReplicateLinkNode *)GetHead();
    if(!rnFollow || !rln)
	return 0;
    else
    {
	if(rnFollow == rln)
	{
	    head = rln->NextNode();
	    delete rln;
	    return 1;
	}
	else
	{
	    ReplicateLinkNode *rn = (ReplicateLinkNode *)rnFollow->NextNode();
	    while(rn)
	    {
		if(rn == rln)
		{
		    rnFollow->next = (SingleLinkNode *)rln->NextNode();
		    delete rln;
		    return 1;
		}
		rnFollow = rn;
		rn = (ReplicateLinkNode *)rn->NextNode();
	    } // end while(rn)
	} // end else
    } // end else
    return 0;
}

BOOLEAN ReplicateList::Remove(MgrNode *rn)
{
    return Remove(FindNode(rn));
}

CmdMgr::CmdMgr()
{
    completeList = new MgrNodeList(completeSE);
    incompleteList = new MgrNodeList(incompleteSE);
//    newList = new MgrNodeList(newSE);
    deleteList = new MgrNodeList(deleteSE);

    mappedWriteList = new DisplayNodeList(mappedWrite);
    mappedViewList = new DisplayNodeList(mappedView);
    closeList = new DisplayNodeList(notMapped);
    replicateList = new ReplicateList();
}

int CmdMgr::ReplicateCmdList(MgrNode *mn)
{
    if(!(replicateList->IsOnList(mn)))
    {
	replicateList->AddNode(mn);
    }
    return 1;
}

/*
void CmdMgr::ModifyCmdList(MgrNode *mn)
{
    mn->ChangeList(mappedWriteList);
}

void CmdMgr::ViewCmdList(MgrNode *mn)
{
    mn->ChangeList(mappedViewList);
}

void CmdMgr::CloseCmdList(MgrNode *mn)
{
    mn->ChangeList(closeList);
}
*/

void CmdMgr::ClearInstances()
{
    completeList->ClearEntries();
    incompleteList->ClearEntries();
    cancelList->ClearEntries();
    deleteList->ClearEntries();
    replicateList->Empty();

//    newList->ClearEntries();
}
			// searches current List for fileId
MgrNode *CmdMgr::StateFindFileId(stateEnum s, int fileId)
{
    switch(s){
	case completeSE:
		return completeList->FindFileId(fileId);
	case incompleteSE:
		return incompleteList->FindFileId(fileId);
	case deleteSE:
		return deleteList->FindFileId(fileId);
	case newSE: // there is no new List
	case noStateSE:
	default:
		G4cout << "ERROR can't find the node containing fileid " <<
			fileId << " from this node\n";
		return 0;
    }
}

MgrNode *CmdMgr::GetHead(stateEnum listType)
{
    switch(listType)
    {
	case completeSE:	// saved complete List
		return (MgrNode *)completeList->GetHead();
	case incompleteSE:	// saved incomplete List
		return (MgrNode *)incompleteList->GetHead();
	case deleteSE:		// delete List
		return (MgrNode *)deleteList->GetHead();
	default:
		return 0;
    }
}

DisplayNode *CmdMgr::GetHead(displayStateEnum listType)
{
    switch(listType)
    {
	case mappedWrite:
		return (DisplayNode *)mappedWriteList->GetHead();

	case mappedView:
		return (DisplayNode *)mappedViewList->GetHead();

	case notMapped:
		return (DisplayNode *)closeList->GetHead();

	case noMapState:
	default:
		return 0;
    }
}

void CmdMgr::ClearEntries(stateEnum listType)
{
    switch(listType)
    {
	case completeSE:	// saved complete List
		completeList->ClearEntries();
		break;
	case incompleteSE:	// saved incomplete List
		incompleteList->ClearEntries();
		break;
	case deleteSE:		// delete List
		deleteList->ClearEntries();
		break;
	default:
		break;
    }
}

void CmdMgr::ClearEntries(displayStateEnum listType)
{
    switch(listType)
    {
	case mappedWrite:
		mappedWriteList->ClearEntries();
		break;
	case mappedView:
		mappedViewList->ClearEntries();
		break;
	case notMapped:
		closeList->ClearEntries();
		break;
	case noMapState:
	default:
		break;
    }
}
