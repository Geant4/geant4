

//



//
// $Id: instmgr.cc,v 1.3 1999-05-21 20:21:09 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/instmgr.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

//////////////////////////////////////////////////////////////////////////////
//
// InstMgr member functions
//
//////////////////////////////////////////////////////////////////////////////

#include <STEPentity.h>
#include <instmgr.h>

///////////////////////////////////////////////////////////////////////////////
//	debug_level >= 2 => tells when a command is chosen
//	debug_level >= 3 => prints values within functions:
///////////////////////////////////////////////////////////////////////////////
static int debug_level = 2;
	// if debug_level is greater than this number then
	// function names get printed out when executed
//static int PrintFunctionTrace = 2;
	// if debug_level is greater than this number then
	// values within functions get printed out
//static int PrintValues = 3;


///////////////////////////////////////////////////////////////////////////////

void 
InstMgr::PrintSortedFileIds()
{
    MgrNode *mn = 0;
    int count = InstanceCount();
    int i = 0;
    for(i = 0; i < count; i++)
      {
	mn = (MgrNode *)((*sortedMaster)[i]);
	G4cout << i << " " << mn->GetFileId() << endl;
      }
}

InstMgr::InstMgr()
{
    maxFileId = -1;

    master = new MgrNodeArray();
    sortedMaster = new MgrNodeArraySorted();
}

///////////////////////////////////////////////////////////////////////////////

void InstMgr::ClearInstances()
{
//    delete master;
//    master = new MgrNodeArray();
    master->DeleteEntries();
    sortedMaster->ClearEntries();
    maxFileId = -1;
}

///////////////////////////////////////////////////////////////////////////////

/**************************************************
 description:
 returns:
    SEVERITY_NULL:        if all instances are complete
    SEVERITY_INCOMPLETE:  if at least one instance is missing a required attribute
    SEVERITY_WARNING:     
    SEVERITY_INPUT_ERROR: if at least one instance can not be fetched
                          from the instance List.
**************************************************/

enum Severity
InstMgr::VerifyInstances(ErrorDescriptor& err) 
{
    int errorCount = 0;
    char errbuf[BUFSIZ];
    
    int n = InstanceCount();
    MgrNode* mn;
    STEPentity* se;
    enum Severity rval = SEVERITY_NULL;
    
    //for each instance on the List,
    //check the MgrNode state.
    //if the state is complete then do nothing
    //if the state is not complete,
    //   then check if it is valid
    //   if it is valid then Increment the warning count,
    //      and set the rval to SEVERITY_INCOMPLETE
    //   if it is not valid, then Increment the error count 
    //      and set the rval to 

    for (int i = 0; i < n; ++i) 
    {
	mn = GetMgrNode(i);
	if (!mn) 
	{
	    ++errorCount;
	    if (errorCount == 1) 
		sprintf(errbuf,
		"VerifyInstances: Unable to verify the following instances: node %d", 
			i);
	    else
		sprintf(errbuf,", node %d",i);

	    err.AppendToDetailMsg(errbuf);
	    rval = SEVERITY_INPUT_ERROR;
	    err.GreaterSeverity(SEVERITY_INPUT_ERROR);
	    continue;
	}
	if (!mn->MgrNodeListMember(completeSE)) 
	{
	    se = mn->GetSTEPentity();
	    if (se->ValidLevel(&err,this,0) < SEVERITY_USERMSG)
	    {
		if (rval > SEVERITY_INCOMPLETE) 
		    rval = SEVERITY_INCOMPLETE;
		++errorCount;
		if (errorCount == 1)
		    sprintf(errbuf,
			    "VerifyInstances: Unable to verify the following instances: #%d",
			    se->FileId());
		else 
		    sprintf(errbuf,", #%d",se->FileId());
		err.AppendToDetailMsg(errbuf);
	    }
	}
    }
    if (errorCount) 
    {
	sprintf(errbuf,
		"VerifyInstances: %d invalid instances in List.\n", 
		errorCount);
	err.AppendToUserMsg(errbuf);
	err.AppendToDetailMsg(".\n");
	err.GreaterSeverity(SEVERITY_INCOMPLETE);
    }

    return rval;
}

///////////////////////////////////////////////////////////////////////////////

MgrNode *InstMgr::FindFileId(int fileId)
{
    int index = sortedMaster->MgrNodeIndex(fileId);
    if(index >= 0)
	return (MgrNode *)(*sortedMaster)[index];
    else
	return (MgrNode *)0;
}

///////////////////////////////////////////////////////////////////////////////

	// get the index into display List given a STEPentity
	//  called by see initiated functions
int InstMgr::GetIndex(MgrNode *mn)
{
    return mn->ArrayIndex();
}

///////////////////////////////////////////////////////////////////////////////

int InstMgr::GetIndex(STEPentity *se)
{
    int fileId = se->FileId();
    return sortedMaster->MgrNodeIndex(fileId);
}

///////////////////////////////////////////////////////////////////////////////

int InstMgr::VerifyEntity(int fileId, const char *expectedType)
{
    MgrNode *mn = FindFileId(fileId);
    if(mn) 
    {
	if(!strcmp(expectedType, mn->GetSTEPentity()->EntityName()))
	    return 2;	// types match
	else
	    return 1;	// possible mismatch depending on descendants
    }
    else
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
//   Append instance to the List of instances.  Checks the file id and 
//   sets it if 1) it is not set already or 2) it already exists in the List.

MgrNode *InstMgr::Append(STEPentity *se, stateEnum listState)
{
    if (debug_level > 3)
	G4cout << "#" << se->FileId() << " append node to InstMgr\n";

    if ((se  -> FileId () == 0)  // no id assigned
	|| FindFileId (se -> FileId ())) // id already in List
      se -> FileId ( NextFileId () );
    else if (se->FileId() > MaxFileId()) maxFileId = se->FileId();

    MgrNode *mn = new MgrNode(se, listState);
    if(listState == noStateSE)
	G4cout << "append to InstMgr **ERROR ** node #" << se->FileId() << 
		" doesn't have state information\n";
    master->Append(mn);
    sortedMaster->Insert(mn);
//PrintSortedFileIds();
    return mn;
}

///////////////////////////////////////////////////////////////////////////////

void InstMgr::Delete(MgrNode *node)
{
	// delete the node from its current state List
    node->Remove();

    int index;    

	// get the index of the node in the sorted master array
    index = sortedMaster->MgrNodeIndex(node->GetFileId());
	// remove the node from the sorted master array
    sortedMaster->Remove(index);

	// get the index into the master array by ptr arithmetic
    index = node->ArrayIndex();
    master->Remove(index);

    delete node;
}

///////////////////////////////////////////////////////////////////////////////

void InstMgr::Delete(STEPentity *se)
{
   Delete (FindFileId (se->FileId()));
}

///////////////////////////////////////////////////////////////////////////////

void InstMgr::ChangeState(MgrNode *node, stateEnum listState)
{
    switch(listState){
	case completeSE:
		if (debug_level > 3)
		    G4cout << "#" << node->GetSTEPentity()->FileId() << 
		    " change node to InstMgr's complete List\n";
		node->ChangeState(listState);
		break;
	case incompleteSE:
		if (debug_level > 3)
		    G4cout << "#" << node->GetSTEPentity()->FileId() << 
		    " change node to InstMgr's incomplete List\n";
		node->ChangeState(listState);
		break;
	case newSE:
		if (debug_level > 3)
		    G4cout << "#" << node->GetSTEPentity()->FileId() << 
		    " change node to InstMgr's new List\n";
		node->ChangeState(listState);
		break;
	case deleteSE:
		if (debug_level > 3)
		    G4cout << "#" << node->GetSTEPentity()->FileId() << 
		    " change node to InstMgr's delete List\n";
		node->ChangeState(listState);
		break;
	case noStateSE:
		G4cout << "#" << node->GetSTEPentity()->FileId() << 
		"ERROR can't change this node state\n";
		node->Remove();
		break;
    }
}

///////////////////////////////////////////////////////////////////////////////

/**************************************************
 description:
    This function returns an integer value indicating
    the number of instances with the given name appearing
    on the instance manager.
**************************************************/
int 
InstMgr::EntityKeywordCount(const char* name) 
{
    int count = 0;
    MgrNode* node;
    STEPentity* se;
    int n = InstanceCount();
    for (int j = 0; j<n; ++j) 
      {
	  node = GetMgrNode(j);
	  se = node->GetSTEPentity();
	  if (!strcmp(se->EntityName(),
		      PrettyTmpName(name)))
	      ++count;
      }
    return count;
}

///////////////////////////////////////////////////////////////////////////////

STEPentity *
InstMgr::GetSTEPentity(int index)
{
    MgrNode *mn = (MgrNode*)(*master)[index];
    if(mn)
	return mn->GetSTEPentity(); 
    else
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

/**************************************************
 description:
    This function returns the first entity (by index)
    on the instance manager, which has the given
    entity name. It starts its search at starting_index,
    and returns ENTITY_NULL if a match does not occur 
    before the last index is reached. This function 
    does not wrap around to search indices before the
    starting_index.
**************************************************/
STEPentity*
InstMgr::GetSTEPentity(const char* entityKeyword, int starting_index) 
{
  // L. Broglia
  G4cout<<"Warning ! This function don`t work correctly...";

  MgrNode *node;
  STEPentity *se;
  
  int count = InstanceCount();
  for (int j = starting_index; j<count; ++j) 
  {
    node = GetMgrNode(j);
    se = node->GetSTEPentity();
   
    if ( ! strcmp( se->EntityName(),PrettyTmpName(entityKeyword)) )
    {
      return se;
    }
  }
  
  return ENTITY_NULL;
}
////////////////////////////////////////////////////////
// L. Broglia : new

int InstMgr::GetIndex(const char* entityKeyword, int starting_index)
{
  // L. Broglia
  // When the entity is founded, put return index

  MgrNode *node;
  STEPentity *se;
  
  int count = InstanceCount();
  const char* tmpname;

  for (int j = starting_index; j<count; ++j) 
  {
    node = GetMgrNode(j);
    se = node->GetSTEPentity();
    tmpname = se->EntityName();

    if ( !strcmp(tmpname, entityKeyword) )
    {
      return j;
    }
  }
  return 0;
}



///////////////////////////////////////////////////////////////////////////////

void *
InstMgr::GetSEE(int index)
{
    MgrNode *mn = (MgrNode*)(*master)[index];
    if(mn)
	return mn->SEE(); 
    else
	return 0;
}
