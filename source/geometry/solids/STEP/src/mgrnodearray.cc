

//



//
// $Id: mgrnodearray.cc,v 1.2 1999-05-21 20:21:10 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/mgrnodearray.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

/*
 * MgrNodeArray - dynamic array object of MgrNodes.
 * the array part of this is copied from Unidraws UArray - dynamic array object
 * Copyright (c) 1990 Stanford University
 */

///////////////////////////////////////////////////////////////////////////////
//	debug_level >= 2 => tells when a command is chosen
//	debug_level >= 3 => prints values within functions:
//	   e.g. 1) entity type List when read
//		2) entity instance List when read
///////////////////////////////////////////////////////////////////////////////
static int debug_level = 1;
	// if debug_level is greater than this number then
	// function names get printed out when executed
static int PrintFunctionTrace = 2;
	// if debug_level is greater than this number then
	// values within functions get printed out
//static int PrintValues = 3;

#include <mgrnodearray.h>
#include <STEPentity.h>

#include <string.h>	// to get bcopy() - ANSI


//////////////////////////////////////////////////////////////////////////////
// class MgrNodeArray member functions
//////////////////////////////////////////////////////////////////////////////

void MgrNodeArray::AssignIndexAddress(int index)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNodeArray::AssignIndexAddress()\n";
    ((MgrNode *)_buf[index])->ArrayIndex(index);
}

MgrNodeArray::~MgrNodeArray()
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArray::~MgrNodeArray()\n";
    DeleteEntries();
}

/*****************************************************************************/

void MgrNodeArray::ClearEntries()
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArray::ClearEntries()\n";
    int i;
    for(i = 0 ; i < _count; i++)
	_buf[i] = 0;
    _count = 0;
}

/*****************************************************************************/

void MgrNodeArray::DeleteEntries()
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArray::DeleteEntries()\n";
    int i;
    for(i = 0 ; i < _count; i++)
	delete ((MgrNode *)_buf[i]);
    _count = 0;
}

/*****************************************************************************/

int MgrNodeArray::Insert(GenericNode* gn, int index)
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArray::Insert()\n";
    int AssignedIndex = GenNodeArray::Insert(gn, index);
    int i;
    for(i = AssignedIndex ; i < _count; i++)
	((MgrNode *)_buf[i])->ArrayIndex(i);
    return AssignedIndex;
}

/*****************************************************************************/

void MgrNodeArray::Remove(int index) {
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArray::Remove()\n";
    if (0 <= index && index < _count) {
	GenNodeArray::Remove(index);
	int i;
	for(i = index; i < _count; i++)
	    ((MgrNode *)_buf[i])->ArrayIndex(i);
    }
}

/*****************************************************************************/

int MgrNodeArray::MgrNodeIndex(int fileId)
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArray::MgrNodeIndex()\n";
    int i;
    for (i = 0; i < _count; ++i) {
        if( ((MgrNode *)_buf[i])->GetSTEPentity()->GetFileId() == fileId) {
            return i;
        }
    }
    return -1;
}

//////////////////////////////////////////////////////////////////////////////
// class MgrNodeArraySorted member functions
//////////////////////////////////////////////////////////////////////////////

int MgrNodeArraySorted::Insert (GenericNode* gn) {
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNodeArraySorted::Insert()\n";

	// since gn is really a MgrNode
    int fileId = ((MgrNode *)gn)->GetSTEPentity()->GetFileId();

    int index = FindInsertPosition(fileId);

    return GenNodeArray::Insert(gn, index);
}

int MgrNodeArraySorted::Index (GenericNode* gn) {
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNodeArraySorted::Index()\n";
	// since gn is really a MgrNode
    return MgrNodeIndex(((MgrNode *)gn)->GetFileId());
}

int MgrNodeArraySorted::Index (GenericNode** gn) {
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "MgrNodeArraySorted::Index()\n";
	// since gn is really a MgrNode
    return MgrNodeIndex(((MgrNode *)(*gn))->GetFileId());
}

void MgrNodeArraySorted::ClearEntries()
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArraySorted::ClearEntries()\n";
    int i;
    for(i = 0 ; i < _count; i++)
	_buf[i] = 0;
    _count = 0;
}

/*****************************************************************************/

void MgrNodeArraySorted::DeleteEntries()
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArraySorted::DeleteEntries()\n";
    int i;
    for(i = 0 ; i < _count; i++)
	delete ((MgrNode *)_buf[i]);
    _count = 0;
}

/*****************************************************************************/

	// the reason this is written this way is because most of the
	// time the file id will be higher than any seen so far and
	// thus the insert position will be at the end
int MgrNodeArraySorted::FindInsertPosition(const int fileId)
{
    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArraySorted::FindInsertPosition()\n";
    int i;
    int curFileId;

    for (i = _count-1; i >= 0; --i) {
	curFileId = ((MgrNode *)_buf[i])->GetSTEPentity()->GetFileId();
        if (curFileId < fileId /*|| curFileId == fileId*/)
	{
            return i + 1;
        }
    }
    return 0;
}

/*****************************************************************************/

int MgrNodeArraySorted::MgrNodeIndex(int fileId) {
// this function assumes that _buf[0] to _buf[_count] ALL point to MgrNodes
// that are sorted by fileId

    if(debug_level >= PrintFunctionTrace)
	G4cout << "MgrNodeArraySorted::MgrNodeIndex()\n";
    int low = 0;
    int high = _count - 1;
    int mid;
    int found = 0;
    int curFileId;

    while(!found && (low <= high))
    {
	mid = (low + high) / 2;
	curFileId = ((MgrNode *)_buf[mid])->GetSTEPentity()->GetFileId();
	if(curFileId == fileId)
	{
	    found = 1;
	}
	else if(curFileId < fileId)
	{
	    low = mid + 1;
	}
	else
	{
	    high = mid - 1;
	}
    }
    if(found)
	return mid;
    return -1;
}
