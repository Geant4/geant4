
/*
* NIST Utils Class Library
* clutils/gennode.cc
* September 1996
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/* $Id: gennode.cc,v 1.1 2000-01-21 13:43:10 gcosmo Exp $  */ 

#include <gennode.h>
#include <gennodelist.h>

#ifdef __O3DB__
#include <OpenOODB.h>
#endif
//////////////////////////////////////////////////////////////////////////////
// class GenericNode inline functions that depend on other classes
// 	inline functions
//
// depends on:
//	void GenNodeList::Append(GenericNode *node) from the gennodelist.h
//////////////////////////////////////////////////////////////////////////////

void GenericNode::Append(GenNodeList *list)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "GenericNode::Append()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "GenericNode::this : '" << this << "'\n";
    list->Append(this);
}

