

//



//
// $Id: gennode_inline.cc,v 1.2 1999-05-21 20:21:08 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef gennode_inline_cc
#define gennode_inline_cc

/*
* NIST Utils Class Library
* clutils/gennode_inline.cc
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

#include <gennode.h>
#include <gennodelist.h>

//////////////////////////////////////////////////////////////////////////////
// class GenericNode inline functions that depend on other classes
// 	inline functions
//
// depends on:
//	void GenNodeList::Append(GenericNode *node) from the gennodelist.h
//////////////////////////////////////////////////////////////////////////////
void GenericNode::Append(GenNodeList *List)
{
//    if(debug_level >= PrintFunctionTrace)
//	G4cout << "GenericNode::Append()\n";
//    if(debug_level >= PrintValues)
//	G4cout << "GenericNode::this : '" << this << "'\n";
    List->Append(this);
}

#endif
