// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: gennode_inline.cc,v 1.1 1999-01-07 16:08:18 gunter Exp $
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
