// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: seeinfodefault.h,v 1.1 1999-01-07 16:08:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef seeinfodefault_h
#define seeinfodefault_h

/*
* NIST STEP Editor Class Library
* cleditor/schemaheader.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

// this is a default seeinfo that does nothing... thus it is not
//	dependent on a user interface toolkit

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

class MgrNode;
class DisplayNode;
class DisplayNodelist;
class STEPentity;

#include <editordefines.h>

class seeInfo : public DisplayNode
{
public:
    seeInfo(MgrNode *node, 
	    STEPentity *se,
	    DisplayNodeList *dnl, displayStateEnum displaySt = mappedWrite);

    void *GetSEE()		{ return see; }
};

inline seeInfo::seeInfo(MgrNode *node, STEPentity *se,
		 DisplayNodeList *dnl, displayStateEnum displaySt)
{
    mn = node; 
    see = 0;
    displayState = displaySt; 
    dnl->Append(this);
}

#endif
