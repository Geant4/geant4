// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: dispnode.cc,v 1.1 1999-01-07 16:08:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Editor Class Library
* cleditor/dispnode.cc
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

#include <gennode.h>
#include <gennodelist.h>

#include <dispnode.h>
#include <dispnodelist.h>

// define this to be the name of the display object
class StepEntityEditor;

// This function needs to be defined outside the SCL libraries.  It needs to do
// two things:
// 1) unmap the StepEntityEditor window if it is mapped.
// 2) delete the StepEntityEditor window
// To see an example of this function used with the Data Probe look in
// ../clprobe-ui/StepEntEditor.cc  Look at DeleteSEE() and ~StepEntityEditor().

// Note: the following function has been redefined for porting on DEC/OSF1 and SUN/OS
// extern void DeleteSEE(StepEntityEditor *se);
void DeleteSEE(StepEntityEditor *se) {
  if (se) se=NULL;
}

DisplayNode::~DisplayNode()
{
    Remove();
    if(see)
    {
	DeleteSEE((StepEntityEditor *)see);
//DAS PORT need the cast from void*	DeleteSEE(see);
    }
}

void DisplayNode::Remove()
{
    GenericNode::Remove();
// DON'T DO THIS!!    displayState = noMapState;
}

int DisplayNode::ChangeState(displayStateEnum s)
{
    displayState = s;
    return 1;
}

int DisplayNode::ChangeList(DisplayNodeList *cmdList)
{
    Remove();
    cmdList->Append(this);
    return 1;
}
