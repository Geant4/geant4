// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: SingleLinkList.cc,v 2.0 1998/07/02 17:04:50 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//

/*
* NIST STEP Core Class Library
* clstepcore/SingleLinkList.cc
* May 1995
* David Sauder
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <SingleLinkList.h>

void
SingleLinkList::AppendNode (SingleLinkNode * item)  
{
    if (head)  {
	tail -> next = item;
	tail = item;
    }
    else head = tail = item;
}

