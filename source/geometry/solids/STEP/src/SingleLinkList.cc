

//



//
// $Id: SingleLinkList.cc,v 1.2 1999-05-21 20:21:05 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

