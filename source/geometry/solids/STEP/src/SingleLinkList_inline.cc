// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: SingleLinkList_inline.cc,v 1.1 1999-01-07 16:08:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/SingleLinkList.inline.cc
* May 1995
* David Sauder
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <SingleLinkList.h>
#ifdef WIN32
#  include "G4ios.hh"
#else
#  include <stream.h>
#endif

SingleLinkNode * 	
SingleLinkNode::NextNode ()  const
{
    return next;
}

SingleLinkList::SingleLinkList ()  
  : head (0), tail (0)
{
}

SingleLinkList::~SingleLinkList ()  
{
    Empty ();
}

void
SingleLinkList::Empty ()  
{
    SingleLinkNode * tmp = head;
    while (tmp)  
      {
	  tmp = head -> NextNode ();
	  delete head;
	  head = tmp;
      }
}

SingleLinkNode *
SingleLinkList::NewNode () 
{
    //  defined in subtypes
    cerr << "\n\n******BUG****** a virtually defined function should \n"
	 << "be called for SingleLinkList::NewNode()\n\n";
    return new SingleLinkNode();
}

SingleLinkNode *
SingleLinkList::GetHead () const
{
    return (head);
}

int SingleLinkList::EntryCount() const
{
    int entryCount = 0;
    SingleLinkNode *entryPtr = head;

    while( entryPtr != 0 )
    {
	entryPtr = entryPtr->NextNode();
	entryCount++;
    }
    return entryCount;
}
