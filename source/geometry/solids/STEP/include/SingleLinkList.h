// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: SingleLinkList.h,v 1.1 1999-01-07 16:08:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef singlelinklist_h
#define	singlelinklist_h

/*
* NIST STEP Core Class Library
* clstepcore/SingleLinkList.h
* May 1995
* David Sauder
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

class SingleLinkNode {
    friend class SingleLinkList;
  protected:

  public:
    SingleLinkNode* next;

  public:
    virtual SingleLinkNode *NextNode () const;
    SingleLinkNode()  { next =0; }
    virtual ~SingleLinkNode() { }
};

class SingleLinkList  {

    // node which represents the value is contained in the subclass
	//  since it may have different types for different lists
    
  protected:
    
    SingleLinkNode *  head;
    SingleLinkNode *  tail;

  public:
    
    virtual SingleLinkNode *NewNode();
    virtual void AppendNode (SingleLinkNode *);

    virtual void Empty ();
    virtual SingleLinkNode * GetHead () const;
    
    int EntryCount() const;

    SingleLinkList ();
    virtual ~SingleLinkList ();
    
};

#endif
