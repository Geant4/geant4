

//



//
// $Id: gennode.h,v 1.3 1999-12-15 18:04:16 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef gennode_h
#define gennode_h

/*
* NIST Utils Class Library
* clutils/gennode.h
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

#include "G4ios.hh"

class GenNodeList;
class MgrNodeList;
class DisplayNodeList;

//////////////////////////////////////////////////////////////////////////////
// GenericNode
// If you delete this object it first removes itself from any List it is in.
//////////////////////////////////////////////////////////////////////////////

class GenericNode
{
friend class GenNodeList;
friend class MgrNodeList;
friend class DisplayNodeList;

protected:
    GenericNode *next;
    GenericNode *prev;
public:
    GenericNode()	{ next = 0; prev = 0; }
    virtual ~GenericNode()	{ Remove(); }
    GenericNode *Next()	{ return next; }
    GenericNode *Prev()	{ return prev; }
    virtual void Append(GenNodeList *List);
    virtual void Remove()
    {
	(next) ? (next->prev = prev) : 0;
	(prev) ? (prev->next = next) : 0;
/*
//	if(next)
//	    next->prev = prev;
//	if(prev)
//	    prev->next = next;
*/
	next = 0;
	prev = 0;

    }
};

//////////////////////////////////////////////////////////////////////////////
// GenericNode inline functions
// these functions don't rely on any inline functions (its own or
//	other classes) that aren't in this file
//////////////////////////////////////////////////////////////////////////////

#endif
