// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPattribute_inline.cc,v 1.1 1999-01-07 16:08:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/STEPattribute.inline.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <STEPattribute.h>
#include <sdai.h>

//  This is needed so that STEPattribute's can be passed as
//  references to inline functions

STEPattribute::STEPattribute (const STEPattribute& a)
: aDesc (a.aDesc), _derive (0) {}

//  INTEGER
STEPattribute::STEPattribute (const class AttrDescriptor& d, SdaiInteger *p)
: aDesc (&d), _derive (0)
{ ptr.i = p; }

//  BINARY
STEPattribute::STEPattribute (const class AttrDescriptor& d, SdaiBinary *p)
: aDesc (&d), _derive (0)
{ptr.b =p;  }

//  STRING
STEPattribute::STEPattribute (const class AttrDescriptor& d, SdaiString *p)
: aDesc (&d), _derive (0)
{ptr.S =p;  }

//  REAL & NUMBER
STEPattribute::STEPattribute (const class AttrDescriptor& d, SdaiReal *p)
: aDesc (&d), _derive (0)
{ ptr.r = p; }

//  REAL_PTR
/*STEPattribute::STEPattribute (const class AttrDescriptor& d, real **p)
: aDesc (&d), _derive (0)
{ ptr.rp = p; }
*/
//  ENTITY
STEPattribute::STEPattribute (const class AttrDescriptor& d, STEPentity* *p)
: aDesc (&d), _derive (0)
{ ptr.c = p; }

//  AGGREGATE
STEPattribute::STEPattribute (const class AttrDescriptor& d, STEPaggregate *p)
: aDesc (&d), _derive (0)
{ ptr.a =  p; }

//  ENUMERATION  and Logical
STEPattribute::STEPattribute (const class AttrDescriptor& d, STEPenumeration *p)
: aDesc (&d), _derive (0)
{ ptr.e = p;  }

//  SELECT
STEPattribute::STEPattribute (const class AttrDescriptor& d, class SdaiSelect *p)
: aDesc (&d), _derive (0)
{ ptr.sh = p;  }

//  UNDEFINED
STEPattribute::STEPattribute (const class AttrDescriptor& d, SCLundefined *p)
: aDesc (&d), _derive (0)
{ ptr.u = p;  }


const s_String 
STEPattribute::Name() const
	{ return aDesc->Name(); }

const s_String 
STEPattribute::TypeName() const
	{ return aDesc->TypeName(); }

const BASE_TYPE 
STEPattribute::Type() const
{
    return aDesc->Type();
}

const BASE_TYPE 
STEPattribute::NonRefType() const
{ 
    return aDesc->NonRefType();
}

const BASE_TYPE 
STEPattribute::BaseType() const
{ 
    return aDesc->BaseType();
}

/*
const EntityDescriptor *
STEPattribute::ReferentEntity() const
{
    return aDesc->ReferentEntity();
}
*/

const TypeDescriptor * 
STEPattribute::ReferentType() const
{
    return aDesc->ReferentType();
}

BOOLEAN 
STEPattribute::Nullable() const
{
    return (aDesc->Optionality().asInt() == T);
}
