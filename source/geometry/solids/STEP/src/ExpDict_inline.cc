// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExpDict_inline.cc,v 1.1 1999-01-07 16:08:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/ExpDict.inline.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#include <ExpDict.h>

///////////////////////////////////////////////////////////////////////////////
// TypeDescriptor functions
///////////////////////////////////////////////////////////////////////////////


TypeDescriptor::TypeDescriptor( ) 
: _name (0), _fundamentalType(UNKNOWN_TYPE), _referentType (0), _description (0)
{
}

TypeDescriptor::TypeDescriptor
(const char * nm, BASE_TYPE ft, const char * d )
:  _name (nm), _fundamentalType (ft), _referentType (0), _description (d)  
{
}

const char * 
TypeDescriptor::BaseTypeName ()  const
{
  return BaseTypeDescriptor () ?  BaseTypeDescriptor () -> Name () : 0;
}

const TypeDescriptor * 
TypeDescriptor::BaseTypeIsA (const TypeDescriptor * td) const
{
  switch (NonRefType ()) {
  case AGGREGATE_TYPE:
    return AggrElemTypeDescriptor () -> IsA (td);
  case ENTITY_TYPE:
  case SELECT_TYPE: 
  default:
    return IsA (td);
  }

}

///////////////////////////////////////////////////////////////////////////////
// AttrDescriptor functions
///////////////////////////////////////////////////////////////////////////////

AttrDescriptor::AttrDescriptor(
		       const char * name,		// i.e. char *
		       const TypeDescriptor *domainType, 
		       LOGICAL optional,	// i.e. F U or T
		       LOGICAL unique,	// i.e. F U or T
		       LOGICAL derived,	// i.e. F U or T
		       const EntityDescriptor & owner )
: _name (name), _domainType (domainType), _optional(optional), _unique(unique), _derived (derived), 

#ifdef __O3DB__
_owner (&owner)
#else
_owner ((EntityDescriptor&)owner)
#endif

{
}

AttrDescriptor::~AttrDescriptor () { }
