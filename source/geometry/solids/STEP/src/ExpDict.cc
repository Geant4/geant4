

//



//
// $Id: ExpDict.cc,v 1.2 1999-05-21 20:20:47 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/ExpDict.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#include <ExpDict.h> 
#include <STEPaggregate.h> 

/*
const TypeDescriptor * const t_INTEGER_TYPE = & TypeDescriptor
                       ("INTEGER",     // Name
		       INTEGER_TYPE, // FundamentalType
		       "INTEGER");   // Description
extern const TypeDescriptor _t_INTEGER_TYPE
                       ("INTEGER",     // Name
		       INTEGER_TYPE, // FundamentalType
		       "INTEGER");   // Description
		       */
/*const TypeDescriptor * const t_INTEGER_TYPE = &_t_INTEGER_TYPE;*/

/*extern const TypeDescriptor _t_REAL_TYPE ("REAL", REAL_TYPE, "Real");*/
/*const TypeDescriptor * const t_REAL_TYPE = &_t_REAL_TYPE;*/

/*extern const TypeDescriptor _t_STRING_TYPE ("STRING", STRING_TYPE, "String");*/
/*const TypeDescriptor * const t_STRING_TYPE = &_t_STRING_TYPE;*/

/*extern const TypeDescriptor _t_BINARY_TYPE ("BINARY", BINARY_TYPE, "Binary") ;*/
/*const TypeDescriptor * const t_BINARY_TYPE = &_t_BINARY_TYPE;*/

/*extern const TypeDescriptor _t_BOOLEAN_TYPE ("BOOLEAN", BOOLEAN_TYPE, "Boolean") ;*/
/*const TypeDescriptor * const t_BOOLEAN_TYPE = &_t_BOOLEAN_TYPE;*/

/*extern const TypeDescriptor _t_LOGICAL_TYPE ("LOGICAL", LOGICAL_TYPE, "Logical") ;*/
/*const TypeDescriptor * const t_LOGICAL_TYPE = &_t_LOGICAL_TYPE;*/
 
/*extern const TypeDescriptor _t_NUMBER_TYPE ("NUMBER", NUMBER_TYPE, "Number") ;*/
/*const TypeDescriptor * const t_NUMBER_TYPE = &_t_NUMBER_TYPE;*/

/*extern const TypeDescriptor _t_GENERIC_TYPE ("GENERIC", GENERIC_TYPE, "Generic") ;*/
/*const TypeDescriptor * const t_GENERIC_TYPE = &_t_GENERIC_TYPE;*/



EnumAggregate * create_EnumAggregate()
{
    return new EnumAggregate; 
}

GenericAggregate * create_GenericAggregate()
{ 
    return new GenericAggregate; 
}

EntityAggregate * create_EntityAggregate()
{
    return new EntityAggregate; 
}

SelectAggregate * create_SelectAggregate()
{
    return new SelectAggregate; 
}

StringAggregate * create_StringAggregate()
{
    return new StringAggregate; 
}

BinaryAggregate * create_BinaryAggregate()
{
    return new BinaryAggregate; 
}

RealAggregate * create_RealAggregate()
{
    return new RealAggregate; 
}

IntAggregate * create_IntAggregate()
{
    return new IntAggregate; 
}

const EntityDescriptor * 
EntityDescItr::NextEntityDesc()
{
    if(cur)
    {
	const EntityDescriptor *ed = cur->EntityDesc();
	cur = (EntityDescLinkNode *)( cur->NextNode() );
	return ed;
    }
    return 0;
} 

const AttrDescriptor * 
AttrDescItr::NextAttrDesc()
{
    if(cur)
    {
	const AttrDescriptor *ad = cur->AttrDesc();
	cur = (AttrDescLinkNode *)( cur->NextNode() );
	return ad;
    }
    return 0;
} 

const InverseAttrDescriptor * 
InverseADItr::NextInverseAttrDesc()
{
    if(cur)
    {
	const InverseAttrDescriptor *iad = cur->InverseAttrDesc();
	cur = (InverseAttrDescLinkNode *)( cur->NextNode() );
	return iad;
    }
    return 0;
} 

const TypeDescriptor * 
TypeDescItr::NextTypeDesc()
{
    if(cur)
    {
	const TypeDescriptor *td = cur->TypeDesc();
	cur = (TypeDescLinkNode *)( cur->NextNode() );
	return td;
    }
    return 0;
} 

///////////////////////////////////////////////////////////////////////////////
// AttrDescriptor functions
///////////////////////////////////////////////////////////////////////////////

const char *
AttrDescriptor::AttrExprDefStr(SCLstring & s) const
{
  s = Name ();
  s.Append (" : ");
  if(_optional.asInt() == T)
    s.Append( "OPTIONAL ");
  if(DomainType())
	s.Append (DomainType()->AttrTypeName());
  return s.chars();
}    

const BASE_TYPE 
AttrDescriptor::BaseType() const
{
    if(_domainType)
	return _domainType->BaseType();
    return UNKNOWN_TYPE;
}

int 
AttrDescriptor::IsAggrType() const
{
    return ReferentType()->IsAggrType();
}

const BASE_TYPE 
AttrDescriptor::AggrElemType() const
{
    if(IsAggrType())
    {
	return ReferentType()->AggrElemType();
    }
    return UNKNOWN_TYPE;
}

const TypeDescriptor *
AttrDescriptor::AggrElemTypeDescriptor() const
{
    if(IsAggrType())
    {
	return ReferentType()->AggrElemTypeDescriptor();
    }
    return 0;
}

const TypeDescriptor * 
AttrDescriptor::NonRefTypeDescriptor() const
{
    if(_domainType)
	return _domainType->NonRefTypeDescriptor();
    return 0;
}

const BASE_TYPE 
AttrDescriptor::NonRefType() const
{
    if(_domainType)
	return _domainType->NonRefType();
    return UNKNOWN_TYPE;
}

const BASE_TYPE 
AttrDescriptor::Type() const
{
    if(_domainType)
	return _domainType->Type();
    return UNKNOWN_TYPE;
}

	// right side of attr def
// NOTE this returns a \'const char * \' instead of an SCLstring
const char *
AttrDescriptor::TypeName() const
{
    if(_domainType)
	return _domainType->AttrTypeName();
    else
	return "";
}

	// an expanded right side of attr def
const char *
AttrDescriptor::ExpandedTypeName(SCLstring & s) const
{
    s.set_null();
    if ((LOGICAL) Derived () == sdaiTRUE) s = "DERIVE  ";
    if(_domainType)
    {
	SCLstring tmp;
	return s.Append (_domainType->TypeString(tmp));
    }
    else
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// InverseAttrDescriptor functions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// EnumDescriptor functions
///////////////////////////////////////////////////////////////////////////////

EnumTypeDescriptor::EnumTypeDescriptor (const char * nm, BASE_TYPE ft, 
					const char * d, EnumCreator f)
: TypeDescriptor (nm, ft, d), CreateNewEnum(f)
{
}

///////////////////////////////////////////////////////////////////////////////
// EntityDescriptor functions
///////////////////////////////////////////////////////////////////////////////

EntityDescriptor::EntityDescriptor ( ) : _abstractEntity("U")
{
    _originatingSchema = 0;
//    _derivedAttr = new StringAggregate;
/*
    _subtypes = 0;
    _supertypes = 0;
    _explicitAttr = 0;
    _inverseAttr = 0;
*/
}

EntityDescriptor::EntityDescriptor (const char * name, // i.e. char *
				    SchemaDescriptor *origSchema, 
				    LOGICAL abstractEntity, // F U or T
				    Creator f
/*
				    EntityDescriptorList *subtypes,
				    EntityDescriptorList *supertypes,
				    AttrDescriptorList *explicitAttr,
				    StringAggregate *derivedAttr,
				    InverseAttrDescriptorList *inverseAttr
*/
				    ) 
: TypeDescriptor (name, ENTITY_TYPE, name), _originatingSchema (origSchema),  
  _abstractEntity(abstractEntity), NewSTEPentity(f)
{
/*
    _subtypes = subtypes;
    _supertypes = supertypes;
    _explicitAttr = explicitAttr;
    _derivedAttr = derivedAttr;
    _inverseAttr = inverseAttr;
*/
}

EntityDescriptor::~EntityDescriptor ()  
{
}

const TypeDescriptor * 
EntityDescriptor::IsA (const TypeDescriptor * td) const 
{ if (td -> NonRefType () == ENTITY_TYPE)
    return IsA ((EntityDescriptor *) td);
  else return 0;
}

const EntityDescriptor * 
EntityDescriptor::IsA (const EntityDescriptor * other)  const {
  const EntityDescriptor * found =0;
  const EntityDescLinkNode * link = (const EntityDescLinkNode *) (GetSupertypes ().GetHead());

  if (this == other) return other;
  else {
    while (link && ! found)  {
      found = link -> EntityDesc() -> IsA (other);
      link = (EntityDescLinkNode *) link -> NextNode ();
    }
  }
  return found;
}
/*
EntityDescriptor::FindLongestAttribute()
{
    AttrDescLinkNode *attrPtr = 
		  (AttrDescLinkNode *)(ed->ExplicitAttr().GetHead());
    while( attrPtr != 0)
    {
	if(attrPtr->AttrDesc()->IsEntityType())
	    maxAttrLen = max(maxAttrLen, 
		    (strlen(attrPtr->AttrDesc()->EntityType()->Name()) +
		      strlen(attrPtr->AttrDesc()->Name()) + 3
		     )
		    );
	else
	    maxAttrLen = max(maxAttrLen, 
	      (strlen(attrPtr->AttrDesc()->DomainType()->NameOrDescription()) +
	       strlen(attrPtr->AttrDesc()->Name()) + 3
	      )
		    );
	attrPtr = (AttrDescLinkNode *)attrPtr->NextNode();
    }
}
*/
///////////////////////////////////////////////////////////////////////////////
// TypeDescriptor functions
///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////
	// This is a fully expanded description of the type.
	// This returns a string like the _description member variable
	// except it is more thorough of a description where possible
	// e.g. if the description contains a TYPE name it will also
	// be explained.
	///////////////////////////////////////////////////////////////////////
const char *
TypeDescriptor::TypeString(SCLstring & s) const
{
    switch(Type())
    {
	  case REFERENCE_TYPE:
		if(Name())
		{
		  s.Append ( "TYPE ");
		  s.Append (Name());
		  s.Append ( " = ");
		}
		if(Description())
		  s.Append ( Description ());
		if(ReferentType())
		{
		  s.Append ( " -- ");
		  SCLstring tmp;
		  s.Append ( ReferentType()->TypeString(tmp));
		}
		return s;

	  case INTEGER_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("Integer");
		break;

	  case STRING_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("String");
		break;

	  case REAL_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("Real");
		break;

	  case ENUM_TYPE:
		s = "Enumeration: ";
		if(Name())
		{
		  s.Append ( "TYPE ");
		  s.Append (Name());
		  s.Append ( " = ");
		}
		if(Description())
		  s.Append ( Description ());
		break;

	  case BOOLEAN_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("Boolean: F, T");
		break;
	  case LOGICAL_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("Logical: F, T, U");
		break;
	  case NUMBER_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("Number");
		break;
	  case BINARY_TYPE:
		s.set_null();
		if(_referentType != 0)
		{
		    s = "TYPE ";
		    s.Append (Name());
		    s.Append ( " = ");
		}
		s.Append("Binary");
		break;
	  case ENTITY_TYPE:
		s = "Entity: ";
		if(Name())
		  s.Append (Name());
		break;
	  case AGGREGATE_TYPE:
	  case ARRAY_TYPE:		// DAS
	  case BAG_TYPE:		// DAS
	  case SET_TYPE:		// DAS
	  case LIST_TYPE:		// DAS
		s = Description();
		if(ReferentType())
		{
		  s.Append ( " -- ");
		  SCLstring tmp;
		  s.Append ( ReferentType()->TypeString(tmp));
		}
		break;
	  case SELECT_TYPE:
		  s.Append ( Description ());
		break;
	  case GENERIC_TYPE:
	  case UNKNOWN_TYPE:
		s = "Unknown";
		break;
    } // end switch
  return s;

}
/* this works
    if( ( (ReferentType() != 0) || (ReferentEntity() != 0) ) && Name())
    {
	strcat(tStr, "TYPE ");
	strcat(tStr, Name());
	strcat(tStr, " = ");
    }
    if(Description())
	strcat(tStr, Description());
    if(ReferentType())
    {
	strcat(tStr, " -- ");
	SCLstring tmp;
	strcat(tStr, ReferentType()->TypeString(tmp));
    }
    else if(ReferentEntity())
    {
	strcat(tStr, " -- ");
	strcat(tStr, "Entity: ");
	strcat(tStr, ReferentEntity()->Name());
    }
    return tStr;
}
*/


const TypeDescriptor * 
TypeDescriptor::IsA (const TypeDescriptor * other)  const {
  if (this == other)  return other;
  return 0;
}

const TypeDescriptor * 
TypeDescriptor::IsA (const char * other) const  {
  if (!Name())  return 0;
  if (!strcmp (Name (), PrettyTmpName (other)))  // this is the type
    return this;
  return (ReferentType () ? ReferentType () -> IsA (other) : 0);
}

	///////////////////////////////////////////////////////////////////////
	// the first BASE_TYPE that is not REFERENCE_TYPE (the first 
	// TypeDescriptor *_referentType that does not have REFERENCE_TYPE 
	// for it's fundamentalType variable).  This would return the same 
	// as BaseType() for fundamental types.  An aggregate type
	// would return AGGREGATE_TYPE then you could find out the type of
	// an element by calling AggrElemType().  Select types
	// would work the same?
	///////////////////////////////////////////////////////////////////////

const BASE_TYPE 
TypeDescriptor::NonRefType() const
{
    const TypeDescriptor *td = NonRefTypeDescriptor();
    if(td)
	return td->FundamentalType();
    return UNKNOWN_TYPE;
}


const TypeDescriptor *
TypeDescriptor::NonRefTypeDescriptor() const
{  
    const TypeDescriptor *td = this;

    while ( td->ReferentType() ) {
      if (td->Type() != REFERENCE_TYPE) 
	  return td;
      td = td->ReferentType();
    }

    return td;
}

	///////////////////////////////////////////////////////////////////////
	// This returns the BASE_TYPE of the first non-aggregate element of 
	// an aggregate
	///////////////////////////////////////////////////////////////////////

TypeDescriptor::IsAggrType() const
{
    switch(NonRefType())
    {
      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
	return 1;

      default:
	return 0;
    }
}

const BASE_TYPE 
TypeDescriptor::AggrElemType() const
{
    const TypeDescriptor *aggrElemTD = AggrElemTypeDescriptor();
    if(aggrElemTD)
    {
	return aggrElemTD->Type();
    }
    return UNKNOWN_TYPE;
}

const TypeDescriptor *
TypeDescriptor::AggrElemTypeDescriptor() const
{
    const TypeDescriptor *aggrTD = NonRefTypeDescriptor();
    const TypeDescriptor *aggrElemTD = aggrTD->ReferentType();
    if(aggrElemTD)
    {
	aggrElemTD = aggrElemTD->NonRefTypeDescriptor();
    }
    return aggrElemTD;
}

	////////////////////////////////////////////////////////////
	// This is the underlying type of this type. For instance:
	// TYPE count = INTEGER;
	// TYPE ref_count = count;
	// TYPE count_set = SET OF ref_count;
	//  each of the above will Generate a TypeDescriptor and for 
	//  each one, BASE_TYPE BaseType() will return INTEGER_TYPE
	//  TypeDescriptor *BaseTypeDescriptor() returns the TypeDescriptor 
	//  for Integer
	////////////////////////////////////////////////////////////

const BASE_TYPE
TypeDescriptor::BaseType() const
{
    const TypeDescriptor *td = BaseTypeDescriptor();
    if(td)
	return td->FundamentalType();
    else
	return ENTITY_TYPE;
}

const TypeDescriptor *
TypeDescriptor::BaseTypeDescriptor() const
{
    const TypeDescriptor *td = this;

    while (td -> ReferentType ()) td = td->ReferentType ();
    return td;
}
#ifdef NOT_YET
///////////////////////////////////////////////////////////////////////////////
// EnumerationTypeDescriptor functions
///////////////////////////////////////////////////////////////////////////////
EnumerationTypeDescriptor::EnumerationTypeDescriptor( ) 
{
    _elements = new StringAggregate; 
}
#endif
///////////////////////////////////////////////////////////////////////////////
// SelectTypeDescriptor functions
///////////////////////////////////////////////////////////////////////////////
const TypeDescriptor *
SelectTypeDescriptor::IsA (const TypeDescriptor * other) const
{  return TypeDescriptor::IsA (other);  }

const TypeDescriptor *
SelectTypeDescriptor::CanBe (const TypeDescriptor * other) const
{
  const TypeDescriptor * found =0;
  TypeDescItr elements (GetElements()) ;
  const TypeDescriptor * td =0;

  if (this == other) return other;
  while (td = elements.NextTypeDesc ())  {
    if (found = (td -> CanBe (other))) return found;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// AggrTypeDescriptor functions
///////////////////////////////////////////////////////////////////////////////

AggrTypeDescriptor::AggrTypeDescriptor( ) : 
		_uniqueElements("UNKNOWN_TYPE")
{
    _bound1 = -1;
    _bound2 = -1;
    _aggrDomainType = 0;
}

AggrTypeDescriptor::AggrTypeDescriptor(SdaiInteger b1, 
						 SdaiInteger b2, 
						 LOGICAL uniqElem, 
						 TypeDescriptor *aggrDomType)
	      : _bound1(b1), _bound2(b2), _uniqueElements(uniqElem)
{
    _aggrDomainType = aggrDomType;
}

AggrTypeDescriptor::~AggrTypeDescriptor()
{
}
