

//



//
// $Id: ExpDict.h,v 1.2 1999-05-21 20:20:29 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef EXPDICT_H
#define EXPDICT_H

/*
* NIST STEP Core Class Library
* clstepcore/ExpDict.h
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 


#ifdef __O3DB__
#include <OpenOODB.h>
#endif

class STEPentity;
typedef  STEPentity * (* Creator) () ;
//class StringAggregate;

#include <SingleLinkList.h>
#include <sdai.h> 

#include <baseType.h>

class SchemaDescriptor;
class AttrDescriptor;
class InverseAttrDescriptor;
class EntityDescriptor;
class TypeDescriptor;
class EnumerationTypeDescriptor;
class AggrTypeDescriptor;
class ArrayTypeDescriptor;
class SetTypeDescriptor;
class ListTypeDescriptor;
class SelectTypeDescriptor;
class StringTypeDescriptor;
class BagTypeDescriptor;
class RealTypeDescriptor;
class EntityDescLinkNode;
class EntityDescriptorList;
class AttrDescLinkNode;
class AttrDescriptorList;
class InverseAttrDescLinkNode;
class InverseAttrDescriptorList;
class TypeDescLinkNode;
class TypeDescriptorList;

/*
**  I tried these variations on the TypeDescriptor to get them to be
*   initialized globally.  I couldn\'t do it.  They are now initialized
*   in the Registry constructor (in Registry.inline.cc

extern const TypeDescriptor t_INTEGER_TYPE;
extern const TypeDescriptor t_REAL_TYPE;
extern const TypeDescriptor t_NUMBER_TYPE;
extern const TypeDescriptor t_STRING_TYPE;
extern const TypeDescriptor t_BINARY_TYPE;
extern const TypeDescriptor t_BOOLEAN_TYPE;
extern const TypeDescriptor t_LOGICAL_TYPE;

#define t_INTEGER_TYPE &_t_INTEGER_TYPE
#define t_REAL_TYPE  &_t_REAL_TYPE
#define t_NUMBER_TYPE &_t_NUMBER_TYPE
#define t_STRING_TYPE &_t_STRING_TYPE
#define t_BINARY_TYPE &_t_BINARY_TYPE
#define t_BOOLEAN_TYPE &_t_BOOLEAN_TYPE
#define t_LOGICAL_TYPE &_t_LOGICAL_TYPE

extern const TypeDescriptor * const t_INTEGER_TYPE;
extern const TypeDescriptor * const t_REAL_TYPE;
extern const TypeDescriptor * const t_NUMBER_TYPE;
extern const TypeDescriptor * const t_STRING_TYPE;
extern const TypeDescriptor * const t_BINARY_TYPE;
extern const TypeDescriptor * const t_BOOLEAN_TYPE;
extern const TypeDescriptor * const t_LOGICAL_TYPE;
*/

extern const TypeDescriptor *  t_INTEGER_TYPE;
extern const TypeDescriptor *  t_REAL_TYPE;
extern const TypeDescriptor *  t_NUMBER_TYPE;
extern const TypeDescriptor *  t_STRING_TYPE;
extern const TypeDescriptor *  t_BINARY_TYPE;
extern const TypeDescriptor *  t_BOOLEAN_TYPE;
extern const TypeDescriptor *  t_LOGICAL_TYPE;

///////////////////////////////////////////////////////////////////////////////
// SchemaDescriptor - a class of this type is generated and contains
// the name of the schema.
///////////////////////////////////////////////////////////////////////////////

class SchemaDescriptor { 

  protected:
	const char *  _name ;
  public:  

	SchemaDescriptor (const char *schemaName ) { _name = schemaName; }
	virtual ~SchemaDescriptor () { }

	const char * Name() const	{ return _name; }
	void Name (const char *  n)  { _name = n; }
};


///////////////////////////////////////////////////////////////////////////////
// EntityDescriptor
// An instance of this class will be generated for each entity type
// found in the schema.  This should probably be derived from the
// CreatorEntry class (see STEPentity.h).  Then the binary tree that the
// current software  builds up containing the entities in the schema
// will be building the same thing but using the new schema info.
// nodes (i.e. EntityDesc nodes) for each entity.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

class EntityDescLinkNode : public SingleLinkNode {

  private:
  protected:
    EntityDescriptor * _entityDesc;

  public:
    EntityDescLinkNode() { _entityDesc = 0; }
    virtual ~EntityDescLinkNode() { }

    EntityDescriptor *EntityDesc() const { return _entityDesc; }
    void EntityDesc(EntityDescriptor *ed) { _entityDesc = ed; }
};

class EntityDescriptorList : public SingleLinkList {

  private:
  protected:
  public:
    EntityDescriptorList()  { }
    virtual ~EntityDescriptorList() { }

    virtual SingleLinkNode * NewNode () { return new EntityDescLinkNode; }

    EntityDescLinkNode * AddNode (EntityDescriptor * ed) { 
	EntityDescLinkNode *node = (EntityDescLinkNode *) NewNode();
	node->EntityDesc(ed);
	SingleLinkList::AppendNode(node);
	return node;
    }
};

class EntityDescItr
{
  protected:
    const EntityDescriptorList &edl;
    const EntityDescLinkNode *cur;

  public:
    EntityDescItr(const EntityDescriptorList &edList) : edl(edList)
    {
	cur = (EntityDescLinkNode *)( edl.GetHead() );
    }
    ~EntityDescItr() { };

    void ResetItr() { cur = (EntityDescLinkNode *)( edl.GetHead() ); }

    const EntityDescriptor * NextEntityDesc();
};

///////////////////////////////////////////////////////////////////////////////

class AttrDescLinkNode : public SingleLinkNode {
  private:
  protected:
    AttrDescriptor *_attrDesc;
  public:
    AttrDescLinkNode()  { _attrDesc = 0; }
    virtual ~AttrDescLinkNode() { }

    const class AttrDescriptor *AttrDesc() const { return _attrDesc; }
    void AttrDesc(AttrDescriptor *ad) { _attrDesc = ad; }
};

class AttrDescriptorList : public SingleLinkList {
  private:
  protected:
  public:
    AttrDescriptorList()  { }
    virtual ~AttrDescriptorList() { }

    virtual SingleLinkNode * NewNode () { return new AttrDescLinkNode; }
    AttrDescLinkNode * AddNode (AttrDescriptor * ad) {
	AttrDescLinkNode *node = (AttrDescLinkNode *) NewNode();
	node->AttrDesc(ad);
	SingleLinkList::AppendNode(node);
	return node;
    }    
};

class AttrDescItr
{
  protected:
    const AttrDescriptorList &adl;
    const AttrDescLinkNode *cur;

  public:
    AttrDescItr(const AttrDescriptorList &adList) : adl(adList)
    {
	cur = (AttrDescLinkNode *)( adl.GetHead() );
    }
    ~AttrDescItr() { };

    void ResetItr() { cur = (AttrDescLinkNode *)( adl.GetHead() ); }

    const AttrDescriptor * NextAttrDesc();
};

///////////////////////////////////////////////////////////////////////////////

class InverseAttrDescLinkNode : public  AttrDescLinkNode {
  private:
  protected:
    class InverseAttrDescriptor *_invAttrDesc;
  public:
    InverseAttrDescLinkNode() { _invAttrDesc = 0; }
    virtual ~InverseAttrDescLinkNode() { }

    const InverseAttrDescriptor *InverseAttrDesc() const { return _invAttrDesc; }
    void InverseAttrDesc(InverseAttrDescriptor *iad) { _invAttrDesc = iad; }
};

class InverseAttrDescriptorList : public  AttrDescriptorList {
  private:
  protected:
  public:
    InverseAttrDescriptorList() { }
    virtual ~InverseAttrDescriptorList() { }
    virtual SingleLinkNode * NewNode () { return new InverseAttrDescLinkNode; }
};

class InverseADItr
{
  protected:
    const InverseAttrDescriptorList &iadl;
    const InverseAttrDescLinkNode *cur;

  public:
    InverseADItr (const InverseAttrDescriptorList &iadList) : iadl(iadList)
    {
	cur = (InverseAttrDescLinkNode *)( iadl.GetHead() );
    }
    ~InverseADItr() { };

    void ResetItr() 
	{ cur = (InverseAttrDescLinkNode *)( iadl.GetHead() ); }

    const InverseAttrDescriptor * NextInverseAttrDesc();
};

///////////////////////////////////////////////////////////////////////////////

class TypeDescLinkNode : public SingleLinkNode {
  private:
  protected:
    TypeDescriptor *_typeDesc;
  public:
    TypeDescLinkNode() { _typeDesc = 0; }
    virtual ~TypeDescLinkNode() { }

    const TypeDescriptor *TypeDesc() const { return _typeDesc; }
    void TypeDesc(TypeDescriptor *td) { _typeDesc = td; }
};

class TypeDescriptorList : public SingleLinkList {
  private:
  protected:
  public:
    TypeDescriptorList() { }
    virtual ~TypeDescriptorList() { }

    virtual SingleLinkNode * NewNode () { return new TypeDescLinkNode; }

    TypeDescLinkNode * AddNode (TypeDescriptor * td) {
	TypeDescLinkNode *node = (TypeDescLinkNode *) NewNode();
	node->TypeDesc(td);
	SingleLinkList::AppendNode(node);
	return node;
    }    
};

class TypeDescItr
{
  protected:
    const TypeDescriptorList &tdl;
    const TypeDescLinkNode *cur;

  public:
    TypeDescItr (const TypeDescriptorList &tdList) : tdl(tdList)
    {
	cur = (TypeDescLinkNode *)( tdl.GetHead() );
    }
    ~TypeDescItr() { };

    void ResetItr() { cur = (TypeDescLinkNode *)( tdl.GetHead() ); }

    const TypeDescriptor * NextTypeDesc();
};

///////////////////////////////////////////////////////////////////////////////
// AttrDescriptor
// An instance of this class will be generated for each attribute for
// an Entity.  They will be pointed to by the EntityTypeDescriptors.
///////////////////////////////////////////////////////////////////////////////

class AttrDescriptor { 

  protected:
	const char *  _name ;	// the attributes name
		// this defines the domain of the attribute
	const TypeDescriptor * _domainType ;
	SdaiLogical _optional;
	SdaiLogical _unique;
	SdaiLogical _derived;
	
#ifdef __O3DB__
	const EntityDescriptor * _owner ;  // the owning entityDescriptor
#else
	const EntityDescriptor & _owner ;  // the owning entityDescriptor
#endif
  public:  

	AttrDescriptor(
		       const char * name,		// i.e. char *
		       const TypeDescriptor *domainType, 
		       LOGICAL optional,	// i.e. F U or T
		       LOGICAL unique,	// i.e. F U or T
		       LOGICAL derived,	// i.e. F U or T
		       const EntityDescriptor & owner 
		   );
	virtual ~AttrDescriptor ();

		// the attribute Express def
	const char *AttrExprDefStr(SCLstring & s) const;

		// left side of attr def
	const char * Name()	const	   { return _name; }
	void     Name (const char *  n) { _name = n; }

	   // BaseType() is the underlying type of this attribute. 
	   // NonRefType() is the first non REFERENCE_TYPE type
	   // e.g. Given attributes of each of the following types 
	   // TYPE count = INTEGER;
	   // TYPE ref_count = count;
	   // TYPE count_set = SET OF ref_count;
	   //  BaseType() will return INTEGER_TYPE for an attr of each type.
	   //  BaseTypeDescriptor() returns the TypeDescriptor for Integer
	   //  NonRefType() will return INTEGER_TYPE for the first two. For an
	   //    attribute of type count_set NonRefType() would return 
	   //    AGGREGATE_TYPE
	   //  NonRefTypeDescriptor() returns the TypeDescriptor for Integer
	   //     for the first two and a TypeDescriptor for an
	   //     aggregate for the last.

	const PrimitiveType BaseType() const;
	const TypeDescriptor *BaseTypeDescriptor() const;

	   // the first PrimitiveType that is not REFERENCE_TYPE (the first 
	   // TypeDescriptor *_referentType that does not have REFERENCE_TYPE 
	   // for it's fundamentalType variable).  This would return the same 
	   // as BaseType() for fundamental types.  An aggregate type
	   // would return AGGREGATE_TYPE then you could find out the type of
	   // an element by calling AggrElemType().  Select types
	   // would work the same?

	const PrimitiveType NonRefType() const;
	const TypeDescriptor *NonRefTypeDescriptor() const;

	int   IsAggrType() const;
	const PrimitiveType	AggrElemType() const;
	const TypeDescriptor *AggrElemTypeDescriptor() const;

		// The type of the attributes TypeDescriptor
	const PrimitiveType Type() const;
	const char * TypeName() const;	// right side of attr def

			// an expanded right side of attr def
	const char *ExpandedTypeName(SCLstring & s) const;

	int RefersToType() const	{ return !(_domainType == 0); }

	const TypeDescriptor * ReferentType() const { return _domainType; }
	const TypeDescriptor * DomainType() const { return _domainType; }
	void DomainType (const TypeDescriptor *td)  { _domainType = td; }
	void ReferentType(const TypeDescriptor *td) { _domainType = td; }

	const SdaiLogical & Optional() const { return _optional; }
	void Optional (SdaiLogical &opt)	{ _optional.put(opt.asInt()); }
	void Optional (LOGICAL opt)	{ _optional.put(opt); }
	void Optional (const char *opt) { _optional.put(opt); }

	const SdaiLogical & Unique() const { return _unique; }
	void Unique (SdaiLogical uniq)	{ _unique.put(uniq.asInt()); }
	void Unique (LOGICAL uniq)	{ _unique.put(uniq); }
	void Unique (const char *uniq)	{ _unique.put(uniq); }

	const SdaiLogical & Derived() const { return _derived; }
	void Derived (SdaiLogical x)	{ _derived.put(x.asInt()); }
	void Derived (LOGICAL x)	{ _derived.put(x); }
	void Derived (const char *x)	{ _derived.put(x); }

	const SdaiLogical & Optionality() const { return _optional; }
	void Optionality (SdaiLogical &opt)  { _optional.put(opt.asInt()); }
	void Optionality (LOGICAL opt)	   { _optional.put(opt); }
	void Optionality (const char *opt) { _optional.put(opt); }

	const SdaiLogical & Uniqueness() const	{ return _unique; }
	void Uniqueness (SdaiLogical uniq)	{ _unique.put(uniq.asInt()); }
	void Uniqueness (LOGICAL uniq)		{ _unique.put(uniq); }
	void Uniqueness (const char *uniq)	{ _unique.put(uniq); }

#ifdef __O3DB__
	const EntityDescriptor  & Owner() const	{ return *_owner; }
#else
	const EntityDescriptor  & Owner() const	{ return _owner; }
#endif
};


///////////////////////////////////////////////////////////////////////////////
// InverseAttrDescriptor
///////////////////////////////////////////////////////////////////////////////

class InverseAttrDescriptor  :    public AttrDescriptor  { 

  protected:
	AttrDescriptor * _inverseAttr ;
  public:  

	InverseAttrDescriptor(
		       const char * name,		// i.e. char *
		       TypeDescriptor *domainType, 
		       LOGICAL optional,	// i.e. F U or T*/
		       LOGICAL unique,	// i.e. F U or T
//		       LOGICAL derived,	// derived will always be F 
		       const EntityDescriptor & owner, 
		       AttrDescriptor *inverseAttr =0
		   ) : AttrDescriptor( name, domainType, optional, unique, 
				       F, owner ), _inverseAttr (inverseAttr)
		{ }
	virtual ~InverseAttrDescriptor () { }
	
	class AttrDescriptor * InverseAttribute() 
                { return _inverseAttr; } 
        void InverseOf (AttrDescriptor * invAttr) 
                { _inverseAttr = invAttr; } 
};

///////////////////////////////////////////////////////////////////////////////
// TypeDescriptor
// This class and the classes inherited from this class are used to describe   
// all types (base types and created types).  There will be an instance of this
// class generated for each type found in the schema.  
// A TypeDescriptor will be generated in three contexts:
// 1) to describe a base type - e.g. INTEGER, REAL, STRING.  There is only one 
//	TypeDescriptor created for each Express base type. Each of these will 
//	be pointed to by several other AttrDescriptors and TypeDescriptors)
// 2) to describe a type created by an Express TYPE statement.
//	e.g. TYPE label = STRING END_TYPE;
//	These TypeDescriptors will be pointed to by other AttrDescriptors (and 
//	TypeDescriptors) representing attributes (and Express TYPEs) that are 
//	of the type created by this Express TYPE.
// 3) to describe a type created in an attribute definition
//	e.g. part_label_grouping : ARRAY [1.10] label;
//	or part_codes : ARRAY [1.10] INTEGER;
//	In this #3 context there will not be a name associated with the type.
//	The TypeDescriptor created in this case will only be pointed to by the
//	single AttrDescriptor associated with the attribute it was created for.
///////////////////////////////////////////////////////////////////////////////

///// _name is the name of the type. 
	// In the case of the TypeDescriptors representing the Express base
	// types this will be the name of the base type.
	// In the case where this TypeDescriptor is representing an Express
	// TYPE it is the LEFT side of an Express TYPE statement (i.e. label 
	// as in TYPE label = STRING END_TYPE;) This name would in turn be 
	// found on the RIGHT side of an Express attribute definition (e.g. 
	// attr defined as part_label : label; )
	// In the case where this TypeDescriptor was generated to describe a
	// type created in an attr definition, it will be a null pointer (e.g
	// attr defined as part_label_grouping : ARRAY [1..10] label)
///// _fundamentalType is the 'type' of the type being represented by 
	//  the TypeDescriptor . i.e. the following 2 stmts
	//  would cause 2 TypeDescriptors to be generated - the 1st having
	//  _fundamentalType set to STRING_TYPE and for the 2nd to 
	//  REFERENCE_TYPE.
	// TYPE label = STRING END_TYPE; 
	// TYPE part_label = label END_TYPE;
	// part_label and label would be the value of the respective
	// _name member variables for the 2 TypeDescriptors.  
///// _referentType will point at another TypeDescriptor furthur specifying 
	//  the type in all cases except when the type is directly 
	//  an enum or select.  i.e. in the following... _referentType for 
	//  the 1st type does not point at anything and for the 2nd it does:
	// TYPE color = ENUMERATION OF (red, blue); END_TYPE;
	// TYPE color_ref = color; END_TYPE;
////// _fundamentalType being REFERENCE_TYPE (as would be the case for 
	// part_label and color_ref above) means that the _referentType 
	// member variable points at a TypeDescriptor representing a type
	// that has been defined in an Express TYPE stmt.
	//  Otherwise _fundamental type reflects 
	//  the type directly as in the type label above.  type label above
	//  has a _referentType that points at a TypeDescriptor for STRING 
	//  described in the next sentence (also see #1 above).
	// A TypeDescriptor would be generated for each of the EXPRESS base 
	// types (int, string, real, etc) having _fundamentalType member 
	// variables set to match the EXPRESS base type being represented.
//////_referentType
	// For the TypeDescriptors describing the EXPRESS base types this will
	// be a null pointer.  For all other TypeDescriptors this will point 
	// to another TypeDescriptor which furthur describes the type. e.g.
	// TYPE part_label = label END_TYPE; TYPE label = STRING END_TYPE; 
	// part_label's _referentType will point to the TypeDescriptor for
	// label.  label's _referentType will point to the TypeDescriptor
	// for STRING. The _fundamentalType for part_label will be
	// REFERENCE_TYPE and for label will be STRING_TYPE.  
	// The _fundamentalType for the EXPRESS base type STRING's
	// TypeDescriptor will be STRING_TYPE.
	// The _referentType member variable will in most cases point to
	// a subtype of TypeDescriptor.
//////_description
	// This is the string description of the type as found in the
	// EXPRESS file. e.g. aggr of [aggr of ...] [List of ...] someType 
	// It is the RIGHT side of an Express TYPE statement 
	// (i.e. LIST OF STRING as in 
	// TYPE label_group = LIST OF STRING END_TYPE;) 
	// It is the same as _name for EXPRESS base types TypeDescriptors (with
	// the possible exception of upper or lower case differences).

class TypeDescriptor { 

  protected:

		// the name of the type (see above)
	const char *  _name ;

		// the type of the type (see above). 
		// it is an enum see file clstepcore/baseType.h
//	BASE_TYPE _fundamentalType ;
	PrimitiveType _fundamentalType ;

		// furthur describes the type (see above)
		// most often (or always) points at a subtype.
	const TypeDescriptor * _referentType ;

		// Express file description (see above)
		// e.g. the right side of an Express TYPE stmt
	const char *  _description ;

  public:  

	TypeDescriptor (const char * nm, PrimitiveType ft, const char * d ); 
	TypeDescriptor ( ); 
	virtual ~TypeDescriptor () { }


		// the name of this type
	const char * Name() const   { return _name; }

		// The name that would be found on the right side of an 
		// attribute definition. In the case of a type defined like
		// TYPE name = STRING END_TYPE; 
		// with attribute definition   employee_name : name;
		// it would be the _name member variable. If it was a type
		// defined in an attribute it will be the _description
		// member variable since _name will be null. e.g. attr. def.
		// project_names : ARRAY [1..10] name;
	const char * AttrTypeName() const {
	    return _name ? _name : _description;
	}	    

		// This is a fully expanded description of the type.
		// This returns a string like the _description member variable
		// except it is more thorough of a description where possible
		// e.g. if the description contains a TYPE name it will also
		// be explained.
	const char *TypeString(SCLstring & s) const;

		// This TypeDescriptor's type
	const PrimitiveType Type() const	{ return _fundamentalType; }
	void  Type(const PrimitiveType type) 	{ _fundamentalType = type; }

	   // This is the underlying Express base type of this type. It will 
	   // be the type of the last TypeDescriptor following the 
	   // _referentType member variable pointers. e.g.
	   // TYPE count = INTEGER;
	   // TYPE ref_count = count;
	   // TYPE count_set = SET OF ref_count;
	   //  each of the above will Generate a TypeDescriptor and for 
	   //  each one, PrimitiveType BaseType() will return INTEGER_TYPE.
	   //  TypeDescriptor *BaseTypeDescriptor() returns the TypeDescriptor 
	   //  for Integer.
	const PrimitiveType       BaseType() const;
	const TypeDescriptor *BaseTypeDescriptor() const;
	const char * BaseTypeName () const;

	   // the first PrimitiveType that is not REFERENCE_TYPE (the first 
	   // TypeDescriptor *_referentType that does not have REFERENCE_TYPE 
	   // for it's fundamentalType variable).  This would return the same 
	   // as BaseType() for fundamental types.  An aggregate type
	   // would return AGGREGATE_TYPE then you could find out the type of
	   // an element by calling AggrElemType().  Select types
	   // would work the same?

	const PrimitiveType	NonRefType() const;
	const TypeDescriptor *NonRefTypeDescriptor() const;

	int   IsAggrType() const;
	const PrimitiveType	AggrElemType() const;
	const TypeDescriptor *AggrElemTypeDescriptor() const;

	const PrimitiveType FundamentalType() const { return _fundamentalType; }
	void FundamentalType (PrimitiveType ftype) { _fundamentalType = ftype; }

		// The TypeDescriptor for the type this type is based on 
	const TypeDescriptor * ReferentType() const { return _referentType; }
	void ReferentType (const TypeDescriptor * rtype) 
			{ _referentType = rtype; }

		// A description of this type's type. Basically you
		// get the right side of a TYPE statement minus END_TYPE.
		// For base type TypeDescriptors it is the same as _name.
	const char * Description() const	{ return _description; }
	void Description (const char * desc) { _description = desc; } 

        virtual const TypeDescriptor * IsA (const TypeDescriptor *) const;
        virtual const TypeDescriptor * BaseTypeIsA (const TypeDescriptor *) const;
        virtual const TypeDescriptor * IsA (const char *) const;
        virtual const TypeDescriptor * CanBe (const TypeDescriptor *n) const
                {  return TypeDescriptor::IsA (n);  }

};

typedef  STEPenumeration * (* EnumCreator) () ;

class EnumTypeDescriptor  :    public TypeDescriptor  { 
  public:
    EnumCreator CreateNewEnum;

    void AssignEnumCreator(EnumCreator f = 0)
    {
	CreateNewEnum = f;
    }

    STEPenumeration *CreateEnum()
    {
	if(CreateNewEnum)
	    return CreateNewEnum();
	else
	    return 0;
    }

    EnumTypeDescriptor ( ) { }
    EnumTypeDescriptor (const char * nm, PrimitiveType ft, const char * d, 
		    EnumCreator f =0 ); 

    virtual ~EnumTypeDescriptor () { }
};


class EntityDescriptor  :    public TypeDescriptor  { 

  protected:
	const SchemaDescriptor * _originatingSchema;
	SdaiLogical _abstractEntity;

	EntityDescriptorList _subtypes;   // OPTIONAL
	EntityDescriptorList _supertypes; // OPTIONAL

	AttrDescriptorList	  _explicitAttr; // OPTIONAL
//	StringAggregate		 * _derivedAttr;  // OPTIONAL  
	InverseAttrDescriptorList _inverseAttr;  // OPTIONAL

  public:
       // pointer to a function that will create a new instance of a STEPentity
        Creator NewSTEPentity;

	EntityDescriptor ( );
	EntityDescriptor (const char * name, // i.e. char *
			  SchemaDescriptor *origSchema, 
			  LOGICAL abstractEntity, // i.e. F U or T
			  Creator f =0		  
			  );

	virtual ~EntityDescriptor ();

	const SchemaDescriptor * OriginatingSchema()  const
			{ return _originatingSchema; }
	void OriginatingSchema (const SchemaDescriptor * os)
			{ _originatingSchema = os; }

	SdaiLogical & AbstractEntity()         { return _abstractEntity; } 
	void AbstractEntity (SdaiLogical &ae)
					   { _abstractEntity.put(ae.asInt()); }
	void AbstractEntity (LOGICAL ae)     { _abstractEntity.put(ae); }
	void AbstractEntity (const char *ae) { _abstractEntity.put(ae); }

	const EntityDescriptorList& Subtypes() const
	  { return _subtypes; } 

	const EntityDescriptorList& Supertypes() const
	  { return _supertypes; }

	const EntityDescriptorList& GetSupertypes()  const 
	  { return _supertypes; }

	const AttrDescriptorList& ExplicitAttr() const
	  { return _explicitAttr; }

//	StringAggregate  & DerivedAttr()	{ return *_derivedAttr; }

	const InverseAttrDescriptorList& InverseAttr() const { return _inverseAttr; }

        virtual const EntityDescriptor * IsA (const EntityDescriptor *) const;
        virtual const TypeDescriptor * IsA (const TypeDescriptor * td) const ;
        virtual const TypeDescriptor * IsA (const char * n) const  
                {  return TypeDescriptor::IsA (n);  }
        virtual const TypeDescriptor * CanBe (const TypeDescriptor *o) const
                {  return o -> IsA (this);  }

    // The following will be used by schema initialization functions

	void AddSubtype(EntityDescriptor *ed) 
		{ _subtypes.AddNode(ed); }

	void AddSupertype(EntityDescriptor *ed) 
		{ _supertypes.AddNode(ed); }

	void AddExplicitAttr(AttrDescriptor *ad)
		{ _explicitAttr.AddNode(ad); }

	void AddInverseAttr(InverseAttrDescriptor *ad)
		{ _inverseAttr.AddNode(ad); }
};


///////////////////////////////////////////////////////////////////////////////
// EnumerationTypeDescriptor
///////////////////////////////////////////////////////////////////////////////
#ifdef NOT_YET
class EnumerationTypeDescriptor  :    public TypeDescriptor  { 

  protected:
	StringAggregate  *_elements ;	  //  of  (null)

  public:  
	EnumerationTypeDescriptor ( );
	virtual ~EnumerationTypeDescriptor () { } 


	StringAggregate & Elements() { return *_elements; }
//	void Elements (StringAggregate  e);
};
#endif

class STEPaggregate;
class EnumAggregate;
class GenericAggregate;
class EntityAggregate;
class SelectAggregate;
class StringAggregate;
class BinaryAggregate;
class RealAggregate;
class IntAggregate;

typedef  STEPaggregate * (* AggregateCreator) () ;
typedef  EnumAggregate * (* EnumAggregateCreator) () ;
typedef  GenericAggregate * (* GenericAggregateCreator) () ;
typedef  EntityAggregate * (* EntityAggregateCreator) () ;
typedef  SelectAggregate * (* SelectAggregateCreator) () ;
typedef  StringAggregate * (* StringAggregateCreator) () ;
typedef  BinaryAggregate * (* BinaryAggregateCreator) () ;
typedef  RealAggregate * (* RealAggregateCreator) () ;
typedef  IntAggregate * (* IntAggregateCreator) () ;

EnumAggregate * create_EnumAggregate();

GenericAggregate * create_GenericAggregate();

EntityAggregate * create_EntityAggregate();

SelectAggregate * create_SelectAggregate();

StringAggregate * create_StringAggregate();

BinaryAggregate * create_BinaryAggregate();

RealAggregate * create_RealAggregate();

IntAggregate * create_IntAggregate();

///////////////////////////////////////////////////////////////////////////////
// AggrTypeDescriptor
// I think we decided on a simplistic representation of aggr. types for now?
// i.e. just have one AggrTypeDesc for Array of [List of] [set of] someType
// the inherited variable _referentType will point to the TypeDesc for someType
// So I don't believe this class was necessary.  If we were to retain
// info for each of the [aggr of]'s in the example above then there would be
// one of these for each [aggr of] above and they would be strung
// together by the _aggrDomainType variables.  If you can make this
// work then go for it.
///////////////////////////////////////////////////////////////////////////////

class AggrTypeDescriptor  :    public TypeDescriptor  { 

  protected:

    SdaiInteger  _bound1 ;
    SdaiInteger  _bound2 ;
    SdaiLogical _uniqueElements ;
    TypeDescriptor * _aggrDomainType ;
    AggregateCreator CreateNewAggr;

  public:  

    void AssignAggrCreator(AggregateCreator f = 0)
    {
	CreateNewAggr = f;
    }

    STEPaggregate *CreateAggregate()
    {
	if(CreateNewAggr)
	    return CreateNewAggr();
	else
	    return 0;
    }

    AggrTypeDescriptor ( ); 
    AggrTypeDescriptor(SdaiInteger b1, SdaiInteger b2, 
		       LOGICAL uniqElem, 
		       TypeDescriptor *aggrDomType);
    AggrTypeDescriptor (const char * nm, PrimitiveType ft, const char * d, 
			AggregateCreator f =0 )
	: TypeDescriptor (nm, ft, d), CreateNewAggr(f) { }
    virtual ~AggrTypeDescriptor ();


    SdaiInteger & Bound1() 		{ return _bound1; }
    void Bound1 (SdaiInteger  b1)    { _bound1 = b1; } 

    SdaiInteger & Bound2() 		{ return _bound2; }
    void Bound2 (SdaiInteger  b2)    { _bound2 = b2; } 

    SdaiLogical& UniqueElements()	{ return _uniqueElements; } 
    void UniqueElements (SdaiLogical &ue) 
					{ _uniqueElements.put(ue.asInt()); } 
    void UniquesElements (LOGICAL ue)     { _uniqueElements.put(ue); }
    void UniqueElements (const char *ue) { _uniqueElements.put(ue); }

    class TypeDescriptor * AggrDomainType()    { return _aggrDomainType; } 
    void AggrDomainType (TypeDescriptor * adt) { _aggrDomainType = adt; }
};

///////////////////////////////////////////////////////////////////////////////
// ArrayTypeDescriptor
///////////////////////////////////////////////////////////////////////////////

class ArrayTypeDescriptor  :    public AggrTypeDescriptor  { 

  protected:
	SdaiLogical _optionalElements ;
  public:  

    ArrayTypeDescriptor ( ) : _optionalElements("UNKNOWN_TYPE") { } 
    ArrayTypeDescriptor (LOGICAL optElem) : _optionalElements(optElem) { } 
    ArrayTypeDescriptor (const char * nm, PrimitiveType ft, const char * d, 
			AggregateCreator f =0 )
	: AggrTypeDescriptor (nm, ft, d, f), _optionalElements("UNKNOWN_TYPE") 
    { }

    virtual ~ArrayTypeDescriptor () {}


    SdaiLogical& OptionalElements()       { return _optionalElements; } 
    void OptionalElements (SdaiLogical &oe) 
				     { _optionalElements.put(oe.asInt()); } 
    void OptionalElements (LOGICAL oe)     { _optionalElements.put(oe); }
    void OptionalElements (const char *oe) { _optionalElements.put(oe); }
};

class ListTypeDescriptor  :    public AggrTypeDescriptor  { 

  protected:
  public:  

/*    void AssignAggrCreator(ListAggregateCreator f = 0)
    {
	CreateNewAggr = f;
    }

    STEPaggregate *CreateListAggregate()
    {
	if(CreateNewAggr)
	    return CreateNewAggr();
	else
	    return 0;
    }
    */
    ListTypeDescriptor ( ) { }
    ListTypeDescriptor (const char * nm, PrimitiveType ft, const char * d, 
			AggregateCreator f =0 )
		: AggrTypeDescriptor (nm, ft, d, f) { }
    virtual ~ListTypeDescriptor () { }

};

class SetTypeDescriptor  :    public AggrTypeDescriptor  { 

  protected:
  public:  

    SetTypeDescriptor ( ) { } 
    SetTypeDescriptor (const char * nm, PrimitiveType ft, const char * d, 
			AggregateCreator f =0 )
		: AggrTypeDescriptor (nm, ft, d, f) { }
    virtual ~SetTypeDescriptor () { }

};

class BagTypeDescriptor  :    public AggrTypeDescriptor  { 

  protected:
  public:  

    BagTypeDescriptor ( ) { }
    BagTypeDescriptor (const char * nm, PrimitiveType ft, const char * d, 
			AggregateCreator f =0 )
		: AggrTypeDescriptor (nm, ft, d, f) { }
    virtual ~BagTypeDescriptor () { }

};

typedef  SdaiSelect * (* SelectCreator) () ;

class SelectTypeDescriptor  :    public TypeDescriptor  { 

  protected:
	TypeDescriptorList _elements ;	  //  of  TYPE_DESCRIPTOR
	int _unique_elements;

  public:  

	SelectCreator CreateNewSelect;

	void AssignSelectCreator(SelectCreator f = 0)
	{
	    CreateNewSelect = f;
	}

	SdaiSelect *CreateSelect()
	{
	    if(CreateNewSelect)
		return CreateNewSelect();
	    else
		return 0;
	}


        SelectTypeDescriptor (int b, const char * nm, PrimitiveType ft, 
			      char * d, SelectCreator f =0 ) 
          : TypeDescriptor (nm, ft, d), 
	  _unique_elements (b), CreateNewSelect(f)
		{ }
	virtual ~SelectTypeDescriptor () { }

	TypeDescriptorList& Elements() { return _elements; }
	const TypeDescriptorList& GetElements() const { return _elements; }
//	void Elements (TypeDescriptorList x);
        int UniqueElements () const {  return _unique_elements; }
        virtual const TypeDescriptor * IsA (const TypeDescriptor *) const;
        virtual const TypeDescriptor * IsA (const char * n) const  
                {  return TypeDescriptor::IsA (n);  }
        virtual const TypeDescriptor * CanBe (const TypeDescriptor *) const;
};

class StringTypeDescriptor  :    public TypeDescriptor  { 

  protected:
	SdaiInteger  _width ;    //  OPTIONAL
	SdaiLogical _fixedSize ;
  public:  

	StringTypeDescriptor ( ) : _fixedSize("UNKNOWN_TYPE") { _width = 0; }
	virtual ~StringTypeDescriptor () { }


	SdaiInteger Width()		{ return _width; }
	void Width (SdaiInteger  w)	{ _width = w; }

	SdaiLogical& FixedSize()		{ return _fixedSize; }
	void FixedSize (SdaiLogical fs)	{ _fixedSize.put(fs.asInt()); }
	void FixedSize (LOGICAL fs)	{ _fixedSize.put(fs); }
	void FixedSize (char * fs)	{ _fixedSize.put(fs); }
};

class RealTypeDescriptor  :    public TypeDescriptor  { 

  protected:
	SdaiInteger  _precisionSpec ;    //  OPTIONAL
  public:  

	RealTypeDescriptor ( ) { _precisionSpec = 0; }
	virtual ~RealTypeDescriptor () { }

	SdaiInteger PrecisionSpec()		{ return _precisionSpec; }
	void PrecisionSpec (SdaiInteger  ps)	{ _precisionSpec = ps; }
};


#endif
