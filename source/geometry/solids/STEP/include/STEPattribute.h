

//



//
// $Id: STEPattribute.h,v 1.3 1999-12-15 14:50:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef STEPATTRIBUTE_H
#define	STEPATTRIBUTE_H	1

/*
* NIST STEP Core Class Library
* clstepcore/STEPattribute.h
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <stdio.h>
#include <errordesc.h>
#include <baseType.h>

// this is used to set a const int Real_Num_Precision 
// in STEPaggregate.cc and STEPattribute.cc
#define REAL_NUM_PRECISION 15

typedef unsigned short BOOLEAN;
typedef double real;  

class InstMgr;
class STEPentity;
class STEPenumeration;
class STEPaggregate;
class SCLundefined;
class SdaiSelect;
class SdaiBinary;

class TypeDescriptor;
class AttrDescriptor;
class EntityDescriptor;

#include "g4std/strstream"
#include <ExpDict.h>

#define s_String	char *

extern int SetErrOnNull(const char *attrValue, ErrorDescriptor *error);
////////////////////
////////////////////

extern Severity 
CheckRemainingInput(G4std::istream &in, ErrorDescriptor *err, 
		    const char *typeName, // used in error message
		    const char *tokenList); // e.g. ",)"

extern STEPentity *
ReadEntityRef(G4std::istream &in, ErrorDescriptor *err, char *tokenList, 
	      InstMgr * instances, int addFileId);

extern STEPentity *
ReadEntityRef(const char * s, ErrorDescriptor *err, char *tokenList, 
	      InstMgr * instances, int addFileId);

extern Severity 
EntityValidLevel(STEPentity *se, 
		 const TypeDescriptor *ed, // entity type that entity se needs 
					   // to match. (this must be an
					   // EntityDescriptor)
		 ErrorDescriptor *err);

extern Severity 
EntityValidLevel(const char *attrValue, // string contain entity ref
		 const TypeDescriptor *ed, // entity type that entity in 
					   // attrValue (if it exists) needs 
					   // to match. (this must be an
					   // EntityDescriptor)
		 ErrorDescriptor *err, InstMgr *im, int clearError);

////////////////////
////////////////////

extern STEPentity *STEPread_reference (const char * s, ErrorDescriptor *err, 
				       InstMgr * instances, int addFileId);
////////////////////

extern int   QuoteInString(G4std::istream& in);

extern void  AppendChar(char c, int& index, char *&s, int& sSize);

extern void 
PushPastString (G4std::istream& in, SCLstring &s, ErrorDescriptor *err);

extern void 
PushPastImbedAggr (G4std::istream& in, SCLstring &s, ErrorDescriptor *err);

extern void 
PushPastAggr1Dim(G4std::istream& in, SCLstring &s, ErrorDescriptor *err);

//extern  Severity ValidateEntityType(STEPentity *se, 
//					const AttrDescriptor *ad, 
//					ErrorDescriptor *error);

class STEPattribute {

    friend G4std::ostream &operator<< ( G4std::ostream&, STEPattribute& );
    friend class STEPentity;
    
  protected:
    ErrorDescriptor _error;
    unsigned int _derive : 1;
    int Derive (unsigned int n =1)  { return _derive =n; }

  public:
    const AttrDescriptor * aDesc;

    // You know which of these to use based on the return value of
    // NonRefType() - see below. BASE_TYPE is defined in baseType.h
    // This variable points to an appropriate member variable in the entity
    // class in the generated schema class library (the entity class is 
    // inherited from STEPentity)
    union  {
	SdaiInteger *i;		// INTEGER_TYPE // SdaiInteger is a long int
	class SdaiString *S;	// STRING_TYPE
	class SdaiBinary *b;	// BINARY_TYPE
	SdaiReal *r;	   // REAL_TYPE and NUMBER_TYPE // SdaiReal is a double
	class STEPentity* *c;	// ENTITY_TYPE
	STEPaggregate *a;	// AGGREGATE_TYPE
	STEPenumeration *e;	// ENUM_TYPE, BOOLEAN_TYPE, and LOGICAL_TYPE
	class SdaiSelect *sh;	// SELECT_TYPE
	SCLundefined *u;	// UNKNOWN_TYPE

	void *p;
	
	} ptr;

  protected:
    char SkipBadAttr(G4std::istream& in, char *StopChars);
    void AddErrorInfo();

  public:

///////////// Read, Write, Assign attr value

    Severity StrToVal(const char *s, InstMgr *instances =0, 
		      int addFileId =0);
    Severity STEPread(G4std::istream& in = G4cin, InstMgr *instances =0, 
		      int addFileId =0);

    const char * asStr(SCLstring &) const; // return the attr value as a string
    void STEPwrite(G4std::ostream& out = G4cout);

    BOOLEAN ShallowCopy(STEPattribute *sa);

    Severity set_null();

////////////// Return info on attr

    BOOLEAN	Nullable() const; // may this attribute be null?
    BOOLEAN	is_null () const; // is this attribute null?
    int 	IsDerived () const  {  return _derive;  }

    const s_String 	Name() const;
    const s_String	TypeName() const;
    const BASE_TYPE	Type() const;
    const BASE_TYPE	NonRefType() const;
    const BASE_TYPE	BaseType() const;

    const TypeDescriptor   *ReferentType() const;

    ErrorDescriptor &Error()	{ return _error; }
    void ClearErrorMsg()	{ _error.ClearErrorMsg(); } 

    Severity ValidLevel (const char *attrValue, ErrorDescriptor *error, 
			     InstMgr *im, int clearError = 1);
  public:

////////////////// Constructors

   STEPattribute (const STEPattribute& a);
   STEPattribute ()  {};

   ~STEPattribute () {}; 

	   //  INTEGER
   STEPattribute (const class AttrDescriptor& d, SdaiInteger *p);
	   //  BINARY
   STEPattribute (const class AttrDescriptor& d, SdaiBinary *p);
	   //  STRING
   STEPattribute (const class AttrDescriptor& d, SdaiString *p);
	   //  REAL & NUMBER
   STEPattribute (const class AttrDescriptor& d, SdaiReal *p);
	   //  ENTITY
   STEPattribute (const class AttrDescriptor& d, STEPentity* *p);
	   //  AGGREGATE
   STEPattribute (const class AttrDescriptor& d, STEPaggregate *p);
	   //  ENUMERATION  and Logical
   STEPattribute (const class AttrDescriptor& d, STEPenumeration *p);
	   //  SELECT
   STEPattribute (const class AttrDescriptor& d, class SdaiSelect *p);
	   //  UNDEFINED
   STEPattribute (const class AttrDescriptor& d, SCLundefined *p);

  friend int operator == (STEPattribute &a1, STEPattribute &a2);
};

#endif
