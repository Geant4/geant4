// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPaggregate.h,v 1.1 1999-01-07 16:08:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef STEPAGGREGATE_H
#define STEPAGGREGATE_H

/*
* NIST STEP Core Class Library
* clstepcore/STEPaggregate.h
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

class InstMgr;
class STEPaggregate;
class Logical;
class Boolean;
class TypeDescriptor;

#include "globals.hh"
#include <errordesc.h>
#include <SingleLinkList.h>
#include <sdai.h>
#include <baseType.h>
#include <SdaiBinary.h>
#include <STEPundefined.h>
#include <scl_string.h>

class STEPentity;
#define 	S_ENTITY_NULL	&NilSTEPentity
extern STEPentity NilSTEPentity;


#define 	AGGR_NULL	&NilSTEPaggregate
extern STEPaggregate NilSTEPaggregate;

typedef unsigned short BOOLEAN;
class SingleLinkNode;

/******************************************************************************
 **
 *****************************************************************************/

typedef STEPaggregate * STEPaggregateH;
class STEPaggregate :  public SingleLinkList 
{
  protected:
    int _null;

  protected:

    virtual Severity ReadValue(istream &in, ErrorDescriptor *err, 
			       const TypeDescriptor *elem_type, 
			       InstMgr *insts, int addFileId =0, 
			       int assignVal =1, int ExchangeFileFormat =1);
  public:

    int is_null() { return _null; } 

    virtual Severity AggrValidLevel(const char *value, ErrorDescriptor *err, 
			    const TypeDescriptor *elem_type, InstMgr *insts, 
			    int optional, char *tokenList, int addFileId =0, 
			    int clearError =0);
    

    virtual Severity AggrValidLevel(istream &in, ErrorDescriptor *err, 
			    const TypeDescriptor *elem_type, InstMgr *insts, 
			    int optional, char *tokenList, int addFileId =0, 
			    int clearError =0);

// INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err =0, 
			      const TypeDescriptor *elem_type =0,
			      InstMgr *insts =0, int addFileId = 0);
    virtual Severity STEPread (istream& in, ErrorDescriptor *err, 
			       const TypeDescriptor *elem_type =0,
			       InstMgr *insts =0, int addFileId =0);
// OUTPUT
    virtual const char *asStr(SCLstring & s) const;
    virtual void STEPwrite  (ostream& out =G4cout) const;

//    SingleLinkNode * GetHead () const
//	{ return (STEPnode *) SingleLinkList::GetHead(); }

    virtual SingleLinkNode * NewNode ();  
    void AddNode (SingleLinkNode *);
    void Empty ();
    
    STEPaggregate ();
    virtual ~STEPaggregate ();

// COPY - defined in subtypes
    virtual STEPaggregate& ShallowCopy (const STEPaggregate&);

};									      

/******************************************************************
 ** Class:  GenericAggregate
 ** Description:  This class supports LIST OF: 
 **	SELECT_TYPE, BINARY_TYPE, GENERIC_TYPE, ENUM_TYPE, UNKNOWN_TYPE type
 ******************************************************************/

class GenericAggregate  :  public STEPaggregate 
{
  public:
    virtual SingleLinkNode *NewNode();
    virtual STEPaggregate& ShallowCopy (const STEPaggregate&);

    GenericAggregate  () {     }
    virtual ~GenericAggregate () {     }
};
typedef  GenericAggregate * GenericAggregateH;

/******************************************************************************
 **
 *****************************************************************************/

class EntityAggregate  :  public  STEPaggregate 
{
  public:
    virtual Severity ReadValue(istream &in, ErrorDescriptor *err, 
			       const TypeDescriptor *elem_type, 
			       InstMgr *insts, int addFileId =0, 
			       int assignVal =1, int ExchangeFileFormat =1);

    virtual  SingleLinkNode *NewNode();
    virtual  STEPaggregate&  ShallowCopy (const STEPaggregate&);

    EntityAggregate ();
    virtual ~EntityAggregate ();
};
typedef   EntityAggregate * EntityAggregateH;

/******************************************************************
 ** Class:  SelectAggregate
 ** Description:  this is a minimal representions for a collection of 
 **	STEPenumerations
 ******************************************************************/
class SelectAggregate  :  public STEPaggregate 
{
  public:
    virtual Severity ReadValue(istream &in, ErrorDescriptor *err, 
			       const TypeDescriptor *elem_type, 
			       InstMgr *insts, int addFileId =0, 
			       int assignVal =1, int ExchangeFileFormat =1);

    virtual SingleLinkNode *NewNode ();
    virtual STEPaggregate& ShallowCopy  (const STEPaggregate&);

    SelectAggregate ();
    virtual ~SelectAggregate ();
};
typedef  SelectAggregate *  SelectAggregateH;

/******************************************************************
 ** Class:  StringAggregate
 ** Description:  This class supports LIST OF STRING type
 ******************************************************************/
class StringAggregate  :  public STEPaggregate 
{
  public:
    virtual SingleLinkNode *NewNode();
    virtual STEPaggregate& ShallowCopy (const STEPaggregate&);

    StringAggregate  () {     }
    virtual ~StringAggregate () {     }
};
typedef  StringAggregate * StringAggregateH;


/******************************************************************
 ** Class:  BinaryAggregate
 ** Description:  This class supports LIST OF BINARY type
 ******************************************************************/
class BinaryAggregate  :  public STEPaggregate 
{
  public:
    virtual SingleLinkNode *NewNode();
    virtual STEPaggregate& ShallowCopy (const STEPaggregate&);

    BinaryAggregate  () {     }
    virtual ~BinaryAggregate () {     }
};
typedef  BinaryAggregate * BinaryAggregateH;

/******************************************************************
 ** Class:  EnumAggregate
 ** Description:  this is a minimal representions for a collection of 
 **	STEPenumerations
 ******************************************************************/
class EnumAggregate  :  public STEPaggregate 
{
  public:
    virtual SingleLinkNode *NewNode ();
    virtual STEPaggregate& ShallowCopy  (const STEPaggregate&);

    EnumAggregate ();
  virtual ~EnumAggregate ();
};
typedef  EnumAggregate *  EnumAggregateH;

class Logicals  : public EnumAggregate  
{
  public:
    virtual SingleLinkNode * NewNode ();  

    Logicals  () { }
    virtual ~Logicals () { }
};
typedef  Logicals *  LogicalsH;
inline Logicals * create_Logicals () { return new Logicals ; }


class Booleans  : public EnumAggregate  
{
  public:
    virtual SingleLinkNode * NewNode ();  

    Booleans  () { }
    virtual ~Booleans () { }
};
typedef  Booleans *  BooleansH;
inline Booleans * create_Booleans () { return new Booleans ; }

class RealAggregate  : public STEPaggregate  {

  public:
    virtual SingleLinkNode *NewNode();
    virtual STEPaggregate& ShallowCopy  (const STEPaggregate&);

    RealAggregate  () {     }
    virtual ~RealAggregate () {     }
};
typedef  RealAggregate *  RealAggregateH;

class IntAggregate  : public STEPaggregate  {

  public:
    virtual SingleLinkNode *NewNode();
    virtual STEPaggregate& ShallowCopy  (const STEPaggregate&);

    IntAggregate  () {     }
    virtual ~IntAggregate () {     }
};
typedef  IntAggregate *  IntAggregateH;

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

class STEPnode :  public SingleLinkNode  {
  protected:
    int _null;

  public:
    int is_null() { return _null; } 
    void set_null() { _null = 1; } 

//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite(SCLstring &s);
    virtual void STEPwrite (ostream& out =G4cout);
};
typedef  STEPnode *  STEPnodeH;

/******************************************************************
 ** Class:  GenericNode
 ** Description:  This class is for the Nodes of GenericAggregates
 ******************************************************************/

class GenericAggrNode  : public STEPnode {
  public:

    SCLundefined value;

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite(SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    GenericAggrNode (const char *str);
    GenericAggrNode (GenericAggrNode& gan);
    GenericAggrNode ();
    ~GenericAggrNode ();

    virtual SingleLinkNode *	NewNode ();

};

/////////////////////////////////////////////////////////////////////////////

class EntityNode  : public STEPnode {
  public:
    STEPentity * node;

// INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);
//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    EntityNode (STEPentity * e);
    EntityNode ()  {    }
    ~EntityNode ()  {    }

    virtual SingleLinkNode *	NewNode ();

    // Calling these funtions is an error.
    Severity StrToVal(const char *s, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return StrToVal(s, err, 0, 0, 0);
    }
    Severity StrToVal(istream &in, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return StrToVal(in, err, 0, 0, 0);
    }

    Severity STEPread(const char *s, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return STEPread(s, err, 0, 0, 0);
    }
    Severity STEPread(istream &in, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return STEPread(in, err, 0, 0, 0);
    }
};

///////////////////////////////////////////////////////////////////////////

/******************************************************************
 ** Class:  SelectNode
 ** Description:  this is a minimal representions for node in lists of
 **	STEPenumerations
 ******************************************************************/

class SelectNode  : public STEPnode {
  public:

    SdaiSelect * node;

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type,
			      InstMgr *insts, int addFileId = 0);
//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    SelectNode (SdaiSelect * s) :  node (s) {    }
    SelectNode ()  {    }
    ~SelectNode ()  {    }

    virtual SingleLinkNode *	NewNode ();

    // Calling these funtions is an error.
    Severity StrToVal(const char *s, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return StrToVal(s, err, 0, 0, 0);
    }
    Severity StrToVal(istream &in, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return StrToVal(in, err, 0, 0, 0);
    }

    Severity STEPread(const char *s, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return STEPread(s, err, 0, 0, 0);
    }
    Severity STEPread(istream &in, ErrorDescriptor *err)
    {
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	    << "\n" << _POC_ "\n";
	return STEPread(in, err, 0, 0, 0);
    }
};

///////////////////////////////////////////////////////////////////////////

/******************************************************************
 ** Class:  StringNode
 ** Description:  This class is for the Nodes of StringAggregates
 ******************************************************************/

class StringNode  : public STEPnode {
  public:

    SdaiString value;

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    StringNode(StringNode& sn);
    StringNode (const char * sStr);
    StringNode ()   { value = 0; }
    ~StringNode ()  { }

    virtual SingleLinkNode *	NewNode ();
};

///////////////////////////////////////////////////////////////////////////

/******************************************************************
 ** Class:  BinaryNode
 ** Description:  This class is for the Nodes of BinaryAggregates
 ******************************************************************/

class BinaryNode  : public STEPnode {
  public:

    SdaiBinary value;

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    BinaryNode(BinaryNode& bn);
    BinaryNode (const char * sStr);
    BinaryNode ()   { value = 0; }
    ~BinaryNode ()  { }

    virtual SingleLinkNode *	NewNode ();
};

/******************************************************************
 ** Class:  EnumNode
 ** Description:  this is a minimal representions for node in lists of
 **	STEPenumerations
 ******************************************************************/

class EnumNode  : public STEPnode {
  public:

    STEPenumeration * node;

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    EnumNode (STEPenumeration * e) :  node (e) {    }
    EnumNode ()  {    }
    ~EnumNode ()  {    }

    virtual SingleLinkNode *	NewNode ();
};

class RealNode  : public STEPnode {
  public:
    SdaiReal value; // double

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    RealNode ()  { value = S_REAL_NULL; }
    ~RealNode ()  {    }

    virtual SingleLinkNode *	NewNode ();
};

class IntNode  : public STEPnode {
  public:
    SdaiInteger value; // long int

  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring & s);
    virtual const char *STEPwrite (SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

//	CONSTRUCTORS
    IntNode ()  { value = S_INT_NULL; }
    ~IntNode ()  {    }

    virtual SingleLinkNode *	NewNode ();
};

/******************************************************************
 **	The following classes are currently stubs
 **							
**/

class Array  :  public STEPaggregate  {	
  public:								      
    int lowerBound;							      
    int upperBound;							      
};
									      
class Bag  :  public STEPaggregate  {								      
  public:								      
    int min_ele;							      
    int max_ele;							      
};
									      
class List :  public STEPaggregate  {								      
    int min_ele;							      
    int max_ele;							      
    List *prev;								      
};
									      
class Set :								      
public STEPaggregate  {								      
    int cnt;								      
};

#endif 
