

//



//
// $Id: Enumeration.h,v 1.2 1999-05-21 20:20:29 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef ENUMERATION_H
#define	ENUMERATION_H 

/*
* NIST STEP Core Class Library
* clstepcore/Enumeration.h
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

#include <Str.h>
#include <errordesc.h>

#define ENUM_NULL -1

class STEPenumeration  {
    friend     ostream &operator<< ( ostream&, const STEPenumeration& );
  protected:
    int v;	//  integer value of enumeration instance 
	//  mapped to a symbolic value in the elements

    int set_value (const char * n);
    int set_value (const int n);
    virtual ~STEPenumeration ()  {
    }

  public:
    
    virtual int no_elements () const =0;
    virtual const char * Name () const =0;
    const char * get_value_at (int n) const { return element_at (n);  }
    virtual const char * element_at (int n) const =0; 

    Severity EnumValidLevel(const char *value, ErrorDescriptor *err,
			    int optional, char *tokenList,
			    int needDelims = 0, int clearError = 1);

    Severity EnumValidLevel(istream &in, ErrorDescriptor *err, 
			    int optional, char *tokenList,
			    int needDelims = 0, int clearError = 1);

    const int asInt () const {	return v;    }
    
    const char * asStr (SCLstring &s) const;
    void STEPwrite (ostream& out = G4cout)  const;
    const char * STEPwrite (SCLstring &s) const;

    Severity StrToVal (const char * s, ErrorDescriptor *err, int optional = 1);
    Severity STEPread(istream& in, ErrorDescriptor *err, int optional = 1);
    Severity STEPread(const char *s, ErrorDescriptor *err, int optional = 1);

    int put (int val)  {  return set_value (val);    }
    int put (const char * n)	{ return set_value (n); }
    int is_null () const	{ return ( v == ENUM_NULL ); }
    int set_null()		{ return  put (ENUM_NULL); }
    STEPenumeration& operator= (const int);
    STEPenumeration& operator= (const STEPenumeration&);
        
    void DebugDisplay (ostream& out =G4cout) const;
  protected:
    Severity ReadEnum(istream& in, ErrorDescriptor *err, int AssignVal = 1,
		      int needDelims = 1);
};


enum LOGICAL { sdaiFALSE, sdaiTRUE, sdaiUNKNOWN };

#define F sdaiFALSE
#define U sdaiUNKNOWN
#define T sdaiTRUE

class Boolean  :
public STEPenumeration  {
  public:
    const char * Name () const  {
	return "Boolean";
    }

    Boolean (char * val =0)  {
	set_value (val);
    }
    Boolean (LOGICAL val);
    virtual ~Boolean ()  {  
    }
    inline virtual int no_elements () const {
	return 2;
    }
  operator LOGICAL () const;
  operator int () const;
  int operator ==( LOGICAL log ) const;

  virtual const char * element_at (int n) const;

}
;

inline Boolean * create_Boolean () { return new Boolean ; }

class Logical  :
public STEPenumeration  {
  public:
    const char * Name () const  {
	return "Logical";
    }

    Logical (char * val =0)  {
	set_value (val);
    }
    Logical (LOGICAL val)   {
	set_value (val);
    }

    // Josh L, 3/24/95
    Logical (const Boolean& boolean) {
      v = boolean.asInt();
    }

    virtual ~Logical ()  {  
    }
    inline virtual int no_elements () const {
	return 3;
    }
  operator LOGICAL () const;

  virtual const char * element_at (int n) const;

}
;

inline Logical * create_Logical () { return new Logical ; }

// Josh L, 3/24/95
// These constants are required for Part 23.  We declare them as statics
// inside a struct in order to prevent name conflicts when linking with
// other libraries.  This technique is discussed in "Effective C++" by 
// Scott Meyers, Addison-Wesley, 1992, pp. 93-95.

// GEANT4 workaround:
// TRUE and FALSE are sometimes defined as preprocessor symbols
// let's get rid of them until the struct has been parsed

struct SDAI {
  static const Logical UNKNOWN;
  static const Boolean SdaiTRUE;    
  static const Boolean SdaiFALSE;  
};

#endif 








