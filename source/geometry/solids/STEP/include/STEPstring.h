

//



//
// $Id: STEPstring.h,v 1.2 1999-05-21 20:20:32 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef STEPSTRING_H
#define	STEPSTRING_H  1

/*
* NIST STEP Core Class Library
* clstepcore/STEPstring.h
* May 1995
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

class ErrorDescriptor;
#include <scl_string.h>
#include <errordesc.h>

#ifndef STRING_DELIM
#define STRING_DELIM '\''
#endif

class SdaiString : public SCLstring {
public:

  //constructor(s) & destructor    
  SdaiString (const char * str = 0, int max =0) : SCLstring (str,max) { }
  SdaiString (const SCLstring& s)   : SCLstring (s) { }
  SdaiString (const SdaiString& s)  : SCLstring (s) { }
  ~SdaiString ()  {  }

//  operators
  SdaiString& operator= (const char* s);

  // format for STEP
  const char * asStr (SCLstring & s) const  {  return s = chars ();  }
  void STEPwrite (ostream& out =G4cout)  const;
  void STEPwrite (SCLstring &s) const;

  Severity StrToVal (const char *s);
  Severity STEPread (istream& in, ErrorDescriptor *err);
  Severity STEPread (const char *s, ErrorDescriptor *err);

 protected:
};

inline
SdaiString& 
SdaiString::operator= (const char* s)
    { SCLstring::operator= (s);
      return *this;  }

#endif
