

//



//
// $Id: SdaiBinary.h,v 1.3 1999-12-15 14:50:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef SDAIBINARY_H
#define	SDAIBINARY_H 1

/*
* NIST STEP Core Class Library
* clstepcore/SdaiBinary.h
* May 1995
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/* $Id: SdaiBinary.h,v */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "g4std/strstream"

class ErrorDescriptor;
#include <scl_string.h>
#include <errordesc.h>

#ifndef BINARY_DELIM
#define BINARY_DELIM '\"'
#endif

class SdaiBinary : public SCLstring
{
  public:

    //constructor(s) & destructor    
    SdaiBinary (const char * str = 0, int max =0) : SCLstring (str,max) { }

//Josh L, 3/28/95
//    SdaiBinary (SCLstring& s)   : SCLstring (s) { }
//    SdaiBinary (SdaiBinary& s)  : SCLstring (s) { }
    SdaiBinary (const SCLstring& s)   : SCLstring (s) { }


    ~SdaiBinary ()  {  }

    //  operators
    SdaiBinary& operator= (const char* s);

    // format for STEP
    const char * asStr () const  {  return chars ();  }
    void STEPwrite (G4std::ostream& out =G4cout)  const;
    const char * STEPwrite (SCLstring &s) const;

    Severity StrToVal (const char *s, ErrorDescriptor *err);
    Severity STEPread (G4std::istream& in, ErrorDescriptor *err);
    Severity STEPread (const char *s, ErrorDescriptor *err);

    Severity BinaryValidLevel (const char *value, ErrorDescriptor *err,
			       int optional, char *tokenList,
			       int needDelims = 0, int clearError = 1);
    Severity BinaryValidLevel (G4std::istream &in, ErrorDescriptor *err, 
			       int optional, char *tokenList,
			       int needDelims = 0, int clearError = 1);

 protected:
  Severity ReadBinary(G4std::istream& in, ErrorDescriptor *err, int AssignVal = 1,
		      int needDelims = 1);
};

inline
SdaiBinary& 
SdaiBinary::operator= (const char* s)
    { SCLstring::operator= (s);
      return *this;  }

#endif
