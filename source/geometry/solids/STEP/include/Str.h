// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Str.h,v 1.1 1999-01-07 16:08:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef STR_H
#define STR_H

/*
* NIST Utils Class Library
* clutils/Str.h
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

#include <ctype.h>

//#include <std.h> // not found in CenterLine C++
// the two includes stdio.h and stdlib.h below are replacing std.h since 
// CenterLine doesn't have std.h

#include <stdio.h> // added to have the definition for BUFSIZE
#include <stdlib.h> 
#include <string.h>
#include <scl_string.h>
#include <errordesc.h>

char ToLower (const char c);
char ToUpper  (const char c);
const char * StrToLower (const char * word, SCLstring &s);
const char * StrToUpper (const char * word, SCLstring &s);
const char * StrToConstant (const char * word, SCLstring &s);
const char * PrettyTmpName (const char * oldname);
char * PrettyNewName (const char * oldname);
char * EntityClassName ( char * oldname);

extern Severity CheckRemainingInput
   (istream &in, ErrorDescriptor *err, 
    const char *typeName, // used in error message
    const char *tokenList); // e.g. ",)"


#endif 
