

//



//
// $Id: scl_string.h,v 1.3 1999-12-15 18:04:17 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef _SCL_STRING_H
#define _SCL_STRING_H

/*
* NIST Utils Class Library
* clutils/scl_string.h
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

/*  point of contact for bug reports  */
#define _POC_  " report problem to dp2@cme.nist.gov "

//#include <std.h> // not found in CenterLine C++
#ifdef __OBJECTCENTER__
// this file is not in gnu C++ but it doesn't seem to be needed.
#include <stdarg.h>
#endif

#include "G4ios.hh"
#include <string.h>

/******************************************************************
 ** Class:  SCLstring
** Description: implements a few basic string handling functions - 
**              hopefully will be replaced by a standard clas
 ** Status:  26-Jan-1994 kc
 ******************************************************************/

class SCLstring   {
protected:
  char * _strBuf;   //  initially empty
  int _strBufSize;  // size of buffer
  int _max_len;  // should be const, but some db\'s don\'t handle that

  int newBufSize (int len) const;

public:
	// returns 1 if _strBuf is a null ptr or if it is an empty string ("")
    int is_null() const;
	// returns 1 if _strBuf is a null ptr, and 0 otherwise
    int is_undefined() const;
	// deletes _strBuf
    void set_undefined() ;
	// sets _strBuf space to be zeroed out
    int set_null() ;
    int StrBufSize() const;
    int Length() const;
    int MaxLength() const	{ return _max_len; }

    const char * chars () const
	{ return _strBuf ? _strBuf : ""; }

    const char * rep() const
	{ return _strBuf; }

//  operators
    SCLstring& operator= (const char*); // must be NULL terminated string
    SCLstring& operator= (char* s) // must be NULL terminated string
	{ operator= ((const char *)s); return *this; }
    operator const char * () const;
    int operator== (const char*) const;
    SCLstring& Append (const char *);
    SCLstring& Append (const char);
    SCLstring& Append (const long int);
    SCLstring& Append (const int);

// Josh L, 4/11/95
    SCLstring& operator=(const SCLstring& s);    

    SCLstring& Append (const double, const int precision =15);
    SCLstring& Prepend (const char *);

//constructor(s) & destructor    
    SCLstring (const char * str = 0, int max =0);
    SCLstring (const SCLstring& s);
    ~SCLstring ();

};

#endif   


