// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPundefined.h,v 1.1 1999-01-07 16:08:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef	STEPUNDEFINED_H
#define	STEPUNDEFINED_H

/*
* NIST STEP Core Class Library
* clstepcore/STEPundefined.h
* May 1995
* KC Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <errordesc.h>
#include <scl_string.h>
#include <read_func.h>

class SCLundefined  {
  protected:
    SCLstring val;
    
  public:
//	INPUT
    virtual Severity StrToVal(const char *s, ErrorDescriptor *err);
    virtual Severity StrToVal(istream &in, ErrorDescriptor *err);

    virtual Severity STEPread(const char *s, ErrorDescriptor *err);
    virtual Severity STEPread(istream &in, ErrorDescriptor *err);

//	OUTPUT
    virtual const char *asStr(SCLstring &s) const;
    virtual const char *STEPwrite(SCLstring &s);
    virtual void 	STEPwrite (ostream& out =G4cout);

    int set_null ();
    int is_null ();
    SCLundefined& operator= (const SCLundefined&); 
    SCLundefined& operator= (const char *str); 
    SCLundefined ();
    virtual ~SCLundefined ();
}
;

#endif
