

//



//
// $Id: STEPundefined.h,v 1.2 1999-05-21 20:20:33 japost Exp $
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
