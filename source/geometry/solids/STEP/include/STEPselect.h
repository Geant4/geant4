// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPselect.h,v 1.1 1999-01-07 16:08:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef _STEPSELECT_H
#define _STEPSELECT_H

/*
* NIST STEP Core Class Library
* clstepcore/STEPselect.h
* May 1995
* Dave Helfrick
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <baseType.h>
#include <scl_string.h>
#include <sdai.h>
#include <errordesc.h>
#include <read_func.h>

class TypeDescriptor;
class SelectTypeDescriptor;
class InstMgr;

/**********
	class definition for the select superclass SdaiSelect.
**********/
class SdaiSelect {
  protected:
        const SelectTypeDescriptor *_type;
        const TypeDescriptor *      underlying_type;
	BASE_TYPE 		    base_type; // used by the subtypes

	SdaiString val;
	ErrorDescriptor _error;
        const TypeDescriptor * SetUnderlyingType (const TypeDescriptor *);

        const TypeDescriptor * CanBe (const char *) const;
        const TypeDescriptor * CanBe (BASE_TYPE) const;
	const TypeDescriptor * CanBe (const TypeDescriptor * td) const;

	virtual const TypeDescriptor * AssignEntity (STEPentity * se) =0;
	virtual SdaiSelect * NewSelect () =0;
  public:
	Severity severity() const;
	Severity severity( Severity );
	const char *Error();
	void Error( char * );
		// clears select's error  
	void ClearError();
		// clears error

  // constructors
        SdaiSelect (const SelectTypeDescriptor * s =0, 
		     const TypeDescriptor * td =0) 
		     : _type (s), underlying_type (td) { }

	virtual ~SdaiSelect ()	{  };

  // from SDAI binding
        SdaiString UnderlyingTypeName () const;
	const TypeDescriptor * CurrentUnderlyingType() const;
	int exists() const;
	void nullify();

	Severity SelectValidLevel(const char *attrValue, ErrorDescriptor *err, 
				  InstMgr *im, int clearError);

  // reading and writing
        const char * STEPwrite(SCLstring& s)  const;
	void STEPwrite (ostream& out =G4cout) const;
        virtual void STEPwrite_content (ostream& out) const =0;
//	char * asStr() const;


	Severity StrToVal(const char *val, const char *selectType, 
			  ErrorDescriptor *err, InstMgr * instances =0);
        virtual Severity StrToVal_content (const char *, 
					   InstMgr * instances =0) =0;

	Severity STEPread(istream& in, ErrorDescriptor *err, 
			  InstMgr * instances = 0, int addFileId =0);

		// abstract function
        virtual Severity STEPread_content (istream& in =cin, 
					   InstMgr * instances =0, 
					   int addFileId =0) =0;

	virtual SdaiSelect& operator =( const SdaiSelect& ) { return *this; } 

	int set_null();
	int is_null();
};	/** end class  **/

typedef SdaiSelect * SdaiSelectH ;

#endif
