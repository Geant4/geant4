// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPentity.h,v 1.1 1999-01-07 16:08:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef STEPENTITY_H
#define	STEPENTITY_H 1

/*
* NIST STEP Core Class Library
* clstepcore/STEPentity.h
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
#include <STEPattributeList.h>

#include <EntityInst.h>

#include <sdai.h>

class STEPattributeList;
class STEPattribute;

#include <ctype.h>
#include <Str.h>

class InstMgr;
class EntityDescriptor;

///////////////////////////////////////////////////////////////////////////////
// STEPentity is like the corresponding SDAI object AppInstance.

class STEPentity : public EntityInstance {
  private:
    int _cur;	// provides a built-in way of accessing attributes in order.

 public:
    STEPattributeList attributes;
    int 	      STEPfile_id;
    ErrorDescriptor   _error;
    SCLstring	      *p21Comment;
					// registry additions
    EntityDescriptor *eDesc;

	// head entity for multiple inheritance.  If it is null then this 
	// STEPentity is not part of a multiply inherited entity.  If it 
	// points to a STEPentity then this STEPentity is part of a mi entity
	// and head points at the root STEPentity of the primary inheritance 
	// path (the one that is the root of the leaf entity).
    STEPentity *headMiEntity;
	// these form a chain of other entity parents for multiple inheritance
    STEPentity *nextMiEntity;

  protected:
    int _complex;

 public:
    STEPentity ();
    STEPentity (int fileid, int complex = 0);
    virtual ~STEPentity();

    int IsComplex() { return _complex; }
    int SetFileId(int fid) { return STEPfile_id = fid; }
    int GetFileId() const  { return STEPfile_id; }
    int FileId (int fid) { return STEPfile_id = fid; }
    int FileId() const  { return STEPfile_id; }

    void AddP21Comment(SCLstring &s, int replace = 1);
    void AddP21Comment(const char *s, int replace = 1);
    void DeleteP21Comment() { delete p21Comment; p21Comment = 0; }

    // guaranteed a string (may be null string)
    const char *P21Comment() 
	{ return ( p21Comment ? p21Comment->chars() : "" ); }
    // returns null if no comment exists
    const char *P21CommentRep() 
	{ return ( p21Comment ? p21Comment->rep() : 0 ); }

    const char *EntityName() const { return eDesc->Name(); }

    virtual Severity ValidLevel(ErrorDescriptor *error, InstMgr *im, 
			int clearError = 1);
    ErrorDescriptor &Error()	{ return _error; }
		// clears entity's error and optionally all attr's errors
    void ClearError(int clearAttrs = 1);
		// clears all attr's errors
    void ClearAttrError();
//    void EnforceOptionality(int on = 1);

    STEPentity *Replicate();

// ACCESS attributes in order.
    int AttributeCount();
    STEPattribute * NextAttribute();
    void ResetAttributes()		{ _cur =0; }
    
// READ
    virtual Severity STEPread(int id, int addFileId, 
			      class InstMgr * instance_set,
			      istream& in =cin);
    virtual void STEPread_error(char c, int index, istream& in);

// WRITE
    virtual void STEPwrite(ostream& out =G4cout, int writeComment = 1);
    virtual const char * STEPwrite(SCLstring &buf);

    void	 STEPwrite_reference (ostream& out =G4cout);
    const char * STEPwrite_reference (SCLstring &buf);

    void beginSTEPwrite(ostream& out =G4cout); // writes out the SCOPE section
    void endSTEPwrite(ostream& out =G4cout);

// MULTIPLE INHERITANCE
    int MultipleInheritance() { return !(headMiEntity == 0); }

    void HeadEntity(STEPentity *se) { headMiEntity = se; }
    STEPentity * HeadEntity() { return headMiEntity; }

    STEPentity *GetNextMiEntity() { return nextMiEntity; }
    STEPentity *GetMiEntity(char *EntityName);
    void AppendMultInstance(STEPentity *se);

 protected:
    STEPattribute * GetSTEPattribute (const char *);
    STEPattribute * MakeDerived (const char *);

    virtual void CopyAs (STEPentity *);
    void PrependEntityErrMsg();
}
;

///////////////////////////////////////////////////////////////////////////////

extern STEPentity NilSTEPentity;
#define ENTITY_NULL	&NilSTEPentity
#define NULL_ENTITY	&NilSTEPentity

typedef STEPentity* STEPentityPtr;
typedef STEPentity* STEPentityH;

extern STEPentity *
ReadEntityRef(istream &in, ErrorDescriptor *err, char *tokenList, 
	      InstMgr * instances, int addFileId);

//typedef  STEPentity * (* Creator) () const;
//typedef  STEPentity * (* Creator) () ;

#endif
