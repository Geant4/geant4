// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPcomplex.h,v 1.1 1999-01-07 16:08:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef STEPCOMPLEX_H
#define STEPCOMPLEX_H

#include <errordesc.h>
#include <sdai.h>
#include <baseType.h>
#include <ExpDict.h>
#include <STEPentity.h>
#include <Registry.h>

class STEPcomplex : public STEPentity {
  public:
    STEPcomplex * sc;
    STEPcomplex * head;
    Registry * _registry;
    int visited; // used when reading (or as you wish?)

  public:
    STEPcomplex(Registry *registry, int fileid);
    STEPcomplex(Registry *registry, const SCLstring **names, int fileid);
    STEPcomplex(Registry *registry, const char **names, int fileid);

    virtual ~STEPcomplex();

    int EntityExists(const char *name);
    STEPcomplex *EntityPart(const char *name);

/*
    // page 241 Stroustrup
    STEPcomplex &operator[](const char *name);
    STEPcomplex &operator[](const int index);
*/

    virtual Severity ValidLevel(ErrorDescriptor *error, InstMgr *im, 
			int clearError = 1);
// READ
    virtual Severity STEPread(int id, int addFileId, 
				   InstMgr * instance_set,
				   istream& in =cin);
    virtual void STEPread_error(char c, int index, istream& in);

// WRITE
    virtual void STEPwrite(ostream& out =G4cout, int writeComment = 1);
    virtual const char * STEPwrite(SCLstring &buf);

    virtual void WriteExtMapEntities(ostream& out =G4cout);
    virtual const char * WriteExtMapEntities(SCLstring &buf);
    virtual void AppendEntity(STEPcomplex *stepc);

  protected:
    virtual void CopyAs (STEPentity *);
    void BuildAttrs(const char *s );
    void AddEntityPart(const char *name);
    void AssignDerives();
    
};

#endif
