

//



//
// $Id: STEPcomplex.h,v 1.3 1999-12-15 14:50:14 gunter Exp $
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
				   G4std::istream& in =G4cin);
    virtual void STEPread_error(char c, int index, G4std::istream& in);

// WRITE
    virtual void STEPwrite(G4std::ostream& out =G4cout, int writeComment = 1);
    virtual const char * STEPwrite(SCLstring &buf);

    virtual void WriteExtMapEntities(G4std::ostream& out =G4cout);
    virtual const char * WriteExtMapEntities(SCLstring &buf);
    virtual void AppendEntity(STEPcomplex *stepc);

  protected:
    virtual void CopyAs (STEPentity *);
    void BuildAttrs(const char *s );
    void AddEntityPart(const char *name);
    void AssignDerives();
    
};

#endif
