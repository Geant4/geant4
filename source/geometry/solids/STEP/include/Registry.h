

//



//
// $Id: Registry.h,v 1.2 1999-05-21 20:20:30 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef _REGISTRY_H
#define _REGISTRY_H

/*
* NIST STEP Core Class Library
* clstepcore/Registry.h
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#include <STEPentity.h>
#include <errordesc.h>
#include <scl_hash.h>
typedef struct Hash_Table * HashTable;

class Registry;
extern char * EntityClassName ( char *);
typedef void (* CF_init) (Registry &);	//  pointer to creation initialization 

class Registry {
  protected:
    HashTable primordialSwamp;	//  dictionary of EntityDescriptors
    HashTable active_schemas;	//  dictionary of SchemaDescriptors
    HashTable active_types;	//  dictionary of TypeDescriptors
    
    int entity_cnt;
    HashEntry 	cur_entity;
    HashEntry 	cur_schema;

  public:
    Registry (CF_init initFunct);
    ~Registry ();
    void DeleteContents ();  // CAUTION: calls delete on all the descriptors 

    const EntityDescriptor* FindEntity (const char *, int check_case =1) const;
    const SchemaDescriptor* FindSchema (const char *, int check_case =1) const;
    const TypeDescriptor*   FindType (const char *, int check_case =1) const;
    
    void 	AddEntity (const EntityDescriptor&);
    void	AddSchema (const SchemaDescriptor&);
    void 	AddType  (const TypeDescriptor&);
    
    void 	RemoveEntity (const char *);
    void	RemoveSchema (const char *);
    void 	RemoveType (const char *);

    int 	GetEntityCnt ();
    void	ResetEntities ();
    const EntityDescriptor *	NextEntity ();
    
    void	ResetSchemas ()  ;
    const SchemaDescriptor *	NextSchema ();

    STEPentity* ObjCreate (const char * nm, int check_case =1) const;
    
};

#endif


