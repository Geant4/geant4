

//



//
// $Id: Registry_inline.cc,v 1.2 1999-05-21 20:20:48 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/Registry.inline.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#include <Registry.h>

 const TypeDescriptor *  t_INTEGER_TYPE;
 const TypeDescriptor *  t_REAL_TYPE;
 const TypeDescriptor *  t_NUMBER_TYPE;
 const TypeDescriptor *  t_STRING_TYPE;
 const TypeDescriptor *  t_BINARY_TYPE;
 const TypeDescriptor *  t_BOOLEAN_TYPE;
 const TypeDescriptor *  t_LOGICAL_TYPE;


/* inline */ 
Registry::Registry (CF_init initFunct) 
: entity_cnt (0)
{ 
    primordialSwamp = HASHcreate (1000);
    active_schemas = HASHcreate (10);
    active_types = HASHcreate (100);

    t_INTEGER_TYPE = new TypeDescriptor("INTEGER",     // Name
		       INTEGER_TYPE, // FundamentalType
		       "INTEGER");   // Description;

    t_REAL_TYPE = new TypeDescriptor("REAL", REAL_TYPE, "Real");

    t_STRING_TYPE = new TypeDescriptor("STRING", STRING_TYPE, "String");

    t_BINARY_TYPE = new TypeDescriptor("BINARY", BINARY_TYPE, "Binary");

    t_BOOLEAN_TYPE = new TypeDescriptor("BOOLEAN", BOOLEAN_TYPE, "Boolean");

    t_LOGICAL_TYPE = new TypeDescriptor("LOGICAL", LOGICAL_TYPE, "Logical");

    t_NUMBER_TYPE = new TypeDescriptor("NUMBER", NUMBER_TYPE, "Number");

/*    t_GENERIC_TYPE = new TypeDescriptor("GENERIC", GENERIC_TYPE, "Generic");*/

    initFunct (*this);  
    HASHlistinit (primordialSwamp, &cur_entity);  // Initialize cur\'s
    HASHlistinit (active_schemas, &cur_schema);
}

/* inline */ 
Registry::~Registry  ()
{
    HASHdestroy (primordialSwamp);
    HASHdestroy (active_schemas);
    HASHdestroy (active_types);
}

void 
Registry::DeleteContents ()
{
  // entities first
  HASHlistinit (primordialSwamp, &cur_entity);
  while (HASHlist (&cur_entity))
    delete (EntityDescriptor *) cur_entity.e->data;

  // schemas
  HASHlistinit (active_schemas, &cur_schema);
  while (HASHlist (&cur_schema))
    delete (SchemaDescriptor *) cur_schema.e->data;

  // types 

}

/* inline */ const EntityDescriptor *
Registry::FindEntity (const char * e, int check_case) const
{
    if (check_case) 
	return
	    (const EntityDescriptor *) HASHfind (primordialSwamp, 
						 (char *)PrettyTmpName (e));
    else return (const EntityDescriptor *) HASHfind (primordialSwamp, (char *) e);
}


/* inline */ const SchemaDescriptor *
Registry::FindSchema (const char * n, int check_case) const
{
    if (check_case) 
       return (const SchemaDescriptor *) HASHfind (primordialSwamp, 
						   (char *)PrettyTmpName (n));
    return (const SchemaDescriptor *) HASHfind (active_schemas, (char *) n);
}

/* inline */ const TypeDescriptor *
Registry::FindType (const char * n, int check_case) const 
{
    if (check_case) 
       return (const TypeDescriptor *) HASHfind (primordialSwamp, 
						 (char *)PrettyTmpName (n));
    return (const TypeDescriptor *) HASHfind (active_types, (char *) n);
}
    
/* inline */ void 	
Registry::AddEntity (const EntityDescriptor& e)
{
    HASHinsert (primordialSwamp, (char *) e.Name (), (EntityDescriptor *) &e);
    ++entity_cnt;
}
  

/* inline */ void	
Registry::AddSchema (const SchemaDescriptor& d)
{
    HASHinsert (active_schemas, (char *) d.Name(), (SchemaDescriptor *) &d);
}

/* inline */ void 	
Registry::AddType  (const TypeDescriptor& d)
{
    HASHinsert (active_types, (char *) d.Name(), (TypeDescriptor *) &d);
}
    
/* inline */ void 	
Registry::RemoveEntity (const char * n)
{
    struct Element tmp;
    tmp.key = (char *) n;
    HASHsearch (primordialSwamp, &tmp, HASH_DELETE) ? --entity_cnt : 0;
    
}

/* inline */ void	
Registry::RemoveSchema (const char * n)
{
    struct Element tmp;
    tmp.key = (char *) n;
    HASHsearch (active_schemas, &tmp, HASH_DELETE);
}

/* inline */ void 	
Registry::RemoveType (const char * n)
{
    struct Element tmp;
    tmp.key = (char *) n;
    HASHsearch (active_types, &tmp, HASH_DELETE);
}

/* inline */ STEPentity *
Registry::ObjCreate (const char * nm, int check_case) const
{	
    const EntityDescriptor *  entd = FindEntity (nm, check_case);
    if (entd) return ((EntityDescriptor *)entd) -> NewSTEPentity ();
    else return ENTITY_NULL;
}    

/* inline */ int
Registry::GetEntityCnt ()  
{
    return entity_cnt;
}

/* inline */ void
Registry::ResetEntities ()  
{
    HASHlistinit (primordialSwamp, &cur_entity);
    
}

/* inline */ const EntityDescriptor *
Registry::NextEntity () 
{
  if (0 == HASHlist (&cur_entity)) return 0;
  return (const EntityDescriptor *) cur_entity.e->data;
}    

/* inline */ void
Registry::ResetSchemas ()  
{
    HASHlistinit (active_schemas, &cur_schema);
}

/* inline */ const SchemaDescriptor *
Registry::NextSchema () 
{
  if (0 == HASHlist (&cur_schema)) return 0;
  return (const SchemaDescriptor *) cur_schema.e->data;
}
