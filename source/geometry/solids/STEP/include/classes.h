

//



//
// $Id: classes.h,v 1.2 1999-05-21 20:20:37 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*
** Fed-x parser output module for generating C++  class definitions
** December  5, 1989
** release 2 17-Feb-1992
** release 3 March 1993
** release 4 December 1993
** K. C. Morris
**
** Development of Fed-x was funded by the United States Government,
** and is not subject to copyright.

*******************************************************************
The conventions used in this binding Follow the proposed specification
for the STEP Standard Data Access Interface as defined in document
N350 ( August 31, 1993 ) of ISO 10303 TC184/SC4/WG7.
*******************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>


#include "express.h"
#include "exppp.h"
#include <dict.h>

#ifdef __CENTERLINE__
#define CONST
#else
#define CONST const
#endif

#define MAX_LEN		240
#define DEBUG		if (0) printf 

#define TD_PREFIX	"t_"
#define ATTR_PREFIX	"a_"
#define ENT_PREFIX	"e_"
#define SCHEMA_PREFIX	"s_"

#define TYPEprefix(t)	     (TYPEis_entity (t) ? ENT_PREFIX : TD_PREFIX)

#define SCHEMA_FILE_PREFIX	"Sdai"
#define TYPE_PREFIX   "Sdai"
#define ENTITYCLASS_PREFIX	TYPE_PREFIX
#define ENUM_PREFIX	"sdai"

#define move(b)		(b = (b + strlen(b)))
#define TYPEtd_name(t,s)	TypeDescriptorName ((t), (s))

typedef  struct file_holder  {
    FILE*  inc;		/*  include file  */
    FILE*  lib;		/*  library file  */
    FILE*  incall;	/*  include file for collecting all include files  */
    FILE*  initall;	/*  for registering all entities from all schemas  */
    FILE*  make;	/*  for indicating schema source files in makefile  */
    FILE*  code;	
    FILE*  init;	/*  contains function to Initialize program
			    to use schema's entities */
    
}  File_holder, FILES;

typedef struct EntityTag *EntityTag;
struct EntityTag {
  /*  these fields are used so that ENTITY types are processed in order
   *  when appearing in differnt schemas   */
  unsigned int started :1; /*  marks the beginning of processing  */
  unsigned int complete :1;  /*  marks the end of processing  */

  Entity superclass;  /*  the entity being used as the supertype
			*  - with multiple inheritance only chose one */
};

Entity ENTITYget_superclass (Entity entity);
Entity ENTITYput_superclass (Entity entity);
int ENTITYhas_explicit_attributes (Entity e);

typedef struct SelectTag *SelectTag;
struct SelectTag {
  /*  these fields are used so that SELECT types are processed in order  */
  unsigned int started :1; /*  marks the beginning of processing  */
  unsigned int complete :1;  /*  marks the end of processing  */
  };

#define String  const char *
const char *	CheckWord( const char *);
const char * 	StrToLower(const char *);
const char *	StrToUpper (const char *);
const char *	FirstToUpper (const char *);
const char *	SelectName (const char *); 
FILE	*FILEcreate(const char *);
void	FILEclose(FILE *);
const char *  ClassName(const char *);
const char *  ENTITYget_classname(Entity);
void	ENTITYPrint(Entity,FILES *,Schema);
const char *	StrToConstant(const char *);
void	TYPEprint_definition(Type, FILES *,Schema);
const char * PrettyTmpName (const char * oldname);
const char * EnumName (const char * oldname)  ;
const char * TypeDescriptorName (Type , Schema );
char * TypeDescription (const Type t);
const char * TypeName (Type t);
const char * AccessType (Type t);
const char * TYPEget_ctype (const Type t);

/*Variable*/
/*VARis_simple_explicit (Variable a)*/
#define VARis_simple_explicit(a)  (!VARis_type_shifter(a))

/*Variable*/
/*VARis_simple_derived (Variable a)*/
#define VARis_simple_derived(a)  (!VARis_overrider(a))

