

//



//
// $Id: sdai.h,v 1.2 1999-05-21 20:20:44 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef SDAI_H
#define	SDAI_H

/*
* NIST STEP Core Class Library
* clstepcore/sdai.h
* May 1995
* David Sauder
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.

* This file is included to support the proposed SDAI/C++ Language Binding
* specified in Annex D of STEP Part 22.

* This file specifies the appropriate naming conventions and NULL
* values for the EXPRESS base types.
*/

/*  */

#ifdef __O3DB__
/*  OpenOODB.h must be the first include file.  */
#include <OpenOODB.h>
#endif

//#include <values.h>
#include <limits.h>

#include <string.h>

//  The values used here to represent NULL values for numeric types come from 
//  the ANSI C standard

//  INTEGER
// may not a good idea because the def SdaiInteger &i may not be legal
//typedef long int SdaiInteger;

#define SdaiInteger long
//#define s_Integer SdaiInteger

// C++ from values.h DAS PORT
//#define S_INT_NULL  MAXLONG 
// C++ from limits.h DAS PORT
#define S_INT_NULL    LONG_MAX 

	//  REAL
//typedef double SdaiReal;
#define SdaiReal double 
//#define	s_Real	SdaiReal

// C++ from values.h DAS PORT
//#  define S_REAL_NULL MINFLOAT
// C++ from limits.h DAS PORT
#define S_REAL_NULL   FLT_MIN

	// NUMBER
	// arbitrary choice by me for number DAS

// C++ from values.h DAS PORT
//#define S_NUMBER_NULL MINFLOAT
// C++ from limits.h DAS PORT
#define S_NUMBER_NULL FLT_MIN

	//  STRING
#define S_STRING_NULL	""
#include <STEPstring.h>

	//  ENTITY
class STEPentity;
extern STEPentity NilSTEPentity;
#define S_ENTITY_NULL	&NilSTEPentity

/******************************************************************************
ENUMERATION
    These types do not use the interface specified.
    The interface used is to enumerated values.  The enumerations are
    mapped in the following way:
    *  enumeration-item-from-schema ==> enumeration item in c++ enum clause
    *  all enumeration items in c++ enum are in upper case 
    *  the value ENUM_NULL is used to represent NULL for all enumerated types 
 *****************************************************************************/

#include <Enumeration.h>

/******************************************************************************
BOOLEAN and LOGICAL

    Logical are implemented as an enumeration in Enumeration.h:
    *  it\'s possible values are sdaiFALSE, sdaiTRUE, sdaiUNKNOWN
    *  or F, T, U

    Boolean are implemented as an enumeration in Enumeration.h:
    *  it\'s possible values are sdaiFALSE, sdaiTRUE
    *  or F, T

******************************************************************************/

typedef Logical  s_Logical;
typedef Logical  SdaiLogical;

typedef Boolean  s_Boolean;
typedef Boolean  SdaiBoolean;

#include <SdaiBinary.h>

/******************************************************************************
AGGREGATE TYPES

    Aggregate types are accessed generically.  (There are not seperate
    classes for the different types of aggregates.)  Aggregates are 
    implemented through the STEPaggregate class.

******************************************************************************/

/******************************************************************************
SELECT

    Selects are represented as subclasses of the SdaiSelectH class in 
    STEPselect.h

******************************************************************************/
#include <STEPselect.h>

#endif
