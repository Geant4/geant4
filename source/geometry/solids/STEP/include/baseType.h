// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: baseType.h,v 1.1 1999-01-07 16:08:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef BASETYPE_H
#define	BASETYPE_H

/*
* NIST STEP Core Class Library
* clstepcore/baseType.h
* May 1995
* David Sauder
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

//	**************  TYPES of attributes

enum PrimitiveType {
	sdaiINTEGER = 0,
	sdaiREAL,
	sdaiBOOLEAN,
	sdaiLOGICAL,
	sdaiSTRING,
	sdaiBINARY,
	sdaiENUMERATION,
	sdaiSELECT,
	sdaiINSTANCE,
	sdaiAGGR,
	sdaiNUMBER,
// The elements defined below are not part of part 23
	ARRAY_TYPE,		// DAS
	BAG_TYPE,		// DAS
	SET_TYPE,		// DAS
	LIST_TYPE,		// DAS
	GENERIC_TYPE,
	REFERENCE_TYPE,
	UNKNOWN_TYPE
 };

// for backwards compatibility with our previous implementation
typedef PrimitiveType BASE_TYPE;

// the previous element types of the enum BASE_TYPE that have been redefined
#define INTEGER_TYPE sdaiINTEGER
#define REAL_TYPE sdaiREAL
#define BOOLEAN_TYPE sdaiBOOLEAN
#define LOGICAL_TYPE sdaiLOGICAL
#define STRING_TYPE sdaiSTRING
#define BINARY_TYPE sdaiBINARY
#define ENUM_TYPE sdaiENUMERATION
#define SELECT_TYPE sdaiSELECT
#define ENTITY_TYPE sdaiINSTANCE
#define AGGREGATE_TYPE sdaiAGGR
#define NUMBER_TYPE sdaiNUMBER

/* not defined in part 23
	ARRAY_TYPE,		// DAS
	BAG_TYPE,		// DAS
	SET_TYPE,		// DAS
	LIST_TYPE,		// DAS
	GENERIC_TYPE,
	REFERENCE_TYPE,
	UNKNOWN_TYPE
*/

#endif 
