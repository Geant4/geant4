// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: editordefines.h,v 1.1 1999-01-07 16:08:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef editordefines_h
#define editordefines_h

/*
* NIST STEP Editor Class Library
* cleditor/editordefines.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

class GenericNode;
class GenNodeList;

class MgrNode;
class MgrNodeList;

class DisplayNode;
class DisplayNodeList;

//////////////////////////////////////////////////////////////////////////////

enum displayStateEnum {
		 mappedWrite,		// has a writable SEE on the screen
		 mappedView,		// has a view only SEE on the screen
		 notMapped,
		 noMapState
		}  ;

//////////////////////////////////////////////////////////////////////////////

enum stateEnum {
		 noStateSE,	// state undefined, not on a List
		 completeSE,	// on saved complete List
		 incompleteSE,	// on saved incomplete List
		 deleteSE,	// on delete List
		 newSE		// on newly created List
		}  ;

/*
   these variable are used by the STEPfile for reading and writing
   files in working session format.
   None of these variable should be set to 'E' as it will disrupt
   the way the read function finds the "ENDSEC;" token.
*/
static const char wsSaveComplete = 'C';
static const char wsSaveIncomplete = 'I';
static const char wsDelete = 'D';
static const char wsNew = 'N';

#endif
