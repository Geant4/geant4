//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: g4xvtdoc.h,v 1.3 2001-07-11 09:56:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/********************************************************************
 *
    Copyright (c) 1996

    File Name:    G4XvtDoc.h
    Author:
    Date:

    Class:        G4XvtDoc
    Inheritance:  CDocument->G4XvtDoc
    Helper(s):

    Purpose:      This is an example of how to override the minimal set of
                  methods of CDocument. It initializes data maintained by
                  the class, creates a window in which to display the data,
                  implements methods for handling requests to access the
                  data, and provides a basic skeleton for persistence.

    Usage:        Override or add methods as necessary.

    Modifications:
 *
 ********************************************************************/

#ifndef G4XvtDoc_H
#define G4XvtDoc_H

#include "XVTPwr.h"
#include "AppDef.h"

#include GUI_i
#include CDocument_i

#ifdef DSP_RELEASE
#include G4XvtDoc_f
#endif

//#include CTypeInfo_i

class G4XvtDoc : public CDocument
{
public:
	G4XvtDoc( CApplication * theApplication = G -> GetApplication( ),
		PWRID theDocumentId = G -> GetId( ) );
	virtual ~G4XvtDoc( void );

	virtual void BuildWindow( void );
	virtual void DoCommand( long theCommand, void * theData = NULL );
	virtual void DoMenuCommand( MENU_TAG theMenuItem,
		BOOLEAN isShiftKey, BOOLEAN isControlKey );

	virtual BOOLEAN DoSave( void );
	virtual BOOLEAN DoOpen( void );

	virtual BOOLEAN Save( void );
	virtual BOOLEAN Open( void );
        
protected:
	// By default, copy and assignment are disallowed

	G4XvtDoc( const G4XvtDoc& theDocument ) : CDocument( theDocument ) { }
	G4XvtDoc& operator=( const G4XvtDoc& theDocument ) { return( *this ); }

private:
	G4DocData itsData;

        WINDOW theWin;

//	PWRClassInfo
};

#endif // G4XvtDoc_H
