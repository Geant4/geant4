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
// $Id: g4xvtgui.h,v 1.3 2001-07-11 09:56:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/********************************************************************
 *
    Copyright (c) 1996

    File Name:    G4XvtGUI.h
    Author:       Simon Prior
    Date:

    Class:        G4XvtGUI
    Inheritance:  CApplication->G4XvtGUI
    Helper(s):    G4XvtDoc, CShellWin

    Purpose:      This is an example of how to override the minimal set of
                  methods of CApplication. It initializes the inherited
                  application, creates a document object to manage the
                  application's data, and implements methods for handling
                  application-level issues.

    Usage:        Override or add methods as necessary.

    Modifications:
 *
 ********************************************************************/

#ifndef G4XvtGUI_H
#define G4XvtGUI_H

#include "XVTPwr.h"
#include "AppDef.h"

#include GUI_i
#include CApplication_i

#ifdef DSP_RELEASE
#include G4XvtGUI_f
#endif

// --- Alert code of the CAttachment frame class use --- //
class CAttachmentFrame;

#include CGlobalClassLib_i
#define ATTACHApp ((CAttachApp*) CObjectRWC::G->GetApplication())

class G4XvtGUI : public CApplication
{
public:

	G4XvtGUI( CApplication::ApplicationByteAware theAppByte = CApplication::AB_SINGLEBYTE );
	virtual ~G4XvtGUI( void );
                
        // --- CAttachApp specific interface --- //
        
        CAttachmentFrame* GetAttachmentFrame() const;

	virtual void StartUp( void );
	virtual void ShutDown( void );
	virtual void SetUpMenus( CMenuBar* theMenuBar );
	virtual void DoCommand( long theCommand, void* theData = NULL );

	virtual BOOLEAN DoNew( void );
	virtual BOOLEAN DoOpen( void );

protected:

	// By default, copy and assignment are disallowed
	G4XvtGUI( const G4XvtGUI& theApplication ) : CApplication( theApplication ) {}
	G4XvtGUI& operator=( const G4XvtGUI& theApplication ) { return *this; }
                
        // --- Class data --- //                

        CAttachmentFrame* itsAttachmentFrame;

private:

	G4XvtGUI1002Data itsData;

//	PWRClassInfo
};

#endif // G4XvtGUI_H
