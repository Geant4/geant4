// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: g4xvtdoc.cxx,v 1.1 1999-01-07 16:05:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*********************************************************************
 *
    Copyright (c) 1996

    File Name:    G4XvtDoc.cxx

    Implements a document class to manage data
 *
 *********************************************************************/

#include "XVTPwr.h"
#include "AppDef.h"

#include GUI_i
#include G4XvtDoc_i

#include G4XvtWin_i

#include NScrollText_i
#include CStatusBar_i

// Define the Run-Time Type Identification for this class
// PWRRegisterClass1(G4XvtDoc, G4XvtDocID, "G4XvtDoc", CDocument)

///////////////////////////////////////////////////////////////////////////////
//
// Document class constructor
//
// Create and Initialize the internal data objects for this class
//
///////////////////////////////////////////////////////////////////////////////
G4XvtDoc::G4XvtDoc( CApplication * theApplication, PWRID theDocumentId )
	: CDocument( theApplication, theDocumentId )
{
	// If the internal data object is derived from CModel, you may also
	// want to instantiate a CController for the data, and take advantage
        // of the automatic data propagation
}

///////////////////////////////////////////////////////////////////////////////
//
// Document class destructor
//
// Delete the internal data objects for this class
//
///////////////////////////////////////////////////////////////////////////////
G4XvtDoc::~G4XvtDoc( void )
{
}

///////////////////////////////////////////////////////////////////////////////
//
// BuildWindow
//
// Create the default view(s) of this data
//
// This method is automatically called by the base class methods
// when new document objects are created or opened.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtDoc::BuildWindow( void )
{
	// Create all window objects associated with this document
	// that have the AutoCreate option specified in XVT-Architect.

	GUIFactory.DoCreateWindows( this, G4Doc, FALSE, &itsData );
}

///////////////////////////////////////////////////////////////////////////////
//
// DoCommand
//
// Handle commands here that have to do with accessing data, or that create,
// destroy or otherwise manage the window objects of the application.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtDoc::DoCommand( long theCommand, void * theData )
{
	switch( theCommand )
	{
	case NULLcmd:
		// stop propagation of "undefined" commands
		break;

	default:
		CDocument::DoCommand( theCommand, theData );
		break;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// DoMenuCommand
//
// Handle menu items here that have to do with accessing data, or that create,
// destroy or otherwise manage the window objects of the application.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtDoc::DoMenuCommand(
	MENU_TAG theMenuItem,
	BOOLEAN  isShiftKey,
	BOOLEAN  isControlKey )
{
	switch( theMenuItem )
	{
	case NULLcmd:
		// stop propagation of "undefined" commands
		break;

	default:
		CDocument::DoMenuCommand( theMenuItem, isShiftKey, isControlKey );
		break;
	}
}
///////////////////////////////////////////////////////////////////////////////
//
// Save
//
// Called by DoSave().  Do the actual work of saving the data.
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtDoc::Save( void )
{
	SetSave( FALSE );
	return( !NeedsSaving( ) );
}

///////////////////////////////////////////////////////////////////////////////
//
// Open
//
// Called by DoOpen().  Do the actual work of restoring the object from
// a data source.
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtDoc::Open( void )
{
	// Add code here to set internal data members of the class based
	// on information read from some data source

	SetSave( FALSE );
	return( !NeedsSaving( ) );
}

///////////////////////////////////////////////////////////////////////////////
//
// DoSave
//
// This method is automatically called when the File Save menu option is
// selected.  Write out the contents of this document object to the path
// specified in the internal FILE_SPEC data member.
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtDoc::DoSave( void )
{
	if( itsXVTFilePointer && itsXVTFilePointer -> name[0] )
		return( Save( ) );
	else
		return( DoSaveAs( ) );
}

///////////////////////////////////////////////////////////////////////////////
//
// DoOpen
//
// This method is automatically called when the File Save menu option is
// selected.  Initialize the contents of this document object from the file
// path in the internal FILE_SPEC data member.  Put up the File Open dialog
// if no file name has been assigned yet.
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtDoc::DoOpen( void )
{
	if( !itsXVTFilePointer )
	{
		itsXVTFilePointer            = new FILE_SPEC;
		itsXVTFilePointer -> type[0] = 0;
		itsXVTFilePointer -> name[0] = 0;
		xvt_fsys_get_default_dir( &itsXVTFilePointer -> dir );
	}

	switch( xvt_dm_post_file_open( itsXVTFilePointer,
		"Enter a new file to open:" ) )
	{
	case FL_OK:
		if( Open( ) )
		{
			BuildWindow( );
			return( TRUE );
		}
		else
		{
			return( FALSE );
		}

	case FL_BAD:
	case FL_CANCEL:
	default:
		return( FALSE );
	}
}
