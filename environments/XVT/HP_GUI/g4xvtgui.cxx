// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: g4xvtgui.cxx,v 1.2 1999-11-11 16:01:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*********************************************************************
 *
    Copyright (c) 1996

    File Name:    G4XvtGUI.cxx

    Implements the application class
 *
 *********************************************************************/

#include "XVTPwr.h"
#include "AppDef.h"

#include GUI_i
#include G4XvtGUI_i
#include CGlobalClassLib_i
#include CTaskWin_i
#include CDocument_i
#include CMenuBar_i
#include "g4rw/ordcltn.h"

#ifdef DSP_RELEASE
#include CTaskWin_f
#endif

// Define the Run-Time Type Identification for this class
// PWRRegisterClass1(G4XvtGUI, G4XvtGUIID, "G4XvtGUI", CApplication)

// --- Include attachment files --- //
#include CToolBarAttachment_i
#include CAttachmentFrame_i
// #include CUserWin_f
#include G4XvtWin_f


///////////////////////////////////////////////////////////////////////////////
//
// Application class constructor
//
// Defines the initial state of the application object.
//
///////////////////////////////////////////////////////////////////////////////
G4XvtGUI::G4XvtGUI( CApplication::ApplicationByteAware theAppByte ) : CApplication( theAppByte )
{
	// Set any System Attribute Table values that
	// need to be set before GUI initialization

#ifdef DSP_RELEASE
	GUIFactory.ConstructApplication( this, G4XvtGUI1002 );
#endif

	// For example...
#if 0
#if (XVTWS == WINWS || XVTWS == NTWS )

	// Follow the MDI look and feel when running on MS-Windows
	xvt_vobj_set_attr( NULL_WIN, ATTR_WIN_MDI, TRUE );

#endif

#if (XVTWS == WINWS || XVTWS == NTWS || XVTWS == PMWS)

	// Allow drawing in the task window:
	xvt_vobj_set_attr(NULL_WIN, ATTR_WIN_PM_DRAWABLE_TWIN, TRUE);

#endif

#if (XVTWS == MTFWS)

	// When running on Motif, don't display the "ghost" window
	xvt_vobj_set_attr( NULL_WIN, ATTR_X_DISPLAY_TASK_WIN, FALSE );

#endif
#endif
}

///////////////////////////////////////////////////////////////////////////////
//
// Application class destructor
//
// Delete the internal data objects for this class
//
///////////////////////////////////////////////////////////////////////////////
G4XvtGUI::~G4XvtGUI( void )
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Startup
//
// Defines the startup behavior for the application - whether an
// initial about box is created, or a login screen is presented,
// or an initial document object and its associated windows are
// created, etc.  Also, Initialize application level objects here.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtGUI::StartUp( void )
{
	// NOTE: The base class method must be called before adding any
	// of your own code.
	CApplication::StartUp();

#ifdef DSP_RELEASE
	GUIFactory.InitializeApplication( this, G4XvtGUI1002 );
#endif

#if 0
#if (XVTWS == WINWS || XVTWS == NTWS || XVTWS == PMWS)
	// To instantiate the views sketched out in the TaskWindow, use
	// the following code:
	CTaskWin* aTaskWindow = G->GetTaskWin();
	GUIFactory.DoCreateViews( aTaskWindow, CTaskWin1004, TRUE );
#endif
#endif

	// Create an instance of its documents and display their initial windows
	DoNew();
}

///////////////////////////////////////////////////////////////////////////////
//
// ShutDown
//
// Defines the termination behavior for the application - typically
// it deletes application level objects that were created in StartUp.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtGUI::ShutDown( void )
{
	// Base class method must be called as the last thing this method does
	CApplication::ShutDown();
}

///////////////////////////////////////////////////////////////////////////////
//
// SetUpMenus
//
// Initialize the application menubar - the one that is
// available when there are no windows with menubars available.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtGUI::SetUpMenus( CMenuBar* /*theMenuBar*/ )
{
	// Example code to Initialize default menus
	// theMenuBar->SetEnabled( M_FILE_NEW, TRUE );
	// theMenuBar->SetEnabled( M_FILE_OPEN, TRUE );
	// theMenuBar->SetEnabled( M_EDIT, FALSE );
	// theMenuBar->SetEnabled( M_FONT_POPUP, FALSE );
}

///////////////////////////////////////////////////////////////////////////////
//
// DoCommand
//
// Handle commands here that have to do with behavior specific
// to the application, or that create, destroy or otherwise
// manage the data (document) objects of the application.
//
///////////////////////////////////////////////////////////////////////////////
void G4XvtGUI::DoCommand( long theCommand, void* theData )
{
	switch ( theCommand )
	{
	case NULLcmd:    // stop propagation of "undefined" commands
		break;

#if 0
#if (XVTWS == MTFWS)
	// if the application is running on motif, Terminate the app when the
	// last window is closed to Follow motif style guidelines
	case NOWINDOWcmd:
		DoClose();          // clean up TaskDoc
		xvt_app_destroy();  // Terminate the application
		break;
#endif
#endif
	default:
		CApplication::DoCommand( theCommand, theData );
		break;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// DoNew
//
// This method is automatically called when the New option on the
// default File menu is selected.  Create all document objects
// that had the AutoCreate option specified in XVT-Architect.
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtGUI::DoNew( void )
{
	RWOrdered anExistingDocList( *itsDocuments );
	GUIFactory.DoCreateDocuments( this, G4XvtGUI1002, FALSE, &itsData );
	RWOrderedIterator nextDoc( *itsDocuments );
	CDocument* aDoc;
	RWOrdered orphanedDocs;
	BOOLEAN allDocsCreated = TRUE;
	while ( (aDoc = (CDocument*)nextDoc()) != NULL )
	if ( !anExistingDocList.find( aDoc ) && !aDoc->DoNew() )
	{
		orphanedDocs.append( aDoc );
		allDocsCreated = FALSE;
	}
	orphanedDocs.clearAndDestroy();
	return allDocsCreated;
}

///////////////////////////////////////////////////////////////////////////////
//
// DoOpen
//
// This method is automatically called when the Open option on the
// default File menu is selected.  Create a new  document object,
// allow the user to select a file, and Initialize the document
// object from the file.
//
///////////////////////////////////////////////////////////////////////////////
BOOLEAN G4XvtGUI::DoOpen( void )
{
	RWOrdered anExistingDocList( *itsDocuments );
	GUIFactory.DoCreateDocuments( this, G4XvtGUI1002, FALSE, &itsData );
	RWOrderedIterator nextDoc( *itsDocuments );
	CDocument* aDoc;
	RWOrdered orphanedDocs;
	BOOLEAN allDocsCreated = TRUE;
	while ( (aDoc = (CDocument*)nextDoc()) != NULL )
	if ( !anExistingDocList.find( aDoc ) && !aDoc->DoOpen() )
	{
		orphanedDocs.append( aDoc );
		allDocsCreated = FALSE;
	}
	orphanedDocs.clearAndDestroy();
	return allDocsCreated;
}


// --- New code for the attachment stuff --- //

///////////////////////////////////////////////////////////////////////////////
//
// GetAttachmentFrame
//
// Returns a pointer to the attacment frame created for the application's 
// task window. A return value of NULL indicates that no attachment frame
// was created.
//
///////////////////////////////////////////////////////////////////////////////
CAttachmentFrame* G4XvtGUI::GetAttachmentFrame() const 
{ 
        return itsAttachmentFrame; 
}

