//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4GMocrenFileViewer.cc,v 1.1 2009-04-01 13:16:11 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Akinori Kimura    March 31, 2009
//


#define __G_ANSI_C__
#define G4GMocrenFile_STRUCTURE_PRIORITY  1.

// #define DEBUG_FR_VIEW

#include "G4ios.hh"
#include <cstdio>
#include <cstring>
#include <cassert>

#include "G4Scene.hh"
#include "G4Vector3D.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#include "G4FRConst.hh"
#include "G4GMocrenFile.hh"
#include "G4GMocrenFileSceneHandler.hh"
#include "G4GMocrenFileViewer.hh"
#include "G4GMocrenMessenger.hh"


//----- constants
const char  FR_ENV_MULTI_WINDOW [] = "G4DAWN_MULTI_WINDOW" ;
const char  FR_ENV_MULTI_WINDOW2[] = "G4GMocrenFile_MULTI_WINDOW" ;

//----- G4GMocrenFileViewer, constructor
G4GMocrenFileViewer::G4GMocrenFileViewer (G4GMocrenFileSceneHandler& sceneHandler,
					  G4GMocrenMessenger & messenger,
					  const G4String& name)
  : G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
    fSceneHandler (sceneHandler),
    fMessenger(messenger)
{
  // Set a g4.gdd-file viewer 
  std::strcpy( fG4GddViewer, "gMocren" ); 
  if( getenv( "G4GMocrenFile_VIEWER" ) != NULL ) {
    std::strcpy( fG4GddViewer, getenv( "G4GMocrenFile_VIEWER" ) ) ;			
  } 

  // string for viewer invocation
  if ( !std::strcmp( fG4GddViewer, "NONE" ) ) {
		
    std::strcpy( fG4GddViewerInvocation, "" );
  } else {

    std::strcpy( fG4GddViewerInvocation, fG4GddViewer );
    std::strcat( fG4GddViewerInvocation, " ");
    std::strcat( fG4GddViewerInvocation, fSceneHandler.GetGddFileName() );
  }

}

//----- G4GMocrenFileViewer, destructor
G4GMocrenFileViewer::~G4GMocrenFileViewer () 
{}

//----- G4GMocrenFileViewer::SetView () 
void G4GMocrenFileViewer::SetView () 
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4GMocrenFileViewer::SetView(): No effects" << G4endl;
#endif 
  // Do nothing, since DAWN is running as a different process.
  // SendViewParameters () will do this job instead.
}


//----- G4GMocrenFileViewer::ClearView()
void
G4GMocrenFileViewer::ClearView( void )
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4GMocrenFileViewer::ClearView (): No effects " << G4endl;
#endif
  if (fSceneHandler.fGddDest.IsOpen()) {
    fSceneHandler.fGddDest.Close();
    // Re-open with same filename...
    fSceneHandler.fGddDest.Open(fSceneHandler.fGddFileName);
    //fSceneHandler.SendStr( FR_G4_GDD_HEADER );
    fSceneHandler.FRflag_in_modeling = false;
    fSceneHandler.FRBeginModeling();
  }
}


//----- G4GMocrenFileViewer::DrawView () 
void G4GMocrenFileViewer::DrawView () 
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4GMocrenFileViewer::DrawView () " << G4endl;
#endif
  //----- 
  fSceneHandler.FRBeginModeling() ;

  //----- Always visit G4 kernel 
  NeedKernelVisit ();
	                           
  //----- Draw
  G4VViewer::ProcessView () ;

} // G4GMocrenFileViewer::DrawView () 



//----- G4GMocrenFileViewer::ShowView()
void G4GMocrenFileViewer::ShowView( void )
{
#if defined DEBUG_FR_VIEW
  G4cerr << "***** G4GMocrenFileViewer::ShowView () " << G4endl;
#endif

  if( fSceneHandler.FRIsInModeling() ) 
    {
      //----- End of modeling
      // !EndModeling, !DrawAll, !CloseDevice,
      // close g4.gdd
      fSceneHandler.FREndModeling();

      //----- Output DAWN GUI file 
      //SendViewParameters(); 

      //----- string for viewer invocation
      if ( !strcmp( fG4GddViewer, "NONE" ) ) {
		
	std::strcpy( fG4GddViewerInvocation, "" );
      } else {

	std::strcpy( fG4GddViewerInvocation, fG4GddViewer );
	std::strcat( fG4GddViewerInvocation, " ");
	std::strcat( fG4GddViewerInvocation, fSceneHandler.GetGddFileName() );
      }


      //----- Invoke DAWN
      /*
	G4cout << G4endl ;
	if( false == G4FRofstream::DoesFileExist( fSceneHandler.GetG4PrimFileName() ) )   
	{
	G4cout << "ERROR: Failed to generate file  ";
	G4cout << fSceneHandler.GetG4PrimFileName() << G4endl;

	} else 	if( strcmp( GetG4PrimViewerInvocation(), "" ) )  
	{
	G4cout << "File  " << fSceneHandler.GetG4PrimFileName() ;
	G4cout << "  is generated." << G4endl;
	G4cout << GetG4PrimViewerInvocation() << G4endl;
	system( GetG4PrimViewerInvocation() );

	} else { // no view, i.e., only file generation
	G4cout << "File  " << fSceneHandler.GetG4PrimFileName() ; 
	G4cout << "  is generated." << G4endl;
	G4cout << "No viewer is invoked." << G4endl;
	}

      */
    }

} // G4GMocrenFileViewer::ShowView()

