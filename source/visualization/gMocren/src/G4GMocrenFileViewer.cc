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
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//


#define __G_ANSI_C__
#define G4GMocrenFile_STRUCTURE_PRIORITY  1.

#include "G4ios.hh"
#include <cstdio>
#include <cstring>
#include <cassert>

#include "G4VisManager.hh"
#include "G4Scene.hh"
#include "G4Vector3D.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#include "G4GMocrenFile.hh"
#include "G4GMocrenFileSceneHandler.hh"
#include "G4GMocrenFileViewer.hh"
#include "G4GMocrenMessenger.hh"


//----- constants

//-- for a debugging
const G4bool GFDEBUG = false;

//----- G4GMocrenFileViewer, constructor
G4GMocrenFileViewer::G4GMocrenFileViewer (G4GMocrenFileSceneHandler& sceneHandler,
					  G4GMocrenMessenger &,
					  const G4String& name)
  : G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
    kSceneHandler (sceneHandler)
{
  // Set a g4.gdd-file viewer 
  std::strncpy( kG4GddViewer, "gMocren", 8);
  if( std::getenv( "G4GMocrenFile_VIEWER" ) != NULL ) {
    char * env = std::getenv( "G4GMocrenFile_VIEWER" );
    G4int len = (G4int)std::strlen(env);
    if(len >= 32) {
      G4Exception("G4GMocrenFileViewer::G4GMocrenFileViewer(*)",
                  "gMocren1000", FatalException,
                  "Invalid length of string set in G4GMocrenFile_VIEWER");
    }
    std::strncpy( kG4GddViewer, env, sizeof(kG4GddViewer) - 1);
    kG4GddViewer[sizeof(kG4GddViewer) - 1] = '\0';
    //std::strcpy( kG4GddViewer, getenv( "G4GMocrenFile_VIEWER" ) ) ;				
  } 

  // string for viewer invocation
  if ( !std::strcmp( kG4GddViewer, "NONE" ) ) {
		
    //std::strcpy( kG4GddViewerInvocation, "" );
    kG4GddViewerInvocation[0] = '\0';
  } else {

    std::strncpy( kG4GddViewerInvocation, kG4GddViewer,
                 sizeof(kG4GddViewerInvocation) - 1);
    kG4GddViewerInvocation[sizeof(kG4GddViewerInvocation) - 1] = '\0';
    G4int n = sizeof(kG4GddViewerInvocation)
            - (G4int)std::strlen(kG4GddViewerInvocation) - 1;
    std::strncat( kG4GddViewerInvocation, " ", n);
    const char * gddfname = kSceneHandler.GetGddFileName();
    G4int len = (G4int)std::strlen(gddfname);
    if(len >= 64) {
      G4Exception("G4GMocrenFileViewer::G4GMocrenFileViewer(*)",
                  "gMocren1001", FatalException,
                  "Invalid length of the GDD file name");
    }
    n = sizeof(kG4GddViewerInvocation)
      - (G4int)std::strlen(kG4GddViewerInvocation) - 1;
    std::strncat( kG4GddViewerInvocation, gddfname, n);
  }

}

//----- G4GMocrenFileViewer, destructor
G4GMocrenFileViewer::~G4GMocrenFileViewer () 
{}

//----- G4GMocrenFileViewer::SetView () 
void G4GMocrenFileViewer::SetView () 
{
  if(GFDEBUG)
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "***** G4GMocrenFileViewer::SetView(): No effects" << G4endl;

  // Do nothing, since DAWN is running as a different process.
  // SendViewParameters () will do this job instead.
}


//----- G4GMocrenFileViewer::ClearView()
void
G4GMocrenFileViewer::ClearView( void )
{
  if(GFDEBUG) {
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
      G4cout << "***** G4GMocrenFileViewer::ClearView (): No effects " << G4endl;
    }
  }  
  //if(kSceneHandler.kGddDest) {
    //kSceneHandler.kGddDest.close();
    // Re-open with same filename...
    //kSceneHandler.kGddDest.open(kSceneHandler.kGddFileName);
    kSceneHandler.kFlagInModeling = false;
    kSceneHandler.GFBeginModeling();
    //}
}


//----- G4GMocrenFileViewer::DrawView () 
void G4GMocrenFileViewer::DrawView () 
{
  if(GFDEBUG)
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "***** G4GMocrenFileViewer::DrawView () " << G4endl;

  //----- 
  kSceneHandler.GFBeginModeling() ;

  //----- Always visit G4 kernel 
  NeedKernelVisit ();
	                           
  //----- Draw
  G4VViewer::ProcessView () ;

} // G4GMocrenFileViewer::DrawView () 



//----- G4GMocrenFileViewer::ShowView()
void G4GMocrenFileViewer::ShowView( void )
{
  if(GFDEBUG)
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "***** G4GMocrenFileViewer::ShowView () " << G4endl;

  if( kSceneHandler.GFIsInModeling() ) 
    {
      //----- End of modeling
      // !EndModeling, !DrawAll, !CloseDevice,
      // close g4.gdd
      kSceneHandler.GFEndModeling();

      //----- Output DAWN GUI file 
      //SendViewParameters(); 

      //----- string for viewer invocation
      if ( !strcmp( kG4GddViewer, "NONE" ) ) {
		
	kG4GddViewerInvocation[0] = '\0';
	//std::strcpy( kG4GddViewerInvocation, "" );
      } else {

        std::strncpy( kG4GddViewerInvocation, kG4GddViewer,
                     sizeof(kG4GddViewerInvocation) - 1);
        kG4GddViewerInvocation[sizeof(kG4GddViewerInvocation) - 1] = '\0';
        G4int n = sizeof(kG4GddViewerInvocation)
                - (G4int)std::strlen(kG4GddViewerInvocation) - 1;
        std::strncat( kG4GddViewerInvocation, " ", n);
        const char * gddfname = kSceneHandler.GetGddFileName();
        G4int len = (G4int)std::strlen(gddfname);
        if(len >= 64) {
          G4Exception("G4GMocrenFileViewer::ShowView()",
                      "gMocren1002", FatalException,
                      "Invalid length of the GDD file name");
        }
        n = sizeof(kG4GddViewerInvocation)
          - (G4int)std::strlen(kG4GddViewerInvocation) - 1;
        std::strncat( kG4GddViewerInvocation, gddfname, n);
      }

    }

} // G4GMocrenFileViewer::ShowView()

