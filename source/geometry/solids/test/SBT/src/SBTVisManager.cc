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
// SBTVisManager
//
// Implementation of visualization manager for SBT
//

#include "SBTVisManager.hh"
#include "G4ModelingParameters.hh"
#include "G4BoundingSphereScene.hh"

//
// The following is borrowed from GEANT4's example "N02"
//
// Supported drivers...
 
// Not needing external packages or libraries...
#include "G4ASCIITree.hh"
#include "G4DAWNFILE.hh"
#include "G4HepRepFile.hh"
#include "G4HepRep.hh"
#include "G4RayTracer.hh"
#include "G4VRML1File.hh"
#include "G4VRML2File.hh"

#ifdef G4VIS_USE_DAWN
#include "G4FukuiRenderer.hh"
#endif
 
#ifdef G4VIS_USE_OPENGLX
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLStoredX.hh"
#endif
 
#ifdef G4VIS_USE_OPENGLWIN32
#include "G4OpenGLImmediateWin32.hh"
#include "G4OpenGLStoredWin32.hh"
#endif
 
#ifdef G4VIS_USE_OPENGLXM
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLStoredXm.hh"
#endif
 
#ifdef G4VIS_USE_OIX
#include "G4OpenInventorX.hh"
#endif

#ifdef G4VIS_USE_OIWIN32
#include "G4OpenInventorWin32.hh"
#endif


#ifdef G4VIS_USE_VRML
#include "G4VRML1.hh"
#include "G4VRML2.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBTVisManager::RegisterGraphicsSystems () {

  // Graphics Systems not needing external packages or libraries...
  RegisterGraphicsSystem (new G4ASCIITree);
  RegisterGraphicsSystem (new G4DAWNFILE);
  RegisterGraphicsSystem (new G4HepRepFile);
  RegisterGraphicsSystem (new G4HepRep);
  RegisterGraphicsSystem (new G4RayTracer);
  RegisterGraphicsSystem (new G4VRML1File);
  RegisterGraphicsSystem (new G4VRML2File);

#ifdef G4VIS_USE_DAWN
  RegisterGraphicsSystem (new G4FukuiRenderer);
#endif

#ifdef G4VIS_USE_OPENGLX
  RegisterGraphicsSystem (new G4OpenGLImmediateX);
  RegisterGraphicsSystem (new G4OpenGLStoredX);
#endif
 
#ifdef G4VIS_USE_OPENGLWIN32
  RegisterGraphicsSystem (new G4OpenGLImmediateWin32);
  RegisterGraphicsSystem (new G4OpenGLStoredWin32);
#endif
 
#ifdef G4VIS_USE_OPENGLXM
  RegisterGraphicsSystem (new G4OpenGLImmediateXm);
  RegisterGraphicsSystem (new G4OpenGLStoredXm);
#endif
 
#ifdef G4VIS_USE_OIX
  RegisterGraphicsSystem (new G4OpenInventorX);
#endif
 
#ifdef G4VIS_USE_OIWIN32
  RegisterGraphicsSystem (new G4OpenInventorWin32);
#endif
 
#ifdef G4VIS_USE_VRML
  RegisterGraphicsSystem (new G4VRML1);
  RegisterGraphicsSystem (new G4VRML2);
#endif
 
  if (fVerbose > 0) {
    G4cout <<
      "\nYou have successfully chosen to use the following graphics systems."
         << G4endl;
    PrintAvailableGraphicsSystems ();
  }
}



//
// BuildFakeWorld
//
// It is a little troubling how difficult it is to get visualization
// working if one doesn't build geometry in the standard manner
// (i.e. have a valid setting at 
//    G4TransportationManager::GetTransportationManager ()
//      -> GetNavigatorForTracking () -> GetWorldVolume ();   )
//
//
// Why is this so? Why shouldn't one be able to use visualization
// in cases like SBT?
//
G4int SBTVisManager::BuildFakeWorld() const
{
        //
        // These are probably leaks...
        //
        G4ModelingParameters::DrawingStyle style = G4ModelingParameters::wf;
        G4ModelingParameters *model = new G4ModelingParameters
                                            (0,      // No default vis attributes.
                                             style,  // Wireframe
                                             true,   // Global culling.
                                             true,   // Cull invisible volumes.
                                             false,  // Density culling.
                                             0.,     // Density (not relevant if density culling false).
                                             true,   // Cull daughters of opaque mothers.
                                             24);    // No of sides (not relevant for this operation).
        SBTFakeModel *fakeModel = new SBTFakeModel(model);

        G4Scene *currentScene = GetCurrentScene();

        if (!currentScene) {
                G4cerr << "Please create a view first" << G4endl;
                return 1;
        }
        
        currentScene->AddRunDurationModel( (G4VModel *)fakeModel );
        return 0;
}
