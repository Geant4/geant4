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
 
#ifdef G4VIS_USE_DAWN
#include "G4FukuiRenderer.hh"
#endif
 
#ifdef G4VIS_USE_DAWNFILE
#include "G4DAWNFILE.hh"
#endif
 
#ifdef G4VIS_USE_OPACS
#include "G4Wo.hh"
#include "G4Xo.hh"
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

#ifdef G4VIS_USE_RAYX
#include "G4RayX.hh"
#endif

#ifdef G4VIS_USE_VRML
#include "G4VRML1.hh"
#include "G4VRML2.hh"
#endif

#ifdef G4VIS_USE_VRMLFILE
#include "G4VRML1File.hh"
#include "G4VRML2File.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBTVisManager::RegisterGraphicsSystems () {

#ifdef G4VIS_USE_DAWN
  RegisterGraphicsSystem (new G4FukuiRenderer);
#endif

#ifdef G4VIS_USE_DAWNFILE
  RegisterGraphicsSystem (new G4DAWNFILE);
#endif

#ifdef G4VIS_USE_OPACS
  RegisterGraphicsSystem (new G4Wo);
  RegisterGraphicsSystem (new G4Xo);
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
 
#ifdef G4VIS_USE_RAYX
  RegisterGraphicsSystem (new G4RayX);
#endif
 
#ifdef G4VIS_USE_VRML
  RegisterGraphicsSystem (new G4VRML1);
  RegisterGraphicsSystem (new G4VRML2);
#endif
 
#ifdef G4VIS_USE_VRMLFILE
  RegisterGraphicsSystem (new G4VRML1File);
  RegisterGraphicsSystem (new G4VRML2File);
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
	G4ModelingParameters *model = new G4ModelingParameters
					    (0,      // No default vis attributes.
					     G4ModelingParameters::wireframe,
					     true,   // Global culling.
					     true,   // Cull invisible volumes.
					     false,  // Density culling.
					     0.,     // Density (not relevant if density culling false).
					     true,   // Cull daughters of opaque mothers.
					     24,     // No of sides (not relevant for this operation).
					     true,   // View geometry.
					     false,  // View hits - not relevant for physical volume model.
					     false); // View digis - not relevant for physical volume model.
	SBTFakeModel *fakeModel = new SBTFakeModel(model);
	
	if (!fpScene) {
		G4cerr << "Please create a view first" << G4endl;
		return 1;
	}
	
	fpScene->AddRunDurationModel( (G4VModel *)fakeModel );
	return 0;
}
