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

// ====================================================================
//
//   H02VisManager.cc
//   $Id: H02VisManager.cc,v 1.1 2002-05-28 14:15:48 murakami Exp $
//
// ====================================================================
#include "H02VisManager.hh"

// supported drivers...
#include "G4ASCIITree.hh"

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

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////
H02VisManager::H02VisManager() 
////////////////////////////
{
}

/////////////////////////////
H02VisManager::~H02VisManager() 
/////////////////////////////
{
}

////////////////////////////////////////////
void H02VisManager::RegisterGraphicsSystems() 
////////////////////////////////////////////
{
  RegisterGraphicsSystem (new G4ASCIITree);

#ifdef G4VIS_USE_DAWN
  G4cout << "*  Visualizaton Driver: Fukui Renderer" << G4endl;
  RegisterGraphicsSystem(new G4FukuiRenderer);
#endif

#ifdef G4VIS_USE_DAWNFILE
  G4cout << "*  Visualizaton Driver: DAWN" << G4endl;
  RegisterGraphicsSystem(new G4DAWNFILE);
#endif

#ifdef G4VIS_USE_OPACS
  RegisterGraphicsSystem(new G4Wo);
  RegisterGraphicsSystem(new G4Xo);
#endif

#ifdef G4VIS_USE_OPENGLX
  RegisterGraphicsSystem(new G4OpenGLImmediateX);
  RegisterGraphicsSystem(new G4OpenGLStoredX);
#endif

#ifdef G4VIS_USE_OPENGLWIN32
  RegisterGraphicsSystem(new G4OpenGLImmediateWin32);
  RegisterGraphicsSystem(new G4OpenGLStoredWin32);
#endif

#ifdef G4VIS_USE_OPENGLXM
  RegisterGraphicsSystem(new G4OpenGLImmediateXm);
  RegisterGraphicsSystem(new G4OpenGLStoredXm);
#endif

#ifdef G4VIS_USE_OIX
  RegisterGraphicsSystem(new G4OpenInventorX);
#endif

#ifdef G4VIS_USE_OIWIN32
  RegisterGraphicsSystem(new G4OpenInventorWin32);
#endif

#ifdef G4VIS_USE_RAYX
  RegisterGraphicsSystem(new G4RayX);
#endif

#ifdef G4VIS_USE_VRML
  G4cout << "*  Visualizaton Driver: VRML1/2" << G4endl;
  RegisterGraphicsSystem(new G4VRML1);
  RegisterGraphicsSystem(new G4VRML2);
#endif

#ifdef G4VIS_USE_VRMLFILE
  G4cout << "*  Visualizaton Driver: VRML1/2 File" << G4endl;
  RegisterGraphicsSystem(new G4VRML1File);
  RegisterGraphicsSystem(new G4VRML2File);
#endif

  if(fVerbose> 0) {
    G4cout << G4endl 
	   << "You have successfully chosen to use "
	   << "the following graphics systems."
	   << G4endl;
    PrintAvailableGraphicsSystems ();
  }
}

