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
// $Id: PhotInVisManager.cc,v 1.4 2006/06/29 16:25:31 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

//#define debug

#ifdef G4VIS_USE

#include "PhotInVisManager.hh"

PhotInVisManager::PhotInVisManager () {}

void PhotInVisManager::RegisterGraphicsSystems()
{
#ifdef debug
  G4cout<<"PhotInVisManager::RegisterGraphicsSystems is called"<<G4endl;
#endif
  // Graphics Systems of not needed external packages or libraries...
  RegisterGraphicsSystem (new G4ASCIITree);
  RegisterGraphicsSystem (new G4DAWNFILE);
  RegisterGraphicsSystem (new G4GAGTree);
  RegisterGraphicsSystem (new G4HepRepFile);
  RegisterGraphicsSystem (new G4HepRep);
  RegisterGraphicsSystem (new G4RayTracer);
  RegisterGraphicsSystem (new G4VRML1File);
  RegisterGraphicsSystem (new G4VRML2File);

  // Graphics systems of needed external packages or libraries...

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

  if (fVerbose > 0)
  {
    G4cout<<"PhotInVisManager::RegisterGraphicsSystems: chosen graphics systems"<<G4endl;
    PrintAvailableGraphicsSystems ();
  }
}

#endif

