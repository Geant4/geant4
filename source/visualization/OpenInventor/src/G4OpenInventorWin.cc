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
// $Id: G4OpenInventorWin.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
// OpenInventor graphics system factory.

#ifdef G4VIS_BUILD_OIWIN32_DRIVER

// this :
#include "G4OpenInventorWin.hh"

#include <Inventor/Win/SoWin.h>

#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorWinViewer.hh"
#include "G4Win32.hh"

G4OpenInventorWin::G4OpenInventorWin ()
:G4OpenInventor("OpenInventorWin","OIWin32",G4VGraphicsSystem::threeD)
,fInited(false)
{
}

void G4OpenInventorWin::Initialize ()
{
  if(fInited) return; //Done

  SetInteractorManager (G4Win32::getInstance());
  GetInteractorManager () -> RemoveDispatcher((G4DispatchFunction)G4Win32::DispatchWin32Event);
  GetInteractorManager () -> AddDispatcher((G4DispatchFunction)SoWin::dispatchEvent);

  HWND toplevel = (HWND)GetInteractorManager()->GetMainInteractor();

  if(!SoWin::getTopLevelWidget()) SoWin::init(toplevel);

  InitNodes();

  fInited = true;
}

G4OpenInventorWin::~G4OpenInventorWin () {}
G4VViewer* G4OpenInventorWin::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
  Initialize();
  G4OpenInventorSceneHandler* pScene = (G4OpenInventorSceneHandler*)&scene;
  return new G4OpenInventorWinViewer (*pScene, name);
}


#endif
