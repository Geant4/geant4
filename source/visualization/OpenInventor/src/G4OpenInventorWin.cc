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
// $Id: G4OpenInventorWin.cc,v 1.4 2004/04/08 10:49:57 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
{
  SetInteractorManager (G4Win32::getInstance());
  GetInteractorManager () -> RemoveDispatcher((G4DispatchFunction)G4Win32::dispatchWin32Event);  
  GetInteractorManager () -> AddDispatcher((G4DispatchFunction)SoWin::dispatchEvent);

  HWND toplevel = (HWND)GetInteractorManager()->GetMainInteractor();

  SoWin::init(toplevel);

  InitNodes();
}

G4OpenInventorWin::~G4OpenInventorWin () {}
G4VViewer* G4OpenInventorWin::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
  G4OpenInventorSceneHandler* pScene = (G4OpenInventorSceneHandler*)&scene;
  return new G4OpenInventorWinViewer (*pScene, name);
}


#endif
