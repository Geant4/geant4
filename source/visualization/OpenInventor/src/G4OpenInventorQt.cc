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
// Frederick Jones TRIUMF 30 NOV 2017

// this :
#include "G4OpenInventorQt.hh"

#include "G4SoQt.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorQtViewer.hh"

#ifndef G4GMAKE
#include "moc_G4OpenInventorQt.cpp"
#endif

#define G4warn G4cout

G4OpenInventorQt::G4OpenInventorQt()
  : G4OpenInventor("OpenInventorQt", "OIQt", G4VGraphicsSystem::threeD),
    fInited(false)
{
}

void G4OpenInventorQt::Initialize()
{
  if(fInited) return; //Done

  // FWJ DEBUG
  //  G4cout << "G4OpenInventorQt: SETINTERACTORMANAGER " << G4SoQt::getInstance()
  //         << G4endl;

  SetInteractorManager(G4SoQt::getInstance());

  // FWJ DEBUG
  // G4cout << "G4OpenInventorQt: GETINTERACTORMANAGER " << 
  //    GetInteractorManager() << G4endl;

  // FWJ for now, create an independent main window
  // NOW should be done in G4SoQt [public G4VInteractorManager]
  // QWidget* mainWin = SoQt::init("Geant4");

  // In parent G4OpenInventor
  InitNodes();

  fInited = true;
}

G4OpenInventorQt::~G4OpenInventorQt()
{
}

G4VViewer* G4OpenInventorQt::CreateViewer(G4VSceneHandler& scene,
                                          const G4String& name)
{
  auto pView = new G4OpenInventorQtViewer(static_cast<G4OpenInventorSceneHandler&>(scene), name);

  if (pView) {
    if (pView->GetViewId() < 0) {
      G4warn <<
      "G4OpenInventorQt::CreateViewer: ERROR flagged by negative"
      " view id in G4OpenInventorQtViewer creation."
      "\n Destroying view and returning null pointer."
      << G4endl;
      delete pView;
      pView = 0;
    }
  }
  if (!pView) {
    G4warn <<
    "G4OpenInventorQt::CreateViewer: ERROR: null pointer on new G4OpenInventorQtViewer."
    << G4endl;
  }

  Initialize();
  
  return pView;
}
