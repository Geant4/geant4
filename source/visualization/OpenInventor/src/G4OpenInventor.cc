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
// $Id: G4OpenInventor.cc 91687 2015-07-31 09:44:16Z gcosmo $
//
#ifdef G4VIS_BUILD_OI_DRIVER

// this :
#include "G4OpenInventor.hh"

#include "HEPVis/nodes/SoBox.h"
#include "HEPVis/nodes/SoTubs.h"
#include "HEPVis/nodes/SoCons.h"
#include "HEPVis/nodes/SoTrd.h"
#include "HEPVis/nodes/SoTrap.h"
#include "HEPVis/nodes/SoMarkerSet.h"
#include "HEPVis/nodes/SoImageWriter.h"
#include "HEPVis/nodekits/SoDetectorTreeKit.h"
#include "HEPVis/actions/SoGL2PSAction.h"
#include "HEPVis/actions/SoCounterAction.h"
#include "HEPVis/actions/SoAlternateRepAction.h"

#include "Geant4_SoPolyhedron.h"

#include "G4OpenInventorSceneHandler.hh"

#include "G4UImanager.hh"
#include "G4UIbatch.hh"

G4OpenInventor::G4OpenInventor (
 const G4String name
,const G4String nickname
,G4VGraphicsSystem::Functionality f
)
:G4VGraphicsSystem(name,nickname,f)
,interactorManager(0)
{
}

G4OpenInventor::~G4OpenInventor () {}

void G4OpenInventor::SetInteractorManager (G4VInteractorManager* im) {
  interactorManager = im;
}
G4VInteractorManager* G4OpenInventor::GetInteractorManager () {
  return interactorManager;
}
G4VSceneHandler* G4OpenInventor::CreateSceneHandler (const G4String& name) {
  Initialize();
  G4VSceneHandler* p = new G4OpenInventorSceneHandler (*this, name);
  return    p;
}

void G4OpenInventor::InitNodes()
{
  SoBox::initClass();
  SoTubs::initClass();
  SoCons::initClass();
  SoTrd::initClass();
  SoTrap::initClass();
  SoDetectorTreeKit::initClass();
  HEPVis_SoMarkerSet::initClass();
  SoImageWriter::initClass();
  Geant4_SoPolyhedron::initClass();

  SoGL2PSAction::initClass();
  SoCounterAction::initClass();
  SoAlternateRepAction::initClass();
}

G4bool G4OpenInventor::IsUISessionCompatible () const
{
  G4bool isCompatible = false;
  G4UImanager* ui = G4UImanager::GetUIpointer();
  G4UIsession* session = ui->GetSession();

  // If session is a batch session, it may be:
  // a) this is a batch job (the user has not instantiated any UI session);
  // b) we are currently processing a UI command, in which case the UI
  //    manager creates a temporary batch session and to find out if there is
  //    a genuine UI session that the user has instantiated we must drill
  //    down through previous sessions to a possible non-batch session.
  while (G4UIbatch* batch = dynamic_cast<G4UIbatch*>(session)) {
    session = batch->GetPreviousSession();
  }

  // OI windows are not appropriate in a batch session.
  if (session) {
    // If non-zero, this is the originating non-batch session
    // The user has instantiated a UI session...
    isCompatible = true;
  }

  return isCompatible;
}

#endif
