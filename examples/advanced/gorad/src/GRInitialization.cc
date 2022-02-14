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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRInitialization.hh
//   Defines the initialization procedure of Gorad and Geant4
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRInitialization.hh"

#include "GRDetectorConstruction.hh"
#include "GRPhysicsList.hh"
#include "GRActionInitialization.hh"
#include "GRPrimGenActionMessenger.hh"
#include "G4GenericMessenger.hh"
#include "GRGeomBiasMessenger.hh"
#include "GRScoreWriter.hh"

#include "G4RunManager.hh"
#include "G4ScoringManager.hh"
#include "G4UIdirectory.hh"
#include "G4UnitsTable.hh"

#include <CLHEP/Units/SystemOfUnits.h>

GRInitialization::GRInitialization(G4int verboseLvl)
: verboseLevel(verboseLvl)
{
  // adding unites
  new G4UnitDefinition("milligray","mGy","Dose",1.e-3*CLHEP::gray);
  new G4UnitDefinition("microgray","muGy","Dose",1.e-6*CLHEP::gray);
  new G4UnitDefinition("nanogray","nGy","Dose",1.e-9*CLHEP::gray);

  G4ScoringManager::GetScoringManager()->SetScoreWriter(new GRScoreWriter());

  messenger = new G4GenericMessenger(this,"/gorad/","GORAD commands");
  auto& initCmd = messenger->DeclareMethod("initialize",
         &GRInitialization::Initialize,"Initialize Gorad and G4RunManager");
  initCmd.SetToBeBroadcasted(false);
  initCmd.SetStates(G4State_PreInit);

  detector = new GRDetectorConstruction();
  physics = new GRPhysicsList();
  actionInitialization = new GRActionInitialization();
  sourceMessenger = new GRPrimGenActionMessenger();
  geomBiasMessenger = new GRGeomBiasMessenger(detector,physics,verboseLvl);
}

GRInitialization::~GRInitialization()
{
  delete geomBiasMessenger;
  delete sourceMessenger;
  delete messenger; 
}

void GRInitialization::Initialize()
{
  auto runManager = G4RunManager::GetRunManager();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(actionInitialization);
  if(verboseLevel>0) G4cout << "GORAD is initialized.........." << G4endl;
  runManager->Initialize();
  sourceMessenger->UpdateParticleList();
}

#include "G4UIExecutive.hh"
#ifdef G4UI_USE_QT
#include "G4UIQt.hh"
#endif

void GRInitialization::SetWindowText(G4UIExecutive* ui)
{
  // If the current GUI is not G4UIQt, do nothing and return.
  if(!(ui->IsGUI())) return;

#ifdef G4UI_USE_QT
  G4UIQt* qt = dynamic_cast<G4UIQt*>(ui->GetSession());
  if(!qt) return;

  qt->SetStartPage(std::string("<table width='100%'><tr><td width='50%'></td><td><div ")+
       "style='color: rgb(140, 31, 31); font-size: xx-large; font-family: Garamond, serif; "+
       "padding-bottom: 0px; font-weight: normal'>GORAD "+
       "</div></td></td></tr></table>"+
       "<p>&nbsp;</p>"+
       "<div><dl>"+
         "<dd><b>Gorad (Geant4 Open-source Radiation Analysis and Design) is meant to be "+
         "a turn-key application for radiation analysis and spacecraft design "+
         "built on top of Geant4. Simulation geometry should be provided in the form of GDML. "+
         "Gorad is developed under the NASA JSC contract NNJ15HK11B."+
         "</dd></dl></div>"+
       "<p>&nbsp;</p>"+
       "<div style='background:#EEEEEE;'><b>Tooltips :</b><ul>"+
         "<li><b>Start an interactive run :</b><br />"+
         "/control/execute <i>run.mac</i><br />"+
         "/run/beamOn <i>number_of_events</i></li></ul></div>"+
       "<div style='background:#EEEEEE;'><b>Documentation :</b><ul>"+
         "<li><i>"+
         "<b>GORAD manual</b> and a sample Orion spacecraft shield geometry can be found at<br />"+
         "<a href='https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesGorad'>"+
         "https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesGorad"+
         "</a></i></li>"+
         "</ul></div>"
       );
#endif
}

