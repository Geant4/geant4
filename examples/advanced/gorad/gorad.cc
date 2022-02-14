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
// gorad.cc
//   Main function of Gorad
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#include "GRInitialization.hh"

int main(int argc,char** argv)
{
  // Construct Run Manager
  auto runManager = G4RunManagerFactory::CreateRunManager();

  // construct Gorad initializer
  auto goradInitialization = new GRInitialization();

  // Visualization manager construction
  auto visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  if ( argc == 2 ) {
    // execute an argument macro file 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    auto UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // start interactive session
    auto ui = new G4UIExecutive(argc, argv); 
    goradInitialization->SetWindowText(ui);
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  delete goradInitialization;
  delete visManager;
  delete runManager;
}

