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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file medical/dna/slowing/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det) :
G4UImessenger(), fpDetector(Det), 
fpTestDir(0), fpDetDir(0), fpMaterCmd(0), fpUpdateCmd(0), fpTrackingCutCmd(0)
{
  fpTestDir = new G4UIdirectory("/slowing/");
  fpTestDir->SetGuidance(" detector control.");
  
  fpDetDir = new G4UIdirectory("/slowing/det/");
  fpDetDir->SetGuidance("detector construction commands");
      
  fpMaterCmd = new G4UIcmdWithAString("/slowing/det/setMat", this);
  fpMaterCmd->SetGuidance("Select material of the world.");
  fpMaterCmd->SetParameterName("choice", false);
  fpMaterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpUpdateCmd = new G4UIcmdWithoutParameter("/slowing/det/update", this);
  fpUpdateCmd->SetGuidance("Update geometry.");
  fpUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fpUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fpUpdateCmd->AvailableForStates(G4State_Idle);

  fpTrackingCutCmd = 
    new G4UIcmdWithADoubleAndUnit("/slowing/det/setTrackingCut",this);
  fpTrackingCutCmd->SetGuidance("Set tracking cut");
  fpTrackingCutCmd->SetParameterName("Cut",false);
  fpTrackingCutCmd->SetRange("Cut>0.");
  fpTrackingCutCmd->SetUnitCategory("Energy");
  fpTrackingCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fpTrackingCutCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fpTestDir;
  delete fpDetDir;
  delete fpMaterCmd;
  delete fpUpdateCmd;
  delete fpTrackingCutCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpMaterCmd) fpDetector->SetMaterial(newValue);
  else if (command == fpUpdateCmd) fpDetector->UpdateGeometry();

  if( command == fpTrackingCutCmd )
   { fpDetector->SetTrackingCut(fpTrackingCutCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
