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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det) :
G4UImessenger(), fpDetector(Det)
{
  fpTestDir = std::make_unique<G4UIdirectory>("/microyz/");
  fpTestDir->SetGuidance(" detector control.");
  
  fpDetDir = std::make_unique<G4UIdirectory>("/microyz/det/");
  fpDetDir->SetGuidance("detector construction commands");
      
  fpTrackingCutCmd = std::make_unique<G4UIcmdWithADoubleAndUnit>("/microyz/det/setTrackingCut",this);
  fpTrackingCutCmd->SetGuidance("Set tracking cut");
  fpTrackingCutCmd->SetParameterName("Cut",false);
  fpTrackingCutCmd->SetRange("Cut>0.");
  fpTrackingCutCmd->SetUnitCategory("Energy");
  fpTrackingCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fpTrackingCutCmd->SetToBeBroadcasted(false);

  fpMaxStepSizeCmd = std::make_unique<G4UIcmdWithADoubleAndUnit>("/microyz/det/setMaxStepSize",this);
  fpMaxStepSizeCmd->SetGuidance("Set maximum step size");
  fpMaxStepSizeCmd->SetParameterName("Size",false);
  fpMaxStepSizeCmd->SetRange("Size>0.");
  fpMaxStepSizeCmd->SetUnitCategory("Length");
  fpMaxStepSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fpMaxStepSizeCmd->SetToBeBroadcasted(false);

  fAddRadius = std::make_unique<G4UIcmdWithADoubleAndUnit>("/microyz/det/Radius",this);
  fpMaxStepSizeCmd->SetGuidance("Set TrackSD radius");
  fAddRadius->SetToBeBroadcasted(false);
  fAddRadius->SetParameterName("Radius",false);
  fAddRadius->SetRange("Radius>0.");
  fAddRadius->SetUnitCategory("Length");
  fAddRadius->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fpTrackingCutCmd.get() )
   { fpDetector->SetTrackingCut(fpTrackingCutCmd->GetNewDoubleValue(newValue));}

  if( command == fpMaxStepSizeCmd.get() )
   { fpDetector->SetMaxStepSize(fpMaxStepSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fAddRadius.get() )
  { fpDetector->SetTrackerSDRadius(fAddRadius->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
