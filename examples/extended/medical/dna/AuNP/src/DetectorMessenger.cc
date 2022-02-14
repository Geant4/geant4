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
/// \file medical/dna/AuNP/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 78723 2014-01-20 10:32:17Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),    
 fNPMaterCmd(0),
 fNReplicaRCmd(0),
 fNReplicaAzmCmd(0),
 fAbsRadiusCmd(0),
 fNPRadiusCmd(0),
 fTrackingCutCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/AuNP/");
  fTestemDir->SetGuidance("Detector control.");
  
  fDetDir = new G4UIdirectory("/AuNP/det/");
  fDetDir->SetGuidance("Detector construction commands");
      
  fNPMaterCmd = new G4UIcmdWithAString("/AuNP/det/setNPMat",this);
  fNPMaterCmd->SetGuidance("Select material of the sphere.");
  fNPMaterCmd->SetParameterName("choice",false);
  fNPMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNPMaterCmd->SetToBeBroadcasted(false);  
  
  fNReplicaRCmd = new G4UIcmdWithAnInteger("/AuNP/det/setNReplicaR",this);
  fNReplicaRCmd->SetGuidance("Set Number of Replica in R direction");
  fNReplicaRCmd->SetParameterName("NR",false);
  fNReplicaRCmd->SetRange("NR>0");
  fNReplicaRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNReplicaRCmd->SetToBeBroadcasted(false);
  
  fNReplicaAzmCmd = new G4UIcmdWithAnInteger("/AuNP/det/setNReplicaAzm",this);
  fNReplicaAzmCmd->SetGuidance("Set Number of Replica in Azimuthal direction");
  fNReplicaAzmCmd->SetParameterName("NAzm",false);
  fNReplicaAzmCmd->SetRange("NAzm>0");
  fNReplicaAzmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNReplicaAzmCmd->SetToBeBroadcasted(false);

  fAbsRadiusCmd = new G4UIcmdWithADoubleAndUnit("/AuNP/det/setAbsRadius",this);
  fAbsRadiusCmd->SetGuidance("Set radius of the absorber");
  fAbsRadiusCmd->SetParameterName("Radius",false);
  fAbsRadiusCmd->SetRange("Radius>0.");
  fAbsRadiusCmd->SetUnitCategory("Length");
  fAbsRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsRadiusCmd->SetToBeBroadcasted(false);

  fNPRadiusCmd = new G4UIcmdWithADoubleAndUnit("/AuNP/det/setNPRadius",this);
  fNPRadiusCmd->SetGuidance("Set radius of the nano particle");
  fNPRadiusCmd->SetParameterName("Radius",false);
  fNPRadiusCmd->SetRange("Radius>0.");
  fNPRadiusCmd->SetUnitCategory("Length");
  fNPRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNPRadiusCmd->SetToBeBroadcasted(false);
      
  fTrackingCutCmd = 
    new G4UIcmdWithADoubleAndUnit("/AuNP/det/setTrackingCut",this);
  fTrackingCutCmd->SetGuidance("Set tracking cut in the absorber");
  fTrackingCutCmd->SetParameterName("Cut",false);
  fTrackingCutCmd->SetRange("Cut>0.");
  fTrackingCutCmd->SetUnitCategory("Energy");
  fTrackingCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTrackingCutCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{

  delete fNPMaterCmd;
  delete fNReplicaRCmd;
  delete fNReplicaAzmCmd;
  delete fAbsRadiusCmd;
  delete fNPRadiusCmd;
  delete fTrackingCutCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fNPMaterCmd )
   { fDetector->SetNPMaterial(newValue);}
   
  if( command == fNReplicaRCmd )
   { fDetector->SetNReplicaR(fNReplicaRCmd->GetNewIntValue(newValue));}

  if( command == fNReplicaAzmCmd )
   { fDetector->SetNReplicaAzm(fNReplicaAzmCmd->GetNewIntValue(newValue));}

  if( command == fAbsRadiusCmd )
   { fDetector->SetAbsRadius(fAbsRadiusCmd->GetNewDoubleValue(newValue));}

  if( command == fNPRadiusCmd )
   { fDetector->SetNPRadius(fNPRadiusCmd->GetNewDoubleValue(newValue));}
   
  if( command == fTrackingCutCmd )
   { fDetector->SetTrackingCut(fTrackingCutCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
