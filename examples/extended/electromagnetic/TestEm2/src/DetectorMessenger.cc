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
/// \file electromagnetic/TestEm2/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc,v 1.8 2006-06-29 16:50:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:fDetector(Det)
{
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance(" detector control.");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");

  fMaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  fMaterCmd->SetGuidance("Select Material.");
  fMaterCmd->SetParameterName("material",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fLBinCmd = new G4UIcmdWith3Vector("/testem/det/setLbin",this);
  fLBinCmd->SetGuidance("set longitudinal bining");
  fLBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  fLBinCmd->SetParameterName("nLtot","dLradl"," ",true);
  fLBinCmd->SetRange("nLtot>=1 && dLradl>0");
  fLBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fRBinCmd = new G4UIcmdWith3Vector("/testem/det/setRbin",this);
  fRBinCmd->SetGuidance("set radial bining");
  fRBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  fRBinCmd->SetParameterName("nRtot","dRradl"," ",true);
  fRBinCmd->SetRange("nRtot>=1 && dRradl>0");
  fRBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);
  fFieldCmd->SetGuidance("Define magnetic field.");
  fFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fFieldCmd->SetParameterName("Bz",false);
  fFieldCmd->SetUnitCategory("Magnetic flux density");
  fFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fLBinCmd;
  delete fRBinCmd;
  delete fFieldCmd;
  delete fUpdateCmd;
  delete fDetDir;
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetMaterial(newValue);}
   
  if( command == fLBinCmd )
   { fDetector->SetLBining(fLBinCmd->GetNew3VectorValue(newValue));}

  if( command == fRBinCmd )
   { fDetector->SetRBining(fRBinCmd->GetNew3VectorValue(newValue));}

  if( command == fFieldCmd )
   { fDetector->SetMagField(fFieldCmd->GetNewDoubleValue(newValue));}

  if( command == fUpdateCmd )
   { fDetector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
