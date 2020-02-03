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
/// \file ExGflashMessenger.cc
/// \brief Implementation of the ExGflashMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExGflashMessenger.hh"
#include "ExGflashDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashMessenger::ExGflashMessenger(ExGflashDetectorConstruction * Det)
:G4UImessenger(),fDetector(Det)
{
  fExGflashDir = new G4UIdirectory("/exgflash/");
  fExGflashDir->SetGuidance(" Gflash example commands.");
  
  fVerbose = new G4UIcmdWithAnInteger("/exgflash/verbose",this);
  fVerbose->SetGuidance("set exglash verbosity");
  fVerbose->SetGuidance("0- silent, 1 - on exit, 2 - run , 3 - event");
  fVerbose->SetParameterName("ver",false);
  fVerbose->SetRange("ver >= 0");
  fVerbose->AvailableForStates(G4State_PreInit,G4State_Idle);
  fVerbose->SetToBeBroadcasted(false);
  
  fDetDir = new G4UIdirectory("/exgflash/det/");
  fDetDir->SetGuidance("detector construction commands");

  fMaterCmd = new G4UIcmdWithAString("/exgflash/det/setMat",this);
  fMaterCmd->SetGuidance("Select Material.");
  fMaterCmd->SetParameterName("material",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fMaterCmd->SetToBeBroadcasted(false);

  fLBinCmd = new G4UIcmdWith3Vector("/exgflash/det/setLbin",this);
  fLBinCmd->SetGuidance("set longitudinal bining");
  fLBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  // must have omitable==true we use 2 values only
  fLBinCmd->SetParameterName("nLtot","dLradl"," ",true);
  fLBinCmd->SetRange("nLtot>=1 && dLradl>0");
  fLBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fLBinCmd->SetToBeBroadcasted(false);

  fRBinCmd = new G4UIcmdWith3Vector("/exgflash/det/setRbin",this);
  fRBinCmd->SetGuidance("set radial bining");
  fRBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  fRBinCmd->SetParameterName("nRtot","dRradl"," ",true);
  fRBinCmd->SetRange("nRtot>=1 && dRradl>0");
  fRBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fRBinCmd->SetToBeBroadcasted(false);    
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashMessenger::~ExGflashMessenger()
{
  delete fRBinCmd;
  delete fLBinCmd;
  delete fMaterCmd;
  delete fDetDir;
  delete fVerbose;
  delete fExGflashDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fMaterCmd )
    { fDetector->SetMaterial(newValue);}
  
  if( command == fVerbose )
    { fDetector->SetVerbose(fVerbose->GetNewIntValue(newValue));}
  
  if( command == fLBinCmd )
    { fDetector->SetLBining(fLBinCmd->GetNew3VectorValue(newValue));}
  
  if( command == fRBinCmd )
    { fDetector->SetRBining(fRBinCmd->GetNew3VectorValue(newValue));}  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
