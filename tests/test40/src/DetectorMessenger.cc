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
// $Id: DetectorMessenger.cc,v 1.2 2006-06-29 22:03:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
:Detector(Det)
{
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");

  MaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  MaterCmd->SetGuidance("Select Material.");
  MaterCmd->SetParameterName("material",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  LBinCmd = new G4UIcmdWith3Vector("/testem/det/setLbin",this);
  LBinCmd->SetGuidance("set longitudinal bining");
  LBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  LBinCmd->SetParameterName("nLtot","dLradl"," ",true);
  LBinCmd->SetRange("nLtot>=1 && dLradl>0");
  LBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RBinCmd = new G4UIcmdWith3Vector("/testem/det/setRbin",this);
  RBinCmd->SetGuidance("set radial bining");
  RBinCmd->SetGuidance("nb of bins; bin thickness (in radl)");
  RBinCmd->SetParameterName("nRtot","dRradl"," ",true);
  RBinCmd->SetRange("nRtot>=1 && dRradl>0");
  RBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  accCmd = new G4UIcmdWith3Vector("/testem/det/acceptance",this);
  accCmd->SetGuidance("set Edep and RMS");
  accCmd->SetGuidance("acceptance values");
  accCmd->SetParameterName("edep","rms","limit",true);
  accCmd->SetRange("edep>0 && edep<1 && rms>0");
  accCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);
  FieldCmd->SetGuidance("Define magnetic field.");
  FieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  FieldCmd->SetParameterName("Bz",false);
  FieldCmd->SetUnitCategory("Magnetic flux density");
  FieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  factoryCmd = new G4UIcmdWithAString("/testem/histo/fileName",this);
  factoryCmd->SetGuidance("set name for the histograms file");

  fileCmd = new G4UIcmdWithAString("/testem/histo/fileType",this);
  fileCmd->SetGuidance("set type (hbook, XML) for the histograms file");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete LBinCmd;
  delete RBinCmd;
  delete accCmd;
  delete FieldCmd;
  delete UpdateCmd;
  delete testemDir;
  delete factoryCmd;
  delete fileCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == MaterCmd )
   { Detector->SetMaterial(newValue);}
   
  if( command == LBinCmd )
   { Detector->SetLBining(LBinCmd->GetNew3VectorValue(newValue));}

  if( command == RBinCmd )
   { Detector->SetRBining(RBinCmd->GetNew3VectorValue(newValue));}

  if( command == accCmd )
   { Detector->SetEdepAndRMS(accCmd->GetNew3VectorValue(newValue));}

  if( command == FieldCmd )
   { Detector->SetMagField(FieldCmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}

  if (command == factoryCmd) Detector->SetHistoName(newValue);

  if (command == fileCmd) Detector->SetHistoType(newValue);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
