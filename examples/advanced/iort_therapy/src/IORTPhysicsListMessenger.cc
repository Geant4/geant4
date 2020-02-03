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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "IORTPhysicsListMessenger.hh"

#include "IORTPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

/////////////////////////////////////////////////////////////////////////////
IORTPhysicsListMessenger::IORTPhysicsListMessenger(IORTPhysicsList* pPhys)
:pPhysicsList(pPhys)
{
  physDir = new G4UIdirectory("/Physics/");
  physDir->SetGuidance("Commands to activate physics models and set cuts");
   
  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setGCut",this);  
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setECut",this);  
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  protoCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setPCut",this);  
  protoCutCmd->SetGuidance("Set positron cut.");
  protoCutCmd->SetParameterName("Pcut",false);
  protoCutCmd->SetUnitCategory("Length");
  protoCutCmd->SetRange("Pcut>0.0");
  protoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setCuts",this);  
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut",false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>0.0");
  allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  allDetectorCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setDetectorCuts",this);  
  allDetectorCmd->SetGuidance("Set cut for all. into Detector");
  allDetectorCmd->SetParameterName("cut",false);
  allDetectorCmd->SetUnitCategory("Length");
  allDetectorCmd->SetRange("cut>0.0");
  allDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  pListCmd = new G4UIcmdWithAString("/Physics/addPhysics",this);  
  pListCmd->SetGuidance("Add physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit, G4State_Idle);  
}

/////////////////////////////////////////////////////////////////////////////
IORTPhysicsListMessenger::~IORTPhysicsListMessenger()
{
  delete gammaCutCmd;
  delete electCutCmd;
  delete protoCutCmd;
  delete allCutCmd;
  delete allDetectorCmd;
  delete pListCmd;
  delete physDir;    
}

/////////////////////////////////////////////////////////////////////////////
void IORTPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  if( command == gammaCutCmd )
   { pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));}
     
  else if( command == electCutCmd )
   { pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));}
     
  else if( command == protoCutCmd )
   { pPhysicsList->SetCutForPositron(protoCutCmd->GetNewDoubleValue(newValue));}

  else if( command == allCutCmd )
    {
      G4double cut = allCutCmd->GetNewDoubleValue(newValue);
      pPhysicsList->SetCutForGamma(cut);
      pPhysicsList->SetCutForElectron(cut);
      pPhysicsList->SetCutForPositron(cut);
    } 
  else if( command == allDetectorCmd)
  {
      G4double cut = allDetectorCmd -> GetNewDoubleValue(newValue);
      pPhysicsList -> SetDetectorCut(cut);
  }
  else if( command == pListCmd )
   { pPhysicsList->AddPhysicsList(newValue);}
}
