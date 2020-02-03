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
/// \file medical/fanoCavity/src/MyKleinNishinaMessenger.cc
/// \brief Implementation of the MyKleinNishinaMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MyKleinNishinaMessenger.hh"

#include "MyKleinNishinaCompton.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyKleinNishinaMessenger::MyKleinNishinaMessenger(MyKleinNishinaCompton* pPhys)
:fKleinNishina(pPhys), fCsFactor(0)
{  
  fCsFactor = new G4UIcmdWithADouble("/testem/phys/crossSectionFactor",this);
  fCsFactor->SetGuidance("multiply Compton cross section");
  fCsFactor->SetParameterName("factor",false);
  fCsFactor->SetRange("factor>=0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyKleinNishinaMessenger::~MyKleinNishinaMessenger()
{
  delete fCsFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyKleinNishinaMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if (command == fCsFactor)
   {fKleinNishina->SetCSFactor(fCsFactor->GetNewDoubleValue(newValue));}        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
