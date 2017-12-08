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
// $Id: G4VEmFluctuationModel.cc 106208 2017-09-20 01:53:57Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEmFluctuationModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.07.2005
//
// Modifications:
//
//
// Class Description: 
//
// Abstract class for interface to simualtion of energy loss fluctuations

// -------------------------------------------------------------------
//

#include "G4VEmFluctuationModel.hh"
#include "G4LossTableManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmFluctuationModel::G4VEmFluctuationModel(const G4String& nam)
  : name(nam) 
{
  fManager = G4LossTableManager::Instance();
  fManager->Register(this);
}

G4VEmFluctuationModel::~G4VEmFluctuationModel() 
{
  fManager->DeRegister(this);
}

void G4VEmFluctuationModel::InitialiseMe(const G4ParticleDefinition*)
{}

void G4VEmFluctuationModel::SetParticleAndCharge(const G4ParticleDefinition*, 
                                                 G4double)
{}

