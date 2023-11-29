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
// GEANT4 Class file
//
// File name:     G4EmStandUtil
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 29.05.2022
//

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmStandUtil.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UrbanFluctuation.hh"
#include "G4LossFluctuationDummy.hh"
#include "G4IonFluctuations.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmFluctuationModel* G4EmStandUtil::ModelOfFluctuations(G4bool isIon)
{
  G4VEmFluctuationModel* f = nullptr;
  auto ftype = G4EmParameters::Instance()->FluctuationType();
  if (ftype == fDummyFluctuation) {
    f = new G4LossFluctuationDummy();
  } else if (isIon) {
    f = new G4IonFluctuations();
  } else if (ftype == fUrbanFluctuation) {
    f = new G4UrbanFluctuation();
  } else {
    f = new G4UniversalFluctuation();
  }
  return f;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
