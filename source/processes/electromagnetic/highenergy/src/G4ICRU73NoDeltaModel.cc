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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ICRU73NoDeltaModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 04.06.2010
//
// Modifications:
//
//
// Class Description:
//
// Ionisation of heavy negatively charged particles

// -------------------------------------------------------------------
//

#include "G4ICRU73NoDeltaModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ICRU73NoDeltaModel::G4ICRU73NoDeltaModel(const G4ParticleDefinition*p,
					   const G4String& nam) :
  G4ICRU73QOModel(p, nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ICRU73NoDeltaModel::~G4ICRU73NoDeltaModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ICRU73NoDeltaModel::ComputeDEDXPerVolume(
                            const G4Material* material,
			    const G4ParticleDefinition* pd,
                            G4double kinEnergy, G4double)
{
  return 
    G4ICRU73QOModel::ComputeDEDXPerVolume(material, pd, kinEnergy, DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ICRU73NoDeltaModel::CrossSectionPerVolume(
                            const G4Material*,
			    const G4ParticleDefinition*,
			    G4double, G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


