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
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#include "G4DNACPA100WaterIonisationStructure.hh"
#include "G4SystemOfUnits.hh"

G4DNACPA100WaterIonisationStructure::G4DNACPA100WaterIonisationStructure(): nLevels(5)
{
  energyConstant.push_back(10.79*eV);
  energyConstant.push_back(13.39*eV);
  energyConstant.push_back(16.05*eV);
  energyConstant.push_back(32.30*eV);
  energyConstant.push_back(539.0*eV);

  UConstant.push_back(61.91*eV);
  UConstant.push_back(59.52*eV);
  UConstant.push_back(48.36*eV);
  UConstant.push_back(70.71*eV);
  UConstant.push_back(796.2*eV);

  nLevels = energyConstant.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACPA100WaterIonisationStructure::~G4DNACPA100WaterIonisationStructure()
{ }
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100WaterIonisationStructure::IonisationEnergy(G4int level)
{
  G4double ionisation = 0.;

  if (level >=0 && level < nLevels) ionisation = energyConstant[level];

  return ionisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNACPA100WaterIonisationStructure::UEnergy(G4int level)
{
  G4double ionisation = 0.;

  if (level >=0 && level < nLevels) ionisation = UConstant[level];

  return ionisation;
}

