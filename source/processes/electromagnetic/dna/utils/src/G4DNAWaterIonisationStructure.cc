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

#include "G4DNAWaterIonisationStructure.hh"
#include "G4SystemOfUnits.hh"

G4DNAWaterIonisationStructure::G4DNAWaterIonisationStructure(): nLevels(5)
{
  energyConstant.push_back(10.79*eV);
  energyConstant.push_back(13.39*eV);
  energyConstant.push_back(16.05*eV);
  energyConstant.push_back(32.30*eV);
  energyConstant.push_back(539.0*eV);

  nLevels = (G4int)energyConstant.size();
}


G4DNAWaterIonisationStructure::~G4DNAWaterIonisationStructure()
{ }
 

G4double G4DNAWaterIonisationStructure::IonisationEnergy(G4int level)
{
  G4double ionisation = 0.;

  if (level >=0 && level < nLevels) ionisation = energyConstant[level];

  return ionisation;
}
