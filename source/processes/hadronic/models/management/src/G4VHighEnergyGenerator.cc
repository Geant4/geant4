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
//      V. Uzhinsky Nov. 2012
//   Added method GetProjectileNucleus for simulation of nucleus-nucleus inter.
//
// G4VHighEnergyGenerator
#include "G4VHighEnergyGenerator.hh"

G4VHighEnergyGenerator::G4VHighEnergyGenerator(const G4String& modelName)
  : G4HadronicInteraction(modelName)
{}

G4VHighEnergyGenerator::~G4VHighEnergyGenerator()
{}

void G4VHighEnergyGenerator::ModelDescription(std::ostream& outFile) const
{
  outFile << " Parton-string models description not written yet \n";
}

G4V3DNucleus * G4VHighEnergyGenerator::GetProjectileNucleus() const
{                                                                     
  G4ExceptionDescription ed;
  ed << "The used HighEnergyGenerator "<<GetModelName()
     <<" cannot manage with a residual projectile nucleus";
  G4Exception("G4VHighEnergyGenerator::GetProjectileNucleus ", "G4had_mod_man",
                FatalException, ed); 
  return 0;
}
