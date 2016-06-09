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
// $Id$
//
// G4VHighEnergyGenerator
#include "G4VHighEnergyGenerator.hh"
#include "G4HadronicException.hh"

G4VHighEnergyGenerator::G4VHighEnergyGenerator(const G4String& modelName)
:   epCheckLevels(DBL_MAX,DBL_MAX)
{
	theGeneratorModelName=modelName;
}


G4VHighEnergyGenerator::~G4VHighEnergyGenerator()
{
}

std::pair<G4double, G4double> G4VHighEnergyGenerator::GetEnergyMomentumCheckLevels() const
{
   return epCheckLevels;
}

void G4VHighEnergyGenerator::SetEnergyMomentumCheckLevels(
									G4double relativeLevel, G4double absoluteLevel)
{
	epCheckLevels.first=relativeLevel;
	epCheckLevels.second=absoluteLevel;
}

void G4VHighEnergyGenerator::ModelDescription(std::ostream& outFile) const
{
  outFile << " Parton-string models description not written yet \n";
}

G4String G4VHighEnergyGenerator::GetModelName() const
{
	return theGeneratorModelName;
}
