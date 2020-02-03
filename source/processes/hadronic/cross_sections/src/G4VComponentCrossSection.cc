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
// File name:    G4VComponentCrossSection
//
// Authors:  G.Folger, V.Ivanchenko, D.Wright
//
// Modifications:
//

#include "G4VComponentCrossSection.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include <iostream>

G4VComponentCrossSection::G4VComponentCrossSection(const G4String& nam) :
  verboseLevel(0),minKinEnergy(0.0),maxKinEnergy(DBL_MAX),name(nam) 
{
  registry = G4CrossSectionDataSetRegistry::Instance();
  registry->Register(this);
}

G4VComponentCrossSection::~G4VComponentCrossSection()
{
  registry->DeRegister(this);
}

G4double 
G4VComponentCrossSection::ComputeQuasiElasticRatio(const G4ParticleDefinition*,
                                                   G4double /*kinEnergy*/, 
						   G4int /*Z*/, G4int /*N*/)
{
  return 0.0;
}

void G4VComponentCrossSection::Description(std::ostream&) const
{}

void 
G4VComponentCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{}

void 
G4VComponentCrossSection::DumpPhysicsTable(const G4ParticleDefinition&)
{}
