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
// Geant4 class : G4LowEHadronElastic
//
// Author : V.Ivanchenko 10 May 2019
//  
//

#include "G4LowEHadronElastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"

G4LowEHadronElastic::G4LowEHadronElastic():G4HadronElastic("hLowEElastic")
{
  plabLowLimit  = 400*CLHEP::MeV;
  plabHighLimit = 2000*CLHEP::MeV;
}

G4LowEHadronElastic::~G4LowEHadronElastic()
{}

G4double 
G4LowEHadronElastic::SampleInvariantT(const G4ParticleDefinition* p, 
				      G4double plab, G4int Z, G4int A)
{
  return (IsResonanseScattering(p, plab, Z, A)) 
    ? G4UniformRand()*pLocalTmax 
    : G4HadronElastic::SampleInvariantT(p, plab, Z, A);
}

G4bool
G4LowEHadronElastic::IsResonanseScattering(const G4ParticleDefinition*, 
				           G4double plab, 
                                           G4int, G4int)
{
  return (plab < plabHighLimit);
}
