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
// $Id: G4DNARuddAngle.cc 68380 2013-03-22 18:39:29Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4DNARuddAngle
//
// Author:        Vladimir Ivantcheko
// 
// Creation date: 23 August 2013
//
// Modifications: 
//
// Class Description: 
//
// Delta-electron Angular Distribution Generation 
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#include "G4DNARuddAngle.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

G4DNARuddAngle::G4DNARuddAngle(const G4String&)
  : G4VEmAngularDistribution("deltaRudd")
{
  fElectron = G4Electron::Electron();
}    

G4DNARuddAngle::~G4DNARuddAngle() 
{}

G4ThreeVector& 
G4DNARuddAngle::SampleDirectionForShell(const G4DynamicParticle* dp,
					G4double secKinetic, G4int, 
					G4int, 
					const G4Material*)
{
  G4double k = dp->GetKineticEnergy();
  G4double cosTheta = 1.0;

  const G4ParticleDefinition* particle = dp->GetDefinition();
  G4double mass = particle->GetPDGMass();

  G4double maximumEnergyTransfer = k;
  if(particle == fElectron) { maximumEnergyTransfer *= 0.5; }
  else if(mass > MeV) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau  = k/mass;
    maximumEnergyTransfer = 2.0*electron_mass_c2*tau*(tau + 2.) /
      (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
    
  }

  if (secKinetic>100*eV && secKinetic <= maximumEnergyTransfer) {
    cosTheta = std::sqrt(secKinetic / maximumEnergyTransfer);
  } else {
    cosTheta = (2.*G4UniformRand())-1.;
  }

  G4double sint = sqrt((1.0 - cosTheta)*(1.0 + cosTheta));
  G4double phi  = twopi*G4UniformRand(); 

  fLocalDirection.set(sint*cos(phi), sint*sin(phi), cosTheta);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

G4ThreeVector& 
G4DNARuddAngle::SampleDirection(const G4DynamicParticle* dp,
				G4double secEkin, G4int Z, 
				const G4Material* mat)
{
  return SampleDirectionForShell(dp, secEkin, Z, 0, mat);
}

void G4DNARuddAngle::PrintGeneratorInformation() const
{} 
