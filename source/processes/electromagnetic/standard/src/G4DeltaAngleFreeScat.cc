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
// File name:     G4DeltaAngleFreeScat
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

#include "G4DeltaAngleFreeScat.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

using namespace std;

G4DeltaAngleFreeScat::G4DeltaAngleFreeScat(const G4String&)
  : G4VEmAngularDistribution("deltaFree")
{}    

G4DeltaAngleFreeScat::~G4DeltaAngleFreeScat() = default;

G4ThreeVector& 
G4DeltaAngleFreeScat::SampleDirection(const G4DynamicParticle* dp,
			      G4double kinEnergyFinal, G4int, 
			      const G4Material*)
{
  G4double deltaMomentum = 
    sqrt(kinEnergyFinal*(kinEnergyFinal + 2*electron_mass_c2));

  G4double costet = kinEnergyFinal*(dp->GetTotalEnergy() + electron_mass_c2) /
    (deltaMomentum * dp->GetTotalMomentum());

  G4double phi = G4UniformRand()*twopi;
  G4double sintet = sqrt((1 - costet)*(1 + costet));
 
  fLocalDirection.set(sintet*cos(phi), sintet*sin(phi), costet);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

void G4DeltaAngleFreeScat::PrintGeneratorInformation() const
{} 
