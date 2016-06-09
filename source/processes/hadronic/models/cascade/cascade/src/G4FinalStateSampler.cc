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
// $Id: G4FinalStateSampler.cc,v 1.2 2009/09/17 18:15:38 dennis Exp $
// GEANT4 tag $Name: geant4-09-03 $
//

#include "G4FinalStateSampler.hh"
#include "Randomize.hh"


std::pair<G4int, G4double>
G4FinalStateSampler::interpolateEnergy(G4double e) const
{
  G4int index = 29;
  G4double fraction = 0.0;

  for (G4int i = 1; i < 30; i++) {
    if (e < energyScale[i]) {
      index = i-1;
      fraction = (e - energyScale[index]) / (energyScale[i] - energyScale[index]);
      break;
    }
  }
  return std::pair<G4int, G4double>(index, fraction);
}


G4int
G4FinalStateSampler::sampleFlat(std::vector<G4double> sigma) const
{
  G4int i;
  G4double sum(0.);
  for (i = 0; i < G4int(sigma.size()); i++) sum += sigma[i];

  G4double fsum = sum*G4UniformRand();
  G4double partialSum = 0.0;
  G4int channel = 0;

  for (i = 0; i < G4int(sigma.size()); i++) {
    partialSum += sigma[i];
    if (fsum < partialSum) {
      channel = i;
      break;
    }
  }

  return channel;
}


void G4FinalStateSampler::CheckQnums(G4FastVector<G4ReactionProduct,256> &vec,
                                G4int &vecLen,
                                G4ReactionProduct &currentParticle,
                                G4ReactionProduct &targetParticle,
                                G4double Q, G4double B, G4double S)
{
  G4ParticleDefinition* projDef = currentParticle.GetDefinition();
  G4ParticleDefinition* targDef = targetParticle.GetDefinition();
  G4double chargeSum = projDef->GetPDGCharge() + targDef->GetPDGCharge();
  G4double baryonSum = projDef->GetBaryonNumber() + targDef->GetBaryonNumber();
  G4double strangenessSum = projDef->GetQuarkContent(3) - 
                            projDef->GetAntiQuarkContent(3) + 
                            targDef->GetQuarkContent(3) -
                            targDef->GetAntiQuarkContent(3);

  G4ParticleDefinition* secDef = 0;
  for (G4int i = 0; i < vecLen; i++) {
    secDef = vec[i]->GetDefinition();
    chargeSum += secDef->GetPDGCharge();
    baryonSum += secDef->GetBaryonNumber();
    strangenessSum += secDef->GetQuarkContent(3) 
                    - secDef->GetAntiQuarkContent(3);
  }

  G4bool OK = true;
  if (chargeSum != Q) {
    G4cout << " Charge not conserved " << G4endl;
    OK = false;
  }
  if (baryonSum != B) {
    G4cout << " Baryon number not conserved " << G4endl;
    OK = false;
  }
  if (strangenessSum != S) {
    G4cout << " Strangeness not conserved " << G4endl;
    OK = false;
  } 

  if (!OK) {
    G4cout << " projectile: " << projDef->GetParticleName() 
           << "  target: " << targDef->GetParticleName() << G4endl;
    for (G4int i = 0; i < vecLen; i++) {
      secDef = vec[i]->GetDefinition();
      G4cout << secDef->GetParticleName() << " " ;
    }
    G4cout << G4endl;
  }

}


const G4double G4FinalStateSampler::energyScale[30] = {
  0.0,  0.01, 0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075, 0.1,
  0.13, 0.18, 0.24,  0.32,  0.42,  0.56,  0.75,  1.0,   1.3,   1.8,
  2.4,  3.2,  4.2,   5.6,   7.5,   10.0,  13.0,  18.0,  24.0, 32.0 };

/* end of file */
