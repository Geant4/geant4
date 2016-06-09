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
// $Id: G4FinalStateExcitationBorn.cc,v 1.2 2007/11/09 20:11:04 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#include "G4FinalStateExcitationBorn.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleMomentum.hh"

G4FinalStateExcitationBorn::G4FinalStateExcitationBorn()
{
  name = "FinalStateExcitationBorn";
  lowEnergyLimit = 7.4 * eV;
  highEnergyLimit = 10 * MeV;
}


G4FinalStateExcitationBorn::~G4FinalStateExcitationBorn()
{ 
  // empty
  // G4DynamicParticle objects produced are owned by client
}
 

const G4FinalStateProduct& G4FinalStateExcitationBorn::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  // Clear previous secondaries, energy deposit and particle kill status
  product.Clear();

  const G4DynamicParticle* particle = track.GetDynamicParticle();

  // Kinetic energy of primary particle
  G4double k = particle->GetKineticEnergy();

  // Select excitation level on the basis of partial excitation cross section
  G4int level = cross.RandomSelect(k);
  // Excitation energy corresponding to the selected level
  G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
  G4double newEnergy = k - excitationEnergy;
  
  if (newEnergy > lowEnergyLimit)
    {
      // Deposit excitation energy locally, modify primary energy accordingly 
      // Particle direction is unchanged
      product.ModifyPrimaryParticle(particle->GetMomentumDirection(),newEnergy);
      product.AddEnergyDeposit(excitationEnergy);
    }
  else
    {
      // Primary particle is killed
      product.KillPrimaryParticle();
    }

  return product;
}


