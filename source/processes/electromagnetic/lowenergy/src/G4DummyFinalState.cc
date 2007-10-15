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
// $Id: G4DummyFinalState.cc,v 1.2 2007-10-15 08:36:35 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Geant4-DNA dummy final state for test purpose
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#include "G4DummyFinalState.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleMomentum.hh"

G4DummyFinalState::G4DummyFinalState()
{
  // These data members will be used in the next implementation iteration, 
  // when the enriched PhysicsModel policy is implemented
  name = "dummyFinalState";
  lowEnergyLimit = 7.4 * eV;
  highEnergyLimit = 10 * MeV;
}


G4DummyFinalState::~G4DummyFinalState()
{ 
  // empty
}
 

const G4FinalStateProduct& G4DummyFinalState::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  // Clear previous secondaries, energy deposit and particle kill status
  product.Clear();

  // Create fake final state consisting of one electron product with 1 keV kinetic energy 
    
  G4double energy = 1. * keV;

  // Kill incident particle if its energy is too low to create secondary
  G4double primaryEnergy = track.GetKineticEnergy();
  if ( primaryEnergy < energy)
    {
      product.KillPrimaryParticle();
      product.AddEnergyDeposit(primaryEnergy);
    }
  else
    {
      // Create a DynamicParticle  
      G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
      G4double initX = 0.; 
      G4double initY = 0.; 
      G4double initZ = 1.;
      G4ParticleMomentum direction(initX,initY,initZ);
      G4DynamicParticle* dynamicParticle = new G4DynamicParticle(electron,direction,energy);
      product.AddSecondary(dynamicParticle);
    }

  return product;
}



