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
// $Id: G4FinalStateChargeIncrease.cc,v 1.2 2007/11/09 20:11:04 pia Exp $
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


#include "G4FinalStateChargeIncrease.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleMomentum.hh"
#include "G4DNAGenericIonsManager.hh"

G4FinalStateChargeIncrease::G4FinalStateChargeIncrease()
{
  name = "ChargeIncrease";
  lowEnergyLimit = 1 * keV;
  highEnergyLimit = 10 * MeV;
}


G4FinalStateChargeIncrease::~G4FinalStateChargeIncrease()
{}
 

const G4FinalStateProduct& G4FinalStateChargeIncrease::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  // Clear previous secondaries, energy deposit and particle kill status
  product.Clear();

  // Primary particle
  product.KillPrimaryParticle();
  product.AddEnergyDeposit(0.);
  G4ParticleDefinition* definition = track.GetDefinition();
 
  // Secondaries
  G4double inK = track.GetDynamicParticle()->GetKineticEnergy();
  
  G4int finalStateIndex = cross.RandomSelect(inK,definition);

  G4int n = NumberOfFinalStates(track.GetDefinition(),finalStateIndex);
  G4double outK = inK - IncomingParticleBindingEnergyConstant(definition,finalStateIndex);

  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4double electronK;
  if (definition == instance->GetIon("hydrogen")) electronK = inK*electron_mass_c2/proton_mass_c2;
  else electronK = inK*electron_mass_c2/(3728*MeV);
  
  if (outK<0)
    {
      G4String message;
      message="G4FinalStateChargeIncrease - Final kinetic energy is below 0! Process ";
    }
  
  product.AddSecondary(new G4DynamicParticle(OutgoingParticleDefinition(definition,finalStateIndex), 
					     track.GetDynamicParticle()->GetMomentumDirection(), 
					     outK));

  n = n - 1;
  
  while (n>0)
    {
      n--;
   
      product.AddSecondary
	(new G4DynamicParticle(G4Electron::Electron(), track.GetDynamicParticle()->GetMomentumDirection(), electronK));
    }
        
  return product;
}

G4int G4FinalStateChargeIncrease::NumberOfFinalStates(G4ParticleDefinition* particleDefinition, 
						      G4int finalStateIndex )
 
{
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();
  if (particleDefinition == instance->GetIon("hydrogen")) return 2;
  if (particleDefinition == instance->GetIon("alpha+")) return 2;
  
  if (particleDefinition == instance->GetIon("helium")) 
    {    if (finalStateIndex==0) return 2;
    return 3;
    }
  return 0;
}


G4ParticleDefinition*  G4FinalStateChargeIncrease::OutgoingParticleDefinition (G4ParticleDefinition* particleDefinition, 
									       G4int finalStateIndex) 
{
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  if (particleDefinition == instance->GetIon("hydrogen")) return G4Proton::Proton(); 
  if (particleDefinition == instance->GetIon("alpha+")) return instance->GetIon("alpha++");

  if (particleDefinition == instance->GetIon("alpha+")) return instance->GetIon("helium");
  {
    if (finalStateIndex==0) return instance->GetIon("alpha+"); 
    return instance->GetIon("alpha++");
  }
  return 0;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateChargeIncrease::IncomingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition, 
									   G4int finalStateIndex )
{
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  if(particleDefinition == instance->GetIon("hydrogen")) return 13.6*eV;
  
  if(particleDefinition == instance->GetIon("alpha+"))
    {
      // Binding energy for    He+ -> He++ + e-    54.509 eV
      // Binding energy for    He  -> He+  + e-    24.587 eV
      return 54.509*eV;
    }
   
  if(particleDefinition == instance->GetIon("helium"))
    {
      // Binding energy for    He+ -> He++ + e-    54.509 eV
      // Binding energy for    He  -> He+  + e-    24.587 eV

      if (finalStateIndex==0) return 24.587*eV;
      return (54.509 + 24.587)*eV;
    }  

  return 0;
}


