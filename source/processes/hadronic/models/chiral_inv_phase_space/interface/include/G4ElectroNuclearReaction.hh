//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4ElectroNuclearReaction_h
#define G4ElectroNuclearReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4ElectroNuclearCrossSection.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4ElectroNuclearReaction : public G4HadronicInteraction
{
  public: 
    virtual ~G4ElectroNuclearReaction()
    {
    }
    
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
    G4ElectroNuclearCrossSection theData;
    G4ParticleChange theResult;
};

inline
G4VParticleChange * G4ElectroNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  const G4ParticleDefinition* aD = aTrack.GetDynamicParticle()->GetDefinition();
  if((aD != G4Electron::ElectronDefinition()) && (aD != G4Positron::PositronDefinition()))
  {
    G4Exception("Called G4ElectroNuclearReaction for particle other than electron or positron");
  }
  
  theResult.Initialize(aTrack);

  const G4ElementTable* aTab = G4Element::GetElementTable();
  G4Element * anElement = 0;
  G4int aZ = static_cast<G4int>(aTargetNucleus.GetZ()+.1);
  for(size_t ii=0; ii<aTab->size(); ii++)
  {
    if ( abs((*aTab)[ii]->GetZ()-aZ) < .1)
    {
      anElement = (*aTab)[ii];
      break;
    }
  }
  if(0==anElement) 
  {
    G4cout << "G4ElectroNuclearReaction::ApplyYourself - trying to react on an element"<<G4endl;
    G4cout << "that is not in the table of elements. Z="<<aTargetNucleus.GetZ()<<G4endl;
    G4Exception("Folding with error.");
  }
  G4double photonEnergy = 10*GeV;  
  G4double xSec;
  while(photonEnergy>3.*GeV)
  {
    xSec = theData.GetCrossSection(aTrack.GetDynamicParticle(), anElement);
    photonEnergy = theData.GetEffectivePhotonEnergy();
  }
  if(aTrack.GetDynamicParticle()->GetKineticEnergy() - photonEnergy < 0)
  {
    G4Exception("G4ElectroNuclearReaction: photonEnergy above electron energy");
  }
  theResult.SetEnergyChange(aTrack.GetDynamicParticle()->GetKineticEnergy() - photonEnergy);
  
  G4ThreeVector photonDirection = aTrack.GetMomentumDirection();
  G4DynamicParticle localGamma(G4Gamma::GammaDefinition(), photonDirection, photonEnergy);
  G4ThreeVector position(0,0,0);
  G4Track localTrack(&localGamma, 0., position);
  G4VParticleChange * result = theModel.ApplyYourself(localTrack, aTargetNucleus, &theResult);
  return result;
}

#endif
