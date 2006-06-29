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
#ifndef G4ParaFissionModel_h
#define G4ParaFissionModel_h

#include "G4CompetitiveFission.hh"
#include "G4ExcitationHandler.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleTable.hh"

// Class Description
// Final state production model for (based on evaluated data
// libraries) description of neutron induced fission below 60 MeV; 
// In case you need the fission fragments, use this model.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.


class G4ParaFissionModel : public G4HadronicInteraction
{
public:
  G4ParaFissionModel()
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 60.*MeV );
  }
  
  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                              G4Nucleus& theNucleus)
  {
    theParticleChange.Clear();
    theParticleChange.SetStatusChange( stopAndKill );
    theParticleChange.SetEnergyChange( 0.0 );
    
    // prepare the fragment

    G4Fragment anInitialState;
    G4double anA = theNucleus.GetN();
    G4double aZ = theNucleus.GetZ();
    G4double nucMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(G4int(aZ) ,G4int(anA));
     
    anA += aTrack.GetDefinition()->GetBaryonNumber();
    aZ += aTrack.GetDefinition()->GetPDGCharge();
     
    G4int numberOfEx = aTrack.GetDefinition()->GetBaryonNumber();
    G4int numberOfCh = G4int(std::abs(aTrack.GetDefinition()->GetPDGCharge()));
    G4int numberOfHoles = 0;
     
    G4ThreeVector exciton3Momentum = aTrack.Get4Momentum().vect();
    G4double compoundMass = aTrack.GetTotalEnergy();
    compoundMass += nucMass;
    compoundMass = std::sqrt(compoundMass*compoundMass - exciton3Momentum*exciton3Momentum);
    G4LorentzVector fragment4Momentum(exciton3Momentum, 
                               std::sqrt(exciton3Momentum.mag2()+compoundMass*compoundMass));
    
    anInitialState.SetA(anA);
    anInitialState.SetZ(aZ);
    anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
    anInitialState.SetNumberOfCharged(numberOfCh);
    anInitialState.SetNumberOfHoles(numberOfHoles);
    anInitialState.SetMomentum(fragment4Momentum);

    // do the fission
    G4FragmentVector * theFissionResult = theFission.BreakUp(anInitialState);
    
    // deexcite the fission fragments and fill result
    std::vector<G4DynamicParticle *> theResult;
    G4int ll = theFissionResult->size();
    for(G4int i=0; i<ll; i++)
    {
      G4ReactionProductVector* theExcitationResult = 0; 
      G4Fragment* aFragment = (*theFissionResult)[i];
      if(aFragment->GetExcitationEnergy()>1.*eV)
      {
	theExcitationResult = theHandler.BreakItUp(*aFragment);
	// add secondaries
	for(G4int j = 0; j < G4int(theExcitationResult->size()); j++)
	{
          G4DynamicParticle* p0 = new G4DynamicParticle;
          p0->SetDefinition( theExcitationResult->operator[](j)->GetDefinition() );
          p0->SetMomentum( theExcitationResult->operator[](j)->GetMomentum() );
          theResult.push_back(p0);
	}
      }
      else
      {
        // add secondary
	G4DynamicParticle* p0 = new G4DynamicParticle;
	p0->SetDefinition(aFragment->GetParticleDefinition());
	p0->SetMomentum(aFragment->GetMomentum().vect());
        theResult.push_back(p0);
      }
    }
    
    // fill particle change
    for(G4int k = 0; k < G4int(theResult.size()); k++)
    {
      theParticleChange.AddSecondary(theResult[k]);
    }
    
    // return
    return &theParticleChange;
  }
private:

  G4CompetitiveFission theFission;
  G4ExcitationHandler theHandler;
  
  G4HadFinalState theParticleChange;
};
#endif
