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
// Class Description - End

class G4ParaFissionModel : public G4HadronicInteraction
{
public:
  G4ParaFissionModel()
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 60.*MeV );
  }
  
  virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
                                              G4Nucleus& theNucleus)
  {
    theParticleChange.Initialize(aTrack);
    theParticleChange.SetStatusChange( fStopAndKill );
    theParticleChange.SetEnergyChange( 0.0 );
    
    // prepare the fragment
    G4Fragment anInitialState;
     G4int anA=theNucleus.GetN();
     G4int aZ=theNucleus.GetZ();
     G4double nucMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ ,anA);
     
     anA += aTrack.GetDynamicParticle()->GetDefinition()->GetBaryonNumber();
     aZ += aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge();
     
     G4int numberOfEx = aTrack.GetDynamicParticle()->GetDefinition()->GetBaryonNumber();
     G4int numberOfCh = abs(aTrack.GetDynamicParticle()->GetDefinition()->GetPDGCharge());
     G4int numberOfHoles = 0;
     G4double exEnergy = 0;
     
     G4ThreeVector exciton3Momentum = aTrack.GetMomentum();
     G4double compoundMass = aTrack.GetTotalEnergy();
     compoundMass += nucMass;
     compoundMass = sqrt(compoundMass*compoundMass - exciton3Momentum*exciton3Momentum);
     G4LorentzVector fragment4Momentum(exciton3Momentum, 
                                      sqrt(exciton3Momentum.mag2()+compoundMass*compoundMass));
    
    anInitialState.SetA(anA);
    anInitialState.SetZ(aZ);
    anInitialState.SetNumberOfCharged(numberOfCh);
    anInitialState.SetNumberOfHoles(numberOfHoles);
    anInitialState.SetNumberOfExcitons(numberOfEx);
    anInitialState.SetMomentum(fragment4Momentum);

    // do the fission
    G4FragmentVector * theFissionResult = theFission.BreakUp(anInitialState);
    
    // deexcite the fission fragments and fill result
    vector<G4DynamicParticle *> theResult;
    G4int ll = theFissionResult->length();
    for(G4int i=0; i<ll; i++)
    {
      G4ReactionProductVector * theExcitationResult = 0; 
      if(theFissionResult->at(i)->GetExcitationEnergy()>1.*eV)
      {
        G4Fragment * aFragment = theFissionResult->at(i);
	G4double exenergy = aFragment->GetExcitationEnergy();
	theExcitationResult = theHandler.BreakItUp(*(theFissionResult->at(i)));
	// add secondaries
	for(G4int j=0; j<theExcitationResult->length(); j++)
	{
          G4DynamicParticle* p0 = new G4DynamicParticle;
          p0->SetDefinition( theExcitationResult->at(j)->GetDefinition() );
          p0->SetMomentum( theExcitationResult->at(j)->GetMomentum() );
          theResult.push_back(p0);
	}
      }
      else
      {
        // add secondary
	G4DynamicParticle* p0 = new G4DynamicParticle;
	p0->SetDefinition(theFissionResult->at(i)->GetParticleDefinition());
	p0->SetMomentum(theFissionResult->at(i)->GetMomentum().vect());
        theResult.push_back(p0);
      }
    }
    
    // fill particle change
    for(G4int k=0; k<theResult.size(); k++)
    {
      theParticleChange.AddSecondary(theResult[k]);
    }
    
    // return
    return &theParticleChange;
  }
private:

  G4CompetitiveFission theFission;
  G4ExcitationHandler theHandler;
  
  G4ParticleChange theParticleChange;
};
#endif
