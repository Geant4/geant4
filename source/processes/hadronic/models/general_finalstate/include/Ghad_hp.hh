#ifndef Ghad_hp_h
#define Ghad_hp_h

#include "G4VIntraNuclearTransportModel.hh"
#include "GhadParticles.hh"
#include "GhadNucleus.hh"
#include "GhadFS.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh"

class Ghad_hp : public G4VIntraNuclearTransportModel
{
  public:
    Ghad_hp() {}
    virtual ~Ghad_hp() {}
    
    virtual G4VParticleChange* 
      ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus) 
      { 
        theResult.Initialize(aTrack);
	theTarget.Initialize(theNucleus);
        G4double maxImpact = theTarget.GetOuterRadius();
	while(theProjectiles.size()<=1) 
	{
	  theProjectiles.Initialize(aTrack, G4ThreeVector(0,0,0));
	  G4Pair<G4double, G4double> impact = theTarget.ChooseImpactXandY(maxImpact);
	  G4ThreeVector aNew(impact.first, impact.second, theProjectiles[0].GetPosition().z());
    	  theProjectiles[0].SetPosition(aNew);
	  Stepping();
	}
	Stepping();
	// fill result ...
        return & theResult;
      }

    virtual G4ReactionProductVector* 
      Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus) 
      { 
	G4int nntry = 0;
        theTarget.Initialize(theNucleus);
        G4double maxImpact = theTarget.GetOuterRadius();
	G4double startingZ = theSecondaries->operator[](0)->GetPosition().z();
	while(theProjectiles.size()+theEscaped.size()<=1) 
	{
	  G4Pair<G4double, G4double> impact = theTarget.ChooseImpactXandY(maxImpact);
	  G4ThreeVector aNew(impact.first, impact.second, startingZ);
    	  theSecondaries->operator[](0)->SetPosition(aNew);
	  theProjectiles.Initialize(theSecondaries);
	  theEscaped.clear();
	  Stepping();
//	  G4cerr << "Try number "<<++nntry<<" "<<theProjectiles.size()<<" "<<theEscaped.size()<<endl;
	}
	G4int nHoles = theTarget.GetNHoles();
	G4int nParti = theProjectiles.size();
	
        G4int anA= theNucleus->GetMassNumber();;
	G4int aZ = theNucleus->GetCharge();
	G4int numberOfCh=theProjectiles.GetNCharged();
	G4LorentzVector part4Momentum(theProjectiles.GetTotal4Momentum());
	G4ThreeVector exciton3Momentum = part4Momentum.vect();
	G4double exEnergy = theProjectiles.GetTotalKinetic();
	
	// get A,Z of residual.
	const std::vector<G4Nucleon *> & theNuc = theTarget.GetNucleons();
	std::vector<G4Nucleon *>::const_iterator theCurrentNucleon;
	for( theCurrentNucleon = theNuc.begin();
	     theCurrentNucleon!= theNuc.end();
	     theCurrentNucleon++ )
        {
          if((*theCurrentNucleon)->AreYouHit()) 
          {
            anA--;
            aZ -= G4int((*theCurrentNucleon)->GetDefinition()->GetPDGCharge());
            exciton3Momentum -= (*theCurrentNucleon)->Get4Momentum().vect();
            exEnergy+=(*theCurrentNucleon)->GetBindingEnergy();
          }
        }  
	G4cout << "exEnergy = "<<exEnergy<<G4endl;
        exEnergy += G4ParticleTable::GetParticleTable()->GetIonTable()
		    ->GetIonMass( static_cast<G4int>(aZ), static_cast<G4int>(anA) );
	G4LorentzVector exciton4Momentum(exciton3Momentum, exEnergy);
	
	G4Fragment anInitialState;
        anInitialState.SetA(anA);
        anInitialState.SetZ(aZ);
        anInitialState.SetNumberOfParticles(nParti);
        anInitialState.SetNumberOfCharged(numberOfCh);
        anInitialState.SetNumberOfHoles(nHoles);
        anInitialState.SetMomentum(exciton4Momentum);
        const G4Fragment aFragment(anInitialState);
        G4ReactionProductVector * aPreResult = theDeExcitation->DeExcite(aFragment);
	
	G4ReactionProductVector * result = new G4ReactionProductVector;
	GhadParticles::iterator it;
	for(it=theEscaped.begin(); it!=theEscaped.end(); it++)
	{
	  G4ReactionProduct * aNew = new G4ReactionProduct(it->GetDefinition());
	  aNew->SetMomentum(it->GetMom().vect());
	  aNew->SetTotalEnergy(it->GetMom().e());
	  result->push_back(aNew);
	}
	theEscaped.clear();
	for(G4int ipre=0; ipre<aPreResult->size(); ipre++)
	{
//          result->push_back(aPreResult->operator[](ipre));            	  
	}
	delete aPreResult;
	return result;
      }
      
  private:
    void Stepping();

  private:
    G4ParticleChange theResult;
    
    GhadParticles theProjectiles;
    GhadParticles theEscaped;
    GhadNucleus theTarget;
};

#endif
