#ifndef A01LEADINGPARTICLEBIASING_Hadronic_HH
#define A01LEADINGPARTICLEBIASING_Hadronic_HH

#include "G4PionZero.hh"
namespace A01LeadingParticleBiasing_Hadronic {

  G4VParticleChange* Biasing(G4GPRProcessWrappers::G4GPRDiscreteDoIt& original, const G4Track& track, const G4Step& step)
  {
    // jane fixme - Different logic from G4HadLeadBias since G4HadronicProcess 
    // kills the incoming and adds it as a secondary

    G4VParticleChange* result(0);
    
    result = original(track, step);
    G4cout << "jane entering hadleadbias "<<track.GetTrackID()<<" "<<result->GetNumberOfSecondaries()<<G4endl;

    // G4cerr << "bias enter"<<G4endl;
    G4int nMeson(0), nBaryon(0), npi0(0), ngamma(0), nLepton(0);
    G4int i(0);
    G4int maxE = -1;
    G4double emax = 0;

    G4int nSecondaries = result->GetNumberOfSecondaries();
    G4int nSecondaryLoop(nSecondaries-1);

    if (0 == nSecondaries) return result;
    
    G4Track* incoming = result->GetSecondary(nSecondaries-1);

    emax = incoming->GetDynamicParticle()->GetKineticEnergy();

    G4cout<<"jane emax "<<emax<<" "<<G4endl;

    for(i=0;i<nSecondaryLoop;i++) {
      G4cout<<"jane secondary energy "<<result->GetSecondary(i)->GetKineticEnergy()<<G4endl;
      if(result->GetSecondary(i)->GetDynamicParticle()->GetKineticEnergy()>emax) {
	maxE = i;
        emax = result->GetSecondary(i)->GetDynamicParticle()->GetKineticEnergy();
      }
    }
    G4cout << "max energy "<<result->GetNumberOfSecondaries()<<" "<<maxE<<" "<<emax<<G4endl;

    //G4cout <<"loop1"<<G4endl;
    for(i=0; i<nSecondaryLoop; i++) {
      const G4DynamicParticle* aSecTrack = result->GetSecondary(i)->GetDynamicParticle();
      if(i==maxE) {
      }
      else if(aSecTrack->GetDefinition()->GetBaryonNumber()!=0) {
	nBaryon++;
      }
      else if(aSecTrack->GetDefinition()->GetLeptonNumber()!=0) {
        nLepton++;
      }
      else if(aSecTrack->GetDefinition()==G4Gamma::Gamma()) {
        ngamma++;
      }
      else if(aSecTrack->GetDefinition()==G4PionZero::PionZero()) {
	npi0++;
      }
      else {
	nMeson++;
      }
    }
    G4cout << "BiasDebug 1 = "<<result->GetNumberOfSecondaries()<<" "
           <<nMeson<<" "<< nBaryon<<" "<< npi0<<" "<< ngamma<<" "<< nLepton<<G4endl;
    G4double mesonWeight = nMeson;
    G4double baryonWeight = nBaryon;
    G4double gammaWeight = ngamma;
    G4double npi0Weight = npi0;
    G4double leptonWeight = nLepton;
    G4int randomMeson = static_cast<G4int>((nMeson+1)*G4UniformRand());
    G4int randomBaryon = static_cast<G4int>((nBaryon+1)*G4UniformRand());
    G4int randomGamma = static_cast<G4int>((ngamma+1)*G4UniformRand());
    G4int randomPi0 = static_cast<G4int>((npi0+1)*G4UniformRand());
    G4int randomLepton = static_cast<G4int>((nLepton+1)*G4UniformRand());
    G4cout<<"jane randdebug "<< randomMeson<<" "<<randomBaryon<<" "<<randomGamma<<" "<<randomPi0<<" "<<randomLepton<<G4endl;
    std::vector<G4Track *> buffer;
    G4int cMeson(0), cBaryon(0), cpi0(0), cgamma(0), cLepton(0);
    for(i=0; i<nSecondaryLoop; i++) {
      G4bool aCatch = false;
      G4double weight = 1;
      G4Track * aSecTrack = new G4Track(*(result->GetSecondary(i)));
      if(i==maxE) {
	aCatch = true;
	weight = 1;
      }
      else if(aSecTrack->GetDynamicParticle()->GetDefinition()->GetBaryonNumber()!=0) {
	if(++cBaryon==randomBaryon) {
	  aCatch = true;
	  weight = baryonWeight;
        }
      }
      else if(aSecTrack->GetDynamicParticle()->GetDefinition()->GetLeptonNumber()!=0) {
	if(++cLepton==randomLepton) {
	  aCatch = true;
	  weight = leptonWeight;
        }
      }
      else if(aSecTrack->GetDynamicParticle()->GetDefinition()==G4Gamma::Gamma()) {
	if(++cgamma==randomGamma) {
          aCatch = true;
          weight = gammaWeight;
	}
      }
      else if(aSecTrack->GetDynamicParticle()->GetDefinition()==G4PionZero::PionZero()) {
	if(++cpi0==randomPi0) {
	  aCatch = true;
	  weight = npi0Weight;
	}
      }
      else {
	if(++cMeson==randomMeson) {
          aCatch = true;
          weight = mesonWeight;
        }
      }
      if(aCatch) {
        buffer.push_back(aSecTrack);
        aSecTrack->SetWeight(aSecTrack->GetWeight()*weight);
      }
      else {
	delete aSecTrack;
      }
    }
    buffer.push_back(new G4Track(*incoming));

    result->SetNumberOfSecondaries(0);
    result->SetNumberOfSecondaries(buffer.size());
    // G4cerr << "pre"<<G4endl;
    for(i=0;i<static_cast<G4int>(buffer.size());i++) {
      G4cout<<"jane lpbnew result: "<<buffer[i]->GetKineticEnergy()<<G4endl;
      result->AddSecondary(buffer[i]);
    }
    // G4cerr << "bias exit"<<G4endl;
    G4cout << "jane leaving hadleadbias "<<result->GetNumberOfSecondaries()<<G4endl;
    
    return result;
  }
  
}

#endif
