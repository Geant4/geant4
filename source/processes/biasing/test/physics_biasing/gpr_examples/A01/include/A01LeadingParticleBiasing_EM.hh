
#ifndef A01LEADINGPARTICLEBIASING_EM_HH
#define A01LEADINGPARTICLEBIASING_EM_HH


namespace A01LeadingParticleBiasing_EM {

  G4VParticleChange* SimpleEM_Conv(G4GPRProcessWrappers::G4GPRDiscreteDoIt& original,
				   const G4Track& track, const G4Step& step)
  {

    G4VParticleChange* particleChange = original(track, step);
    
    G4int nSecondaries = particleChange->GetNumberOfSecondaries();

    G4cout<<"jane leading particle biasing nsecondaries "<<track.GetDefinition()->GetParticleName()<<" "<<original.GetIdentifier()<<" "<<nSecondaries<<" "<<G4endl;

    assert (particleChange->GetTrackStatus() == fStopAndKill);
    assert (nSecondaries == 2);

    G4Track* secondary0 = particleChange->GetSecondary(0);
    G4Track* secondary1 = particleChange->GetSecondary(1);

    G4Track* lowerEnergy(0);
    G4Track* higherEnergy(0);

    if (secondary0->GetKineticEnergy() < secondary1->GetKineticEnergy()) {
      lowerEnergy = secondary0;
      higherEnergy = secondary1;
    }
    else {
      lowerEnergy = secondary1;
      higherEnergy = secondary0;
    }
    
    G4double eLow = lowerEnergy->GetKineticEnergy();
    G4double eHigh = higherEnergy->GetKineticEnergy();
    G4double eTot = eHigh+eLow;

    G4cout<<"jane em lpb "<<eLow<<" "<<eHigh<<G4endl;

    G4double fraction = eLow/eTot;
      
    G4double rand = G4UniformRand();
    G4double weight(0);
    G4Track* result(0);

    if (rand < fraction) {
      result = new G4Track(*lowerEnergy);
      weight = eTot/eLow;
      G4cout<<"jane lpb selected lower"<<G4endl;
    }
    else {
      result = new G4Track(*higherEnergy);
      weight = eTot/eHigh;
      G4cout<<"jane lpb selected higher"<<G4endl;
    }
    
    result->SetWeight(weight);
    particleChange->SetNumberOfSecondaries(1);
    particleChange->SetSecondaryWeightByProcess(true);
    particleChange->AddSecondary(result);

    nSecondaries = particleChange->GetNumberOfSecondaries();
    G4cout<<"jane new secondaires "<<nSecondaries<<G4endl;
    return particleChange;

  }

  G4VParticleChange* SimpleEM(G4GPRProcessWrappers::G4GPRDiscreteDoIt& original,
			      const G4Track& track, const G4Step& step)
  {    
    G4cout<<"jane lead particle selection - original functor id "<<original.GetIdentifier()<<" "<<track.GetVolume()->GetName()<<G4endl;

    G4VParticleChange* particleChange = original(track, step);

    G4cout<<"jane in lead particle  "<<particleChange->GetNumberOfSecondaries()<<" "<<track.GetVolume()->GetName()<<G4endl;
    
    G4int nSecondaries = particleChange->GetNumberOfSecondaries();
    if ((nSecondaries == 0) || (nSecondaries == 1 && (particleChange->GetTrackStatus() == fStopAndKill))) return particleChange;

    assert (nSecondaries == 1);

    //Dummy step to extract info - jane fixme - better way to do this ? too heavy.
    G4Step dummyStep;
    G4Track* tmpTrack = step.GetTrack();
    dummyStep.InitializeStep(tmpTrack);
    
    particleChange->UpdateStepForPostStep(&dummyStep);
    
    G4double energyParent = dummyStep.GetPostStepPoint()->GetKineticEnergy();
    G4double energyDaughter = particleChange->GetSecondary(0)->GetKineticEnergy();

    G4double rand = G4UniformRand();
    G4double eTot = energyParent + energyDaughter;

    if (energyParent < energyDaughter) {
      G4double fraction = energyParent/eTot;

      if (rand < fraction) {
	// Keep the parent, kill the daughter

	G4double weight = eTot/energyParent;
	particleChange->SetNumberOfSecondaries(0);
	// jane fixme - does this work yet ?
	particleChange->ProposeParentWeight(weight);
      }
      else {
	// Kill the parent, keep the daughter
	particleChange->ProposeTrackStatus(fStopAndKill);

	G4double weight = eTot/energyDaughter;
	particleChange->GetSecondary(0)->SetWeight(weight);
      }
    }
    else {
      G4double fraction = energyDaughter/eTot;

      if (rand < fraction) {
	// Keep the daughter, kill the parent
	G4double weight = eTot/energyDaughter;
	particleChange->ProposeTrackStatus(fStopAndKill);
	particleChange->GetSecondary(0)->SetWeight(weight);
      }
      else {
	// Keep the parent, kill the daughter
	G4double weight = eTot/energyParent;
	particleChange->SetNumberOfSecondaries(0);
	// jane fixme - does this work yet ?
	particleChange->ProposeParentWeight(weight);
      }
    }
  
    return particleChange;
  }
}

#endif
