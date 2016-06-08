// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#include "G4PreCompoundModel.hh"


const G4PreCompoundModel & G4PreCompoundModel::operator=(const G4PreCompoundModel &right)
{
    G4Exception("G4PreCompoundModel::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundModel::operator==(const G4PreCompoundModel &right) const
{
  return false;
}

G4bool G4PreCompoundModel::operator!=(const G4PreCompoundModel &right) const
{
  return true;
}



// Additional Declarations

G4VParticleChange * G4PreCompoundModel::ApplyYourself(const G4Track & thePrimary,
						      G4Nucleus & theNucleus)
{
  theResult.Initialize(thePrimary);
  
  // prepare fragment
  G4Fragment anInitialState;
  G4int anA=G4int(theNucleus.GetN());
  anA += thePrimary.GetDynamicParticle()->GetDefinition()->GetBaryonNumber();
  anInitialState.SetA(anA);
  
  G4int aZ=G4int(theNucleus.GetZ());
  aZ += G4int(thePrimary.GetDynamicParticle()->GetDefinition()->GetPDGCharge());
  anInitialState.SetZ(aZ);
  
  
  // Number of Excitons
  anInitialState.SetNumberOfExcitons(thePrimary.GetDynamicParticle()->GetDefinition()->GetBaryonNumber());
  
  // Number of Charged
  anInitialState.SetNumberOfCharged(thePrimary.GetDynamicParticle()->GetDefinition()->GetPDGCharge());
  
  // Number of Holes 
  anInitialState.SetNumberOfHoles(0);
  
  // pre-compound nucleus energy.
  G4double anEnergy = 0;
  G4double nucleusMass =  G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theNucleus.GetZ(),
                                                                                         theNucleus.GetN());
  anEnergy =  nucleusMass + thePrimary.GetTotalEnergy();
  
  // Momentum
  G4ThreeVector p = thePrimary.GetDynamicParticle()->Get4Momentum().vect();

  // 4-momentum
  G4LorentzVector momentum(p, anEnergy);
  anInitialState.SetMomentum(momentum);
  
  
  
  // call excitation handler
  const G4Fragment aFragment(anInitialState);
  G4ReactionProductVector * result = DeExcite(aFragment);
  
  // fill particle change
  theResult.SetStatusChange(fStopAndKill);
  theResult.SetNumberOfSecondaries(result->length());
  for(G4int i=0; i<result->length(); i++)
    {
      G4DynamicParticle * aNew = 
         new G4DynamicParticle(result->at(i)->GetDefinition(),
                               result->at(i)->GetTotalEnergy(),
                               result->at(i)->GetMomentum());
      delete result->at(i);
      theResult.AddSecondary(aNew);
    }
  delete result;
  
  //return the filled particle change
  return &theResult;
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector* G4PreCompoundModel::DeExcite(const G4Fragment & theInitialState) const
{
  
	G4ReactionProductVector * Result = new G4ReactionProductVector;
  
	// Copy of the initial state 
	G4Fragment aFragment(theInitialState);
  
	// Main loop. It is performed until equilibrium deexcitation.
	for (;;) {
		G4PreCompoundEmission aEmission;
		// Initialize fragment according with the nucleus parameters
		aEmission.Initialize(aFragment);

		// Equilibrium exciton number
		G4double EquilibriumExcitonNumber = 
			sqrt(1.19*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*aFragment.GetA()*
															aFragment.GetExcitationEnergy()+0.5);
															
		// Loop for transitions, it is performed while there are preequilibrium transitions.
		G4bool ThereIsTransition = false;
		do {
			if (aFragment.GetNumberOfExcitons() < EquilibriumExcitonNumber) {
				if (aFragment.GetNumberOfParticles() < 1) {
					aFragment.SetNumberOfHoles(aFragment.GetNumberOfHoles()+1);
					aFragment.SetNumberOfExcitons(aFragment.GetNumberOfExcitons()+2);       
				}
				
				G4double TotalEmissionProbability = aEmission.GetTotalProbability(aFragment);
      	
				// Check if number of excitons is greater than 0
				// else perform equilibrium emission
				if (aFragment.GetNumberOfExcitons() <= 0) {
	  				// Perform Equilibrium Emission
	  				PerformEquilibriumEmission(aFragment,Result);
	  				return Result;
				}
	
				G4PreCompoundTransitions aTransition(aFragment);
	
				// Sum of transition probabilities
				G4double TotalTransitionProbability = aTransition.GetTotalProbability();
	
				// Sum of all probabilities
				G4double TotalProbability = TotalEmissionProbability + TotalTransitionProbability;
	
				// Select subprocess
				if (G4UniformRand() > TotalEmissionProbability/TotalProbability) {
	  				// It will be transition to state with a new number of excitons
	  				ThereIsTransition = true;
					
					// Perform the transition
					aFragment = aTransition.PerformTransition(aFragment);
				} else {
	  				// It will be fragment emission
	  				ThereIsTransition = false;

	  				// Perform the emission and Add emitted fragment to Result
	  				Result->insert(aEmission.PerformEmission(aFragment));
				}
      	} else {
				// Perform Equilibrium Emission
				PerformEquilibriumEmission(aFragment,Result);
				return Result;
      	}
    	} while (ThereIsTransition);   // end of do loop
  	} // end of for (;;) loop
}




void G4PreCompoundModel::PerformEquilibriumEmission(const G4Fragment & aFragment,
						    G4ReactionProductVector * Result) const 
{
 G4ReactionProductVector * theEquilibriumResult;
 theEquilibriumResult = GetExcitationHandler()->BreakItUp(aFragment);
  
  while (theEquilibriumResult->entries() > 0) Result->insert(theEquilibriumResult->removeFirst());
  
  delete theEquilibriumResult;
  return;
}



