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
//
// $Id: G4PreCompoundModel.cc,v 1.16 2001/10/04 20:00:30 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
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
  
  
  // Number of Excited Particles
  anInitialState.SetNumberOfParticles(thePrimary.GetDynamicParticle()->GetDefinition()->GetBaryonNumber());
  
  // Number of Charged Excited Particles
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
  theResult.SetNumberOfSecondaries(result->size());
  for(unsigned int i=0; i<result->size(); i++)
    {
      G4DynamicParticle * aNew = 
	new G4DynamicParticle(result->operator[](i)->GetDefinition(),
			      result->operator[](i)->GetTotalEnergy(),
			      result->operator[](i)->GetMomentum());
      delete result->operator[](i);
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
  
  if (aFragment.GetA() < 5) {
    G4ReactionProduct * theRP = new G4ReactionProduct(G4ParticleTable::GetParticleTable()->
						      GetIon(aFragment.GetZ(),aFragment.GetA(),
							     aFragment.GetExcitationEnergy()));
    theRP->SetMomentum(aFragment.GetMomentum().vect());
    theRP->SetTotalEnergy(aFragment.GetMomentum().e());	  
    Result->push_back(theRP);
    return Result;
  }
  
  G4PreCompoundEmission aEmission(theInitialState);
	
  // Main loop. It is performed until equilibrium deexcitation.
  for (;;) {
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
//  	if (aFragment.GetNumberOfParticles() < 1) {
//  	  aFragment.SetNumberOfHoles(aFragment.GetNumberOfHoles()+1);
//  	  aFragment.SetNumberOfParticles(aFragment.GetNumberOfParticles()+1);       
//  	}
				
	G4double TotalEmissionProbability = aEmission.GetTotalProbability(aFragment);
      	
	// Check if number of excitons is greater than 0
        // else perform equilibrium emission
	if (aFragment.GetNumberOfExcitons() <= 0) {
	  // Perform Equilibrium Emission
#ifdef debug
	  CheckConservation(theInitialState,aFragment,Result);
#endif
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
	  Result->push_back(aEmission.PerformEmission(aFragment));
	}
      } else {
	// Perform Equilibrium Emission
#ifdef debug
	CheckConservation(theInitialState,aFragment,Result);
#endif
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
  
  while (theEquilibriumResult->size() > 0) 
  {
    Result->push_back(*theEquilibriumResult->begin());
    theEquilibriumResult->erase(theEquilibriumResult->begin());
  }
  delete theEquilibriumResult;
  return;
}


#ifdef debug
void G4PreCompoundModel::CheckConservation(const G4Fragment & theInitialState,
					   const G4Fragment & aFragment,
					   G4ReactionProductVector * Result) const
{
  G4double ProductsEnergy = aFragment.GetMomentum().e();
  G4ThreeVector ProductsMomentum = aFragment.GetMomentum();
  G4int ProductsA = G4int(aFragment.GetA());
  G4int ProductsZ = G4int(aFragment.GetZ());
  for (G4int h = 0; h < Result->entries(); h++) {
    ProductsEnergy += Result->at(h)->GetTotalEnergy();
    ProductsMomentum += Result->at(h)->GetMomentum();
    ProductsA += G4int(Result->at(h)->GetDefinition()->GetBaryonNumber());
    ProductsZ += G4int(Result->at(h)->GetDefinition()->GetPDGCharge());
  }

  if (ProductsA != theInitialState.GetA()) {
    G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundModel.cc: Barionic Number Conservation test for just preequilibrium fragments" 
	   << G4endl; 
    G4cout << "Initial A = " << theInitialState.GetA() 
	   << "   Fragments A = " << ProductsA << "   Diference --> " 
	   << theInitialState.GetA() - ProductsA << G4endl;
  }
  if (ProductsZ != theInitialState.GetZ()) {
    G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundModel.cc: Charge Conservation test for just preequilibrium fragments" 
	   << G4endl; 
    G4cout << "Initial Z = " << theInitialState.GetZ() 
	   << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	   << theInitialState.GetZ() - ProductsZ << G4endl;
  }
  if (abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundModel.cc: Energy Conservation test for just preequilibrium fragments" 
	   << G4endl; 
    G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	   << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	   << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
  } 
  if (abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
      abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
      abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PreCompoundModel.cc: Momentum Conservation test for just preequilibrium fragments" 
	   << G4endl; 
    G4cout << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	   << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	   << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV" << G4endl;
  }
  return;
}

#endif
