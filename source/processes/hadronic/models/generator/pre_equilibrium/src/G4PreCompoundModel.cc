// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PreCompoundModel.cc,v 1.6 1999/06/23 09:49:06 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// by V. Lara

#include "G4PreCompoundModel.hh"



G4PreCompoundModel::G4PreCompoundModel(G4ExcitationHandler * const value):
  G4VPreCompoundModel(value)
{
  // neutron 
  theChannels.insert(new G4PreCompoundNeutron());
  // proton
  theChannels.insert(new G4PreCompoundProton());
  // deuterium
  theChannels.insert(new G4PreCompoundDeuteron());
  // triton
  theChannels.insert(new G4PreCompoundTriton());
  // helium3
  theChannels.insert(new G4PreCompoundHe3());
  // alpha
  theChannels.insert(new G4PreCompoundAlpha());
}



G4PreCompoundModel::~G4PreCompoundModel()
{
  theChannels.clearAndDestroy();
}


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
  G4int anA=theNucleus.GetN();
  anA += thePrimary.GetDynamicParticle()->GetDefinition()->GetBaryonNumber();
  anInitialState.SetA(anA);
  
  G4int aZ=theNucleus.GetZ();
  aZ += thePrimary.GetDynamicParticle()->GetDefinition()->GetPDGCharge();
  anInitialState.SetZ(aZ);
  
  
  // Nucleus mass
  //  G4double nucleusMass = 
  //    (theNucleus.GetN()-theNucleus.GetZ())*G4Neutron::Neutron()->GetPDGMass()
  //    + theNucleus.GetZ()*G4Proton::Proton()->GetPDGMass() 
  //    - G4NucleiPropertiesTable::GetBindingEnergy(theNucleus.GetN() , theNucleus.GetZ());
  G4double nucleusMass =  G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theNucleus.GetZ(),
											 theNucleus.GetN());
   
   
  // Excitation Energy
  G4double anEnergy = 0;
  anEnergy =  nucleusMass + thePrimary.GetTotalEnergy();
  // anEnergy += -aZ*G4Proton::Proton()->GetPDGMass() 
  //   - (anA-aZ)*G4Neutron::Neutron()->GetPDGMass()
  //   -G4NucleiPropertiesTable::GetBindingEnergy(anA,aZ);
  anEnergy -= G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ,anA);
  //  anInitialState.SetExcitationEnergy(anEnergy);
  
  // Number of Excitons
  anInitialState.SetNumberOfExcitons(thePrimary.GetDynamicParticle()->GetDefinition()->GetBaryonNumber());
  
  // Number of Charged
  anInitialState.SetNumberOfCharged(thePrimary.GetDynamicParticle()->GetDefinition()->GetPDGCharge());
  
  // Number of Holes 
  anInitialState.SetNumberOfHoles(0);
  
  // Momentum
  G4ThreeVector p = thePrimary.GetDynamicParticle()->Get4Momentum().vect();
  G4LorentzVector momentum(p, sqrt(p.mag2()+(anEnergy+nucleusMass) * (anEnergy+nucleusMass)) );
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
  // result = GetExcitationHandler()->BreakItUp(aFragment);
  
  G4Fragment aFragment(theInitialState);
  
  
  // Main loop. It is performed until equilibrium deexcitation.
  for (;;) {
    
    // Compute atomic numbers and charges for rest nuclei
    for (G4int i = 0; i < NumberOfPossibleFragments; i++) {
      theChannels(i)->Init(aFragment);
    }
    
    // Equilibrium exciton number
    G4double EquilibriumExcitonNumber = sqrt(1.19*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*
					     aFragment.GetA()*aFragment.GetExcitationEnergy()/MeV+0.5);
    
    // Loop for transitions, it is performed while there are preequilibrium transitions.
    G4bool ThereIsTransition = false;
    do {
      if (aFragment.GetNumberOfExcitons() < EquilibriumExcitonNumber) {
	if (aFragment.GetNumberOfParticles() < 1) {
	  aFragment.SetNumberOfHoles(aFragment.GetNumberOfHoles()+1);
	  aFragment.SetNumberOfExcitons(aFragment.GetNumberOfExcitons()+2);       
	}
	
	G4double TotalEmissionProbability = 0.0;
	G4int i;
        for (i = 0; i < NumberOfPossibleFragments; i++) {
	  theChannels(i)->CalcExcitonLevelDensityRatios(aFragment.GetNumberOfParticles()+
							aFragment.GetNumberOfHoles(),
							aFragment.GetNumberOfParticles());	
	  theChannels(i)->CalcCondensationProbability(aFragment.GetA());
	  // Calculate emission probailities
	  if (aFragment.GetNumberOfParticles() <= theChannels(i)->GetA()-0.01)
	    // if number of particles less than a fragment atomic number 
	    // set probability to emit a fragment 0
	    theChannels(i)->SetEmissionProbability(0.0);
	  else if (aFragment.GetNumberOfExcitons() <= theChannels(i)->GetA()+0.01 && 
		   aFragment.GetNumberOfExcitons() != 1)
	    theChannels(i)->SetEmissionProbability(0.0);    	  
	  else if (aFragment.GetNumberOfCharged() <= theChannels(i)->GetZ()-0.01) 
	    // if number of charged particles (protons) is less than charge of fragment
	    // set probability to emit a fragment 0
	    theChannels(i)->SetEmissionProbability(0.0);
	  else if (theChannels(i)->GetMaximalKineticEnergy() <= 0.0) 
	    // if the energy threshold for emitted fragment is less or equal 0
	    // set probability to emit a fragment 0
	    theChannels(i)->SetEmissionProbability(0.0);
	  else
	    // Compute total (integrated over kinetic energy) emission 
	    // probability of a fragment and
	    // Summing channel emission probabilities
	    TotalEmissionProbability += theChannels(i)->CalcEmissionProbability(aFragment);
	}
      	
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
	  G4int deltaN = aTransition.GetDeltaNExciton();
	  aFragment.SetNumberOfExcitons(aFragment.GetNumberOfExcitons()+deltaN);
	  aFragment.SetNumberOfHoles(aFragment.GetNumberOfHoles()+deltaN/2);
	  // With weight Z/A, number of charged particles is decreased on +1
	  if ((deltaN > 0 || aFragment.GetNumberOfCharged() > 0) &&
              (G4UniformRand() <= aFragment.GetZ()/aFragment.GetA()))
	    aFragment.SetNumberOfCharged(aFragment.GetNumberOfCharged()+deltaN/2);
	} else {
	  // It will be fragment emission
	  ThereIsTransition = false;
	  G4double * running = new G4double[NumberOfPossibleFragments];
	  running[0] = theChannels(0)->GetEmissionProbability();
	  for (i = 1; i < NumberOfPossibleFragments; i++) 
	    running[i]=running[i-1]+theChannels(i)->GetEmissionProbability();
	  
	  // Choose an emission channel
	  G4double ChoosedChannel = G4UniformRand()*TotalEmissionProbability;
	  G4int aChannel = -1;
	  for (i = 0; i < NumberOfPossibleFragments; i++) {
	    if (ChoosedChannel <= running[i]) {
	      aChannel = i;
	      break;
	    }
	  }
	  delete [] running;
	  
	  // Compute Kinetic Energy of emitted fragment
	  G4double KineticEnergyOfEmittedFragment = 
	    theChannels(aChannel)->GetKineticEnergy(aFragment);
	  

	  // G4cout << "Kinetic energy of Emitted fragment " << KineticEnergyOfEmittedFragment << endl;
	  
	  // Sample Fermi momentum of emitted fragment
	  static const G4double FermiMaxMom = 250.0; // MeV
	  G4ThreeVector FermiMomentum(IsotropicRandom3Vector(FermiMaxMom*pow(G4UniformRand(),1./3.)));

	  G4ThreeVector P12(FermiMomentum + 
			    // (1/#Particles before emission)
			    ( 1.0/G4double(aFragment.GetNumberOfParticles()) )*
			    aFragment.GetMomentum().vect()
			    );
	  


	  G4double p = sqrt(KineticEnergyOfEmittedFragment*(KineticEnergyOfEmittedFragment+
							    2.0*theChannels(aChannel)->GetNuclearMass()));
	  G4ParticleMomentum momentum;


	  if (aFragment.GetMomentum().boostVector().mag2() > 1.e-7) {
	    // sample a non-isotropic random vector
	    G4double CosTheta = sqrt(G4UniformRand());
	    G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
	    G4double Phi = twopi*G4UniformRand();
	    momentum = G4ParticleMomentum(p*cos(Phi)*SinTheta,
					  p*sin(Phi)*SinTheta,
					  p*CosTheta);
	    momentum = RotateMomentum(P12,aFragment.GetMomentum().boostVector(),momentum);
	    
	  } else {
	    
	    momentum = IsotropicRandom3Vector(p);

	  } 

	  
	  G4LorentzVector EmittedMomentum(momentum,
					  sqrt(momentum.mag2()+
					       theChannels(aChannel)->GetNuclearMass() *
					       theChannels(aChannel)->GetNuclearMass() )
					  );

	  // Excitation energy
	  //                    check that Excitation energy is > 0
	  G4double CheckU = theChannels(aChannel)->GetMaximalKineticEnergy() - 
	    KineticEnergyOfEmittedFragment + 
	    theChannels(aChannel)->GetCoulombBarrier();
	  if (CheckU < 0.0)
	    G4Exception("G4PreCompoundModel::DeExcite: Excitation energy less than 0! ");
	  

	  // Update nucleus parameters
	  // Number of excitons
	  aFragment.SetNumberOfExcitons(aFragment.GetNumberOfExcitons()-
					G4int(theChannels(aChannel)->GetA()));
	  // Number of charges
	  aFragment.SetNumberOfCharged(aFragment.GetNumberOfCharged()-
				       G4int(theChannels(aChannel)->GetZ()));

	  // Atomic number
	  aFragment.SetA(theChannels(aChannel)->GetRestA());
	  
	  // Charge
	  aFragment.SetZ(theChannels(aChannel)->GetRestZ());
	  
	  //	  aFragment.SetExcitationEnergy(CheckU);
 	  G4LorentzVector RestMomentum(-momentum,
 				       sqrt(momentum.mag2()+ 
 					    (theChannels(aChannel)->GetRestNuclearMass()+CheckU)*
 					    (theChannels(aChannel)->GetRestNuclearMass()+CheckU)
					    )
				       );
	  
	  // Perform Lorentz boosts
	  EmittedMomentum.boost(aFragment.GetMomentum().boostVector());
	  RestMomentum.boost(aFragment.GetMomentum().boostVector());
	  
	  // Update nucleus momentum
	  aFragment.SetMomentum(RestMomentum);
	  
	  // Set emitted fragment momentum
	  theChannels(aChannel)->SetMomentum(EmittedMomentum);
	  

	  // Add emitted fragment to Result
	  G4DynamicParticle MyDP = theChannels(aChannel)->GetDynamicParticle();
	  G4ReactionProduct * theNew = new G4ReactionProduct(MyDP.GetDefinition());
          theNew->SetMomentum(MyDP.GetMomentum());
          theNew->SetTotalEnergy(MyDP.Get4Momentum().e());
//	  delete MyDP;
	  Result->insert(theNew);
	}
      } else {
	// Perform Equilibrium Emission
	PerformEquilibriumEmission(aFragment,Result);
	return Result;
      }
    } while (ThereIsTransition);   // end of do loop
  } // end of for (;;) loop
}



G4ThreeVector G4PreCompoundModel::IsotropicRandom3Vector(G4double Magnitude) const
  // Create a unit vector with a random direction isotropically distributed
{

  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ThreeVector Vector(Magnitude*cos(Phi)*SinTheta,
                       Magnitude*sin(Phi)*SinTheta,
                       Magnitude*CosTheta);

  return Vector;
		
}


void G4PreCompoundModel::PerformEquilibriumEmission(const G4Fragment & aFragment,
						    G4ReactionProductVector * Result) const 
{

  
 G4ReactionProductVector * theEquilibriumResult;
 theEquilibriumResult = GetExcitationHandler()->BreakItUp(aFragment);
  
  while (theEquilibriumResult->entries() > 0) 
    Result->insert(theEquilibriumResult->removeFirst());
  
  delete theEquilibriumResult;
  
  return;
}



G4ParticleMomentum G4PreCompoundModel::RotateMomentum(G4ParticleMomentum Pa,
						      G4ParticleMomentum V,
						      G4ParticleMomentum P) const 
{
  G4ParticleMomentum U = Pa.unit();
  
  G4double Alpha1 = U * V;
  
  G4double Alpha2 = sqrt(V.mag2() - Alpha1*Alpha1);

  G4ThreeVector N = (1./Alpha2)*U.cross(V);

  G4ParticleMomentum RotatedMomentum(
				     ( (V.x() - Alpha1*U.x())/Alpha2 ) * P.x() + N.x() * P.y() + U.x() * P.z(),
				     ( (V.y() - Alpha1*U.y())/Alpha2 ) * P.x() + N.y() * P.y() + U.y() * P.z(),
				     ( (V.z() - Alpha1*U.z())/Alpha2 ) * P.x() + N.z() * P.y() + U.z() * P.z()
				     );
  return RotatedMomentum;
}
