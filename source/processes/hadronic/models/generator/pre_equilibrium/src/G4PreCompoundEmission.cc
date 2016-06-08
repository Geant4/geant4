#include "G4PreCompoundEmission.hh"

const G4PreCompoundEmission & G4PreCompoundEmission::operator=(const G4PreCompoundEmission &right)
{
    G4Exception("G4PreCompoundEmission::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundEmission::operator==(const G4PreCompoundEmission &right) const
{
  return false;
}

G4bool G4PreCompoundEmission::operator!=(const G4PreCompoundEmission &right) const
{
  return true;
}


G4ReactionProduct * G4PreCompoundEmission::PerformEmission(G4Fragment & aFragment)
{
	// Choose a Fragment for emission
	G4VPreCompoundFragment * theFragment = theFragmentsVector.ChooseFragment();

	// Kinetic Energy of emitted fragment
	G4double KineticEnergyOfEmittedFragment = theFragment->GetKineticEnergy(aFragment);
	  	  
	// Sample the Fermi momentum of emitted fragment
	static const G4double FermiMaxMom = 250.0*MeV;
	G4ThreeVector FermiMomentum(IsotropicRandom3Vector(FermiMaxMom*pow(G4UniformRand(),1./3.)));

	// Get the fragment momentum
	G4ThreeVector P12(aFragment.GetMomentum().vect());	
	// Share the fragment momentum between the particles system
	P12 *= 1.0/G4double(aFragment.GetNumberOfParticles());
	// Add the Fermi momentum
	P12 += FermiMomentum;
			
	// Calculate the momentum magnitude of emitted fragment 	
	G4double EmittedMass = theFragment->GetNuclearMass();
	G4double p = sqrt(KineticEnergyOfEmittedFragment*(KineticEnergyOfEmittedFragment+2.0*EmittedMass));
	// And sample a direction for it
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
	
	// Now we can calculate the four momentum 
	G4LorentzVector EmittedMomentum(momentum,sqrt(momentum.mag2()+EmittedMass*EmittedMass));
	

	// Excitation energy
	G4double anU = theFragment->GetMaximalKineticEnergy() - KineticEnergyOfEmittedFragment + 
						  theFragment->GetCoulombBarrier();
	
	// check that Excitation energy is > 0
	if (anU < 0.0) G4Exception("G4PreCompoundModel::DeExcite: Excitation energy less than 0!");

	// Update nucleus parameters
	// Number of excitons
	aFragment.SetNumberOfExcitons(aFragment.GetNumberOfExcitons()-
															G4int(theFragment->GetA()));
	// Number of charges
	aFragment.SetNumberOfCharged(aFragment.GetNumberOfCharged()-
				 	      								G4int(theFragment->GetZ()));

	// Atomic number
	aFragment.SetA(theFragment->GetRestA());
	  
	// Charge
	aFragment.SetZ(theFragment->GetRestZ());


	// Calculate the residual Fragment momentum
	G4double ResidualMass = theFragment->GetRestNuclearMass()+anU;
	G4LorentzVector RestMomentum(-momentum,sqrt(momentum.mag2()+ ResidualMass*ResidualMass));	  

	// Perform Lorentz boosts
   EmittedMomentum.boost(aFragment.GetMomentum().boostVector());
	RestMomentum.boost(aFragment.GetMomentum().boostVector());


	// Update nucleus momentum
	aFragment.SetMomentum(RestMomentum);
	
	// Set emitted fragment momentum
	theFragment->SetMomentum(EmittedMomentum);	

	G4DynamicParticle MyDP = theFragment->GetDynamicParticle();
	G4ReactionProduct * theNew = new G4ReactionProduct(MyDP.GetDefinition());
	theNew->SetMomentum(MyDP.GetMomentum());
	theNew->SetTotalEnergy(MyDP.Get4Momentum().e());

	return theNew;
}


G4ThreeVector G4PreCompoundEmission::IsotropicRandom3Vector(G4double Magnitude) const
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


G4ParticleMomentum G4PreCompoundEmission::RotateMomentum(G4ParticleMomentum Pa,
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
