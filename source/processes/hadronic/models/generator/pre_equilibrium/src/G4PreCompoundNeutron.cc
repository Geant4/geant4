// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#include "G4PreCompoundNeutron.hh"


G4double G4PreCompoundNeutron::ProbabilityDistributionFunction(const G4double & eKin,
																				  const G4Fragment & aFragment)
{
	const G4double r0 = 1.5*fermi;
	// g = 0.595*a*A
	const G4double g = 0.595*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*GetRestA(); 
	G4double Alpha = 0.76+2.2/pow(GetRestA(),1.0/3.0); 
	G4double Beta = (2.12/pow(GetRestA(),2.0/3.0)-0.05)*MeV/Alpha;
	
	G4double Probability = 2.0/(pi*hbarc*hbarc*hbarc) * GetReducedMass() * Alpha * 
										r0 * r0 * pow(GetRestA(),2.0/3.0) *
										GetExcitonLevelDensityRatio()/(g*aFragment.GetExcitationEnergy()) *
										pow((1.0 - (eKin+GetBindingEnergy())/aFragment.GetExcitationEnergy()),
	 									(aFragment.GetNumberOfExcitons()-2.0))*(eKin + Beta); 

// 	G4double Probability = 0.000234*r0*r0*pow(GetRestA(),2.0/3.0)*Alpha*GetExcitonLevelDensityRatio()/
//     (SingleParticleLevelDensity*aFragment.GetExcitationEnergy())*
//     pow((1.0 - (eKin+GetBindingEnergy())/aFragment.GetExcitationEnergy()),
// 	 		(aFragment.GetNumberOfExcitons()-2.0))*(eKin + Beta); 
			
	return Probability;
}


G4double G4PreCompoundNeutron::GetKineticEnergy(const G4Fragment & aFragment)
{
	G4double Beta = (2.12/pow(GetRestA(),2.0/3.0)-0.05)*MeV/(0.76+2.2/pow(GetRestA(),1.0/3.0));
					 
	G4double T = aFragment.GetNumberOfParticles() + aFragment.GetNumberOfHoles() - GetA() - 1.0;
	G4double R2 = GetMaximalKineticEnergy();
	G4double R1 = R2 + GetCoulombBarrier();

	G4double E = 0.0;

	if (T <= -0.1) {
		 E = R1; 
	} else if (T <= 0.1) {
		E = -Beta + sqrt(Beta*Beta + (G4UniformRand()*(R2*R2 + 2.0*Beta*R2)));
	} else {
		G4double E1 = (R1 - Beta*T)/(T + 1.0);
		G4double T3 = 0.0;
		do {
			E = GetCoulombBarrier()+G4UniformRand()*R2;
			G4double T1 = (E + Beta)/(E1 + Beta);
			G4double T2 = (R1 - E)/(R1 - E1);
			T3 = T1*pow(T2,T);
		} while (G4UniformRand() > T3);
	}
	return E;
}
