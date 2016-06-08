// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#include "G4PreCompoundProton.hh"


G4double G4PreCompoundProton::ProbabilityDistributionFunction(const G4double & eKin,
																				 const G4Fragment & aFragment)
{
  const G4double r0 = 1.5*fermi;
  // g = 0.595*a*A;
  const G4double g = 0.595*G4PreCompoundParameters::GetAddress()->GetLevelDensity()*GetRestA();
//  G4double R0J=1.2;

	G4double aZ = G4double(GetRestZ());
	G4double C = 0.0;
	if (aZ >= 70) {
		C = 0.10;
	} else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
	}                                         


	G4double Probability = 2.0/(pi*hbarc*hbarc*hbarc) * GetReducedMass() * 
										(1.0 + C) * r0 * r0 * pow(GetRestA(),2.0/3.0) *
                            GetExcitonLevelDensityRatio()/(g*aFragment.GetExcitationEnergy()) *
										pow(1.0 - (eKin+GetBindingEnergy())/aFragment.GetExcitationEnergy(),
										(aFragment.GetNumberOfExcitons()-2.0))*
										(eKin - GetCoulombBarrier());
										
//   G4double Probability = 0.000234*r0*r0*pow(GetRestA(),2.0/3.0)*R0J*GetExcitonLevelDensityRatio()/
//     (SingleParticleLevelDensity*aFragment.GetExcitationEnergy())*
//     pow(1.0 - (eKin+GetBindingEnergy())/aFragment.GetExcitationEnergy(),
// 	(aFragment.GetNumberOfExcitons()-2.0))*(eKin - GetCoulombBarrier());
	
	return Probability;
}




G4double G4PreCompoundProton::GetKineticEnergy(const G4Fragment & aFragment)
{
	G4double DJ = - GetCoulombBarrier();

	G4double T = aFragment.GetNumberOfParticles() + aFragment.GetNumberOfHoles() - GetA() - 1.0;
	G4double R2 = GetMaximalKineticEnergy();
	G4double R1 = R2 + GetCoulombBarrier();
	
	G4double E = 0.0; 

	if (T <= -0.1) {
  		E = R1;
	} else if (T <= 0.1) {
		E = sqrt(G4UniformRand())*R2 + GetCoulombBarrier();
	} else {
		G4double E1 = (R1 - DJ*T)/(T + 1.0);
		G4double T3 = 0.0;
		do {
			E = GetCoulombBarrier() + G4UniformRand()*R2;
			G4double T1 = (E + DJ)/(E1 + DJ);
			G4double T2 = (R1 - E)/(R1 - E1);
			T3 = T1*pow(T2,T);
		} while (G4UniformRand() > T3);
	}
	return E;
}


G4double G4PreCompoundProton::GetBarrierPenetrationFactor(const G4double aZ) const
{
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.80;
	} else {
		K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
	}                                  

	return K;
}




