// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara
// Corrections by V. Krylov

#include "G4VPreCompoundIon.hh"  

G4double G4VPreCompoundIon::
ProbabilityDistributionFunction(const G4double & eKin,
				const G4Fragment & aFragment)
{
  const G4double r0 = 1.5; // fm
  const G4double SingleParticleLevelDensity = 
    0.595*G4PreCompoundParameters::GetAddress()->GetLevelDensity(); // AC

  G4double R0J = 1.1;
  G4double exEnergy = aFragment.GetExcitationEnergy()/MeV;
  G4double probA = GetCondensationProbability()*R0J*0.104/
    (r0*pow(GetRestA(),1.0/3.0)*sqrt(GetA()*exEnergy));
  G4double probB =     GetExcitonLevelDensityRatio()*
    ( (eKin-GetCoulombBarrier())/exEnergy );
  G4double ratio = (eKin+GetBindingEnergy() )/exEnergy;
  G4double exponent = GetRestA()-1.5;
  if ( exponent>100. && ratio<1. ) return 0.;
  G4double probC = pow( ratio, exponent );
  G4double probD = pow( 1.0 - ratio,
			aFragment.GetNumberOfExcitons()-GetA()-1.0 ) ;
	
//   return GetCondensationProbability()*R0J*0.104/                 
//          (r0*pow(GetRestA(),1.0/3.0)*sqrt(GetA()*exEnergy))*   
//          GetExcitonLevelDensityRatio()*
//          ( (eKin-GetCoulombBarrier())/exEnergy )*
//          pow( ( (eKin+GetBindingEnergy() )/exEnergy), GetRestA()-1.5)*
//          pow(1.0 - (eKin + GetBindingEnergy())/exEnergy  ,
// 	          aFragment.GetNumberOfExcitons()-GetA()-1.0 ) ;

  G4double prob = probA*probB*probC*probD;
  if (prob < 1.e-100) return 0.;
  else return prob;

  // Corrections in return statemet by V. Krylov:
  //    - GetA() and GetRestA() were intechanged

}


G4double G4VPreCompoundIon::GetKineticEnergy(const G4Fragment & aFragment)
{
  G4double DJ = - GetCoulombBarrier();

  G4double T = aFragment.GetNumberOfParticles() + aFragment.GetNumberOfHoles() - GetA() - 1.0;
  G4double R2 = GetMaximalKineticEnergy();
  G4double R1 = R2 + GetCoulombBarrier();
	

  if (T <= -0.1) return R1;
  else if (T <= 0.1) {
    G4double E1 = R1;
    G4double E = 0.0;
    G4double T3 = 0.0;
    do {
      G4double PJ1 = GetA() - 1.5;
      G4double AbsBindingE = abs(GetBindingEnergy());
      if (GetBindingEnergy() <= 0.0 && AbsBindingE > GetCoulombBarrier()) 
	E = AbsBindingE + G4UniformRand()*(aFragment.GetExcitationEnergy()/MeV);
      else
	E = GetCoulombBarrier() + G4UniformRand()*R2;	
      T3 = pow((E+GetBindingEnergy())/(E1+GetBindingEnergy()),PJ1)*
	((E+DJ)/(E1+DJ));
    } while (G4UniformRand() > T3);
    return E;
  } else {
    G4double PJ1 = GetA() - 1.5;
    G4double ES = (aFragment.GetExcitationEnergy()/MeV)*(GetA()-0.5)+
      ((aFragment.GetExcitationEnergy()/MeV)-R2)*(aFragment.GetNumberOfParticles()+
						  aFragment.GetNumberOfHoles()-2.5);
    G4double E1 = (ES + sqrt(ES*ES-((aFragment.GetExcitationEnergy()/MeV)-R2)*(GetA()-1.5)*
			     (aFragment.GetNumberOfParticles()+aFragment.GetNumberOfHoles()-1.5)*
			     4.0*(aFragment.GetExcitationEnergy()/MeV)))/
      ((aFragment.GetNumberOfParticles()+aFragment.GetNumberOfHoles()-1.5)*2.0)
      - (aFragment.GetExcitationEnergy()/MeV) + R1;
    G4double E = 0.0;
    G4double T3 = 0.0;
    do {
      if (GetBindingEnergy() <= 0.0 && abs(GetBindingEnergy()) > GetCoulombBarrier()) 
	E = abs(GetBindingEnergy()) + G4UniformRand()*(aFragment.GetExcitationEnergy()/MeV);
      else 
	E = GetCoulombBarrier()+G4UniformRand()*R2;
      T3 = (pow((E+GetBindingEnergy())/(E1+GetBindingEnergy()),PJ1)*
	    ((E+DJ)/(E1+DJ))) * pow((R1-E)/(R1-E1),T);
    } while (G4UniformRand() > T3);
    return E;
  }
}
