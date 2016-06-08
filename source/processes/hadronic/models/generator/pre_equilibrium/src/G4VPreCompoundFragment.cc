// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara
// corrections by V. Krylov
 
#include "G4VPreCompoundFragment.hh"

//G4VPreCompoundFragment::G4VPreCompoundFragment(): theA(0.0),theZ(0.0),
//  theRestNucleusA(0.0),theRestNucleusZ(0.0),CoulombBarrier(0.0),MaximalKineticEnergy(-1.0),
//  ExcitonLevelDensityRatio(0.0),EmissionProbability(0.0),CondensationProbability(0.0),
//  Momentum(0.0,0.0,0.0,0.0)
//{
//}


G4VPreCompoundFragment::G4VPreCompoundFragment
(const G4VPreCompoundFragment & right)
{
  theA = right.theA;
  theZ = right.theZ;
  theRestNucleusA = right.theRestNucleusA;
  theRestNucleusZ = right.theRestNucleusZ;
  CoulombBarrier = right.CoulombBarrier;
  MaximalKineticEnergy = right.MaximalKineticEnergy;
  ExcitonLevelDensityRatio = right.ExcitonLevelDensityRatio;
  EmissionProbability = right.EmissionProbability;
  CondensationProbability = right.CondensationProbability;
  Momentum = right.Momentum;
}


G4VPreCompoundFragment::G4VPreCompoundFragment(const G4double anA, const G4double aZ):
  theA(anA),theZ(aZ),theRestNucleusA(0.0),theRestNucleusZ(0.0),CoulombBarrier(0.0),
  MaximalKineticEnergy(-1.0),ExcitonLevelDensityRatio(0.0),EmissionProbability(0.0),
  CondensationProbability(0.0),Momentum(0.0,0.0,0.0,0.0)
{}



G4VPreCompoundFragment::~G4VPreCompoundFragment()
{
}


const G4VPreCompoundFragment & G4VPreCompoundFragment::operator=
(const G4VPreCompoundFragment & right)
{
  if (this != &right) {
    theA = right.theA;
    theZ = right.theZ;
    theRestNucleusA = right.theRestNucleusA;
    theRestNucleusZ = right.theRestNucleusZ;
    CoulombBarrier = right.CoulombBarrier;
    MaximalKineticEnergy = right.MaximalKineticEnergy;
    ExcitonLevelDensityRatio = right.ExcitonLevelDensityRatio;
    EmissionProbability = right.EmissionProbability;
    CondensationProbability = right.CondensationProbability;
    Momentum = right.Momentum;
  }
  return *this;
}

G4int G4VPreCompoundFragment::operator==(const G4VPreCompoundFragment & right) const
{
  return (this == (G4VPreCompoundFragment *) &right);
}

G4int G4VPreCompoundFragment::operator!=(const G4VPreCompoundFragment & right) const
{
  return (this != (G4VPreCompoundFragment *) &right);
}


G4std::ostream& operator << (G4std::ostream &out, const G4VPreCompoundFragment &theFragment)
{
  out << &theFragment;
  return out; 
}


G4std::ostream& operator << (G4std::ostream &out, const G4VPreCompoundFragment *theFragment)
{
  long old_floatfield = out.setf(0,G4std::ios::floatfield);

  out 
    << "PreCompound Model Emitted Fragment: A = " << G4std::setprecision(3) << theFragment->theA 
    << ", Z = " << G4std::setprecision(3) << theFragment->theZ;
  out.setf(G4std::ios::scientific,G4std::ios::floatfield);
//   out
//     << ", U = " << theFragment->theExcitationEnergy/MeV 
//     << " MeV" << G4endl
//     << "          P = (" 
//     << theFragment->theMomentum.x()/MeV << ","
//     << theFragment->theMomentum.y()/MeV << ","
//     << theFragment->theMomentum.z()/MeV 
//     << ") MeV   E = " 
//     << theFragment->theMomentum.t()/MeV << " MeV";

  out.setf(old_floatfield,G4std::ios::floatfield);

  return out;
    
}


void G4VPreCompoundFragment::Init(const G4Fragment & aFragment)
{
  
  SetRestA(aFragment.GetA());
  SetRestZ(aFragment.GetZ());

  if ((GetRestA() < GetRestZ()) ||
      (GetRestA() < GetA()) ||
	  (GetRestZ() < GetZ())) {
    // for to be sure that emission probability will be 0.
    SetMaximalKineticEnergy(0.0);
    return;
  }
      // Compute nuclear radius (needed to calculate Coulomb barrier)
  G4double NuclearRadius = 2.173*
    (1.0+0.006103*GetZ()*GetRestZ())/
    (1.0+0.009443*GetZ()*GetRestZ());
      // Calculate Coulomb barrier
  SetCoulombBarrier(CalcCoulombBarrier(NuclearRadius));
  
  // Compute Binding Energies for fragments (needed to separate a fragment from the nucleus)
  
  SetBindingEnergy(G4NucleiProperties::GetMassExcess(GetA(),GetZ())/MeV+
		   G4NucleiProperties::GetMassExcess(GetRestA(),GetRestZ())/MeV-
		   G4NucleiProperties::GetMassExcess(aFragment.GetA(),aFragment.GetZ())/MeV);
      
  // Compute Maximal Kinetic Energy which can be carried by fragments after separation
  SetMaximalKineticEnergy((aFragment.GetExcitationEnergy()/MeV)-
			  (GetBindingEnergy()+
					       GetCoulombBarrier()));
  
  
}

G4double G4VPreCompoundFragment::CalcCoulombBarrier(const G4double & NucRad)
  // Calculation of Coulomb potential energy (barrier) for outgoing particles (in MeV)
{
  // for neutron 
  G4double Barrier;
  if (GetZ() == 0) 
  {
   Barrier = 0.0;
  }
  else
  {
    Barrier = (1.44/NucRad)*
              ((GetZ()*GetRestZ())/
               (pow(GetA(),1.0/3.0)+pow(GetRestA(),1.0/3.0)));
  }
  //  return Barrier;
  return Barrier;
}


G4double G4VPreCompoundFragment::
CalcEmissionProbability(const G4Fragment & aFragment)
{
  if (GetMaximalKineticEnergy() <= 0.0) return 0.0;

  // Coulomb barrier for fragment "index" (for nucleon it is 0) is the lower limit 
  // of integration over kinetic energy
  G4double LowerLimit = GetCoulombBarrier();
  
  // Excitation energy of nucleus after fragment emission is the upper limit
  // of integration over kinetic energy
  G4double UpperLimit = aFragment.GetExcitationEnergy()/MeV - GetBindingEnergy();
  
  return EmissionProbability = IntegrateEmissionProbability(LowerLimit,UpperLimit,aFragment);

  // Correction by V. Krylov:
  //      Emission probability was stored in a local variable but not in
  //      data member EmissionProbability, then the returned probability 
  //      by GetEmissionProbability() method was always 0.0

}

G4double G4VPreCompoundFragment::
IntegrateEmissionProbability(const G4double & Low, const G4double & Up,
			     const G4Fragment & aFragment)
{
  static const G4double w[8] = {0.1012285363,
				0.2223810345,
				0.3137066459,
				0.3626837834,
				0.3626837834,
				0.3137066459,
				0.2223810345,
				0.1012285363};

  static const G4double FIKS[8] = {0.9602898565,
				   0.7966664774,
				   0.5255324099,
				   0.1834346425,
				   -0.1834346425,
				   -0.5255324099,
				   -0.7966664774,
				   -0.9602898565};

  G4double Total = 0.0;
  for (G4int i = 0; i < 8; i++) {
    G4double KineticE = ((Up-Low)*FIKS[i]+(Up+Low))/2.0;
    Total += w[i]*ProbabilityDistributionFunction(KineticE, aFragment)
                 *(Up-Low)/2.0;
  }
  return Total;
}


