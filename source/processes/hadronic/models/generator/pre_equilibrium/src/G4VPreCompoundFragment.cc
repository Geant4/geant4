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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPreCompoundFragment.cc,v 1.8.2.1 2001/06/28 19:13:35 gunter Exp $
// GEANT4 tag $Name:  $
//
// by V. Lara
 
#include "G4VPreCompoundFragment.hh"


G4VPreCompoundFragment::G4VPreCompoundFragment(const G4VPreCompoundFragment & right)
{
  theA = right.theA;
  theZ = right.theZ;
  theRestNucleusA = right.theRestNucleusA;
  theRestNucleusZ = right.theRestNucleusZ;
  theCoulombBarrier = right.theCoulombBarrier;
  theCoulombBarrierPtr = right.theCoulombBarrierPtr;
  theMaximalKineticEnergy = right.theMaximalKineticEnergy;
  theExcitonLevelDensityRatio = right.theExcitonLevelDensityRatio;
  theEmissionProbability = right.theEmissionProbability;
  theCondensationProbability = right.theCondensationProbability;
  theMomentum = right.theMomentum;
  theFragmentName = right.theFragmentName;
}


G4VPreCompoundFragment::G4VPreCompoundFragment(const G4double anA,
					       const G4double aZ, G4VCoulombBarrier* aCoulombBarrier):
  theA(anA),theZ(aZ), theCoulombBarrierPtr(aCoulombBarrier),
  theRestNucleusA(0.0),theRestNucleusZ(0.0),theCoulombBarrier(0.0),
  theMaximalKineticEnergy(-1.0),theExcitonLevelDensityRatio(0.0),theEmissionProbability(0.0),
  theCondensationProbability(0.0),theMomentum(0.0,0.0,0.0,0.0),theFragmentName("No Name")
{}

G4VPreCompoundFragment::G4VPreCompoundFragment(const G4double anA,
					       const G4double aZ, G4VCoulombBarrier* aCoulombBarrier,
					       const G4String & aName):
  theA(anA),theZ(aZ), theCoulombBarrierPtr(aCoulombBarrier),
  theRestNucleusA(0.0),theRestNucleusZ(0.0),theCoulombBarrier(0.0),
  theMaximalKineticEnergy(-1.0),theExcitonLevelDensityRatio(0.0),theEmissionProbability(0.0),
  theCondensationProbability(0.0),theMomentum(0.0,0.0,0.0,0.0),theFragmentName(aName)
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
    theCoulombBarrier = right.theCoulombBarrier;
    theCoulombBarrierPtr = right.theCoulombBarrierPtr;
    theMaximalKineticEnergy = right.theMaximalKineticEnergy;
    theExcitonLevelDensityRatio = right.theExcitonLevelDensityRatio;
    theEmissionProbability = right.theEmissionProbability;
    theCondensationProbability = right.theCondensationProbability;
    theMomentum = right.theMomentum;
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
  //     << " MeV" << endl
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
  
  theRestNucleusA = aFragment.GetA() - theA;
  theRestNucleusZ = aFragment.GetZ() - theZ;

  if ((theRestNucleusA < theRestNucleusZ) ||
      (theRestNucleusA < theA) ||
      (theRestNucleusZ < theZ)) {
    // In order to be sure that emission probability will be 0.
    theMaximalKineticEnergy = 0.0;
    return;
  }


  // Calculate Coulomb barrier
  theCoulombBarrier = theCoulombBarrierPtr->
    GetCoulombBarrier(theRestNucleusA,theRestNucleusZ,
		      aFragment.GetExcitationEnergy());
  
  // Compute Binding Energies for fragments 
  // (needed to separate a fragment from the nucleus)
  
  theBindingEnergy = G4NucleiProperties::GetMassExcess(theA,theZ) +
    G4NucleiProperties::GetMassExcess(theRestNucleusA,theRestNucleusZ) -
    G4NucleiProperties::GetMassExcess(aFragment.GetA(),aFragment.GetZ());
      
  // Compute Maximal Kinetic Energy which can be carried by fragments after separation
  //  theMaximalKineticEnergy = aFragment.GetExcitationEnergy() -
  //    (theBindingEnergy + theCoulombBarrier);
  
  G4double m = aFragment.GetMomentum().m();
  G4double rm = GetRestNuclearMass();
  G4double em = GetNuclearMass();
  theMaximalKineticEnergy = ((m - rm)*(m + rm) + em*em)/(2.0*m) - em - theCoulombBarrier;

  return;
}


G4double G4VPreCompoundFragment::CalcEmissionProbability(const G4Fragment & aFragment)
{
  if (GetMaximalKineticEnergy() <= 0.0) return 0.0;

  // Coulomb barrier is the lower limit 
  // of integration over kinetic energy
  G4double LowerLimit = theCoulombBarrier;
  
  // Excitation energy of nucleus after fragment emission is the upper limit
  // of integration over kinetic energy
  G4double UpperLimit = aFragment.GetExcitationEnergy() - theBindingEnergy;
  
  theEmissionProbability = IntegrateEmissionProbability(LowerLimit,UpperLimit,aFragment);
  
  return theEmissionProbability;
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
				0.1012285363
  };

  static const G4double FIKS[8] = { 0.9602898565,
				    0.7966664774,
				    0.5255324099,
				    0.1834346425,
				   -0.1834346425,
				   -0.5255324099,
				   -0.7966664774,
				   -0.9602898565
  };

  G4double Total = 0.0;
  for (G4int i = 0; i < 8; i++) {
    G4double KineticE = ((Up-Low)*FIKS[i]+(Up+Low))/2.0;
    Total += w[i]*ProbabilityDistributionFunction(KineticE, aFragment)*(Up-Low)/2.0;
  }
  return Total;
}





