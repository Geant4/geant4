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
// by V. Lara
 
#include "G4HETCFragment.hh"
#include "G4PreCompoundParameters.hh"

G4HETCFragment::
G4HETCFragment(const G4HETCFragment & right) :
  G4VPreCompoundFragment(right)
{
}

G4HETCFragment::
G4HETCFragment(const G4double anA,
	       const G4double aZ, 
	       G4VCoulombBarrier* aCoulombBarrier,
	       const G4String & aName) :
  G4VPreCompoundFragment(anA,aZ,aCoulombBarrier,aName)
{}

G4HETCFragment::~G4HETCFragment()
{
}

const G4HETCFragment & G4HETCFragment::
operator= (const G4HETCFragment & right)
{
  if (&right != this) this->G4VPreCompoundFragment::operator=(right);
  return *this;
}

G4int G4HETCFragment::operator==(const G4HETCFragment & right) const
{
  return G4VPreCompoundFragment::operator==(right);
}

G4int G4HETCFragment::operator!=(const G4HETCFragment & right) const
{
  return G4VPreCompoundFragment::operator!=(right);
}



G4double G4HETCFragment::
CalcEmissionProbability(const G4Fragment & aFragment)
{
  if (GetEnergyThreshold() <= 0.0) 
  {
      theEmissionProbability = 0.0;
      return 0.0;
  }    
  // Coulomb barrier is the lower limit 
  // of integration over kinetic energy
  G4double LowerLimit = theCoulombBarrier;
  
  // Excitation energy of nucleus after fragment emission is the upper limit
  // of integration over kinetic energy
  G4double UpperLimit = this->GetMaximalKineticEnergy();
  
  theEmissionProbability = IntegrateEmissionProbability(LowerLimit,UpperLimit,aFragment);
    
  return theEmissionProbability;
}

G4double G4HETCFragment::
IntegrateEmissionProbability(const G4double & Low, const G4double & Up,
			     const G4Fragment & aFragment)
{
  
  if ( !IsItPossible(aFragment) ) return 0.0;

  const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();
    
  G4double U = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;
  G4double Pb = P - GetA();
  G4double Nb = Pb + H;
  if (Nb <= 0.0) return 0.0;
  
  G4double A = (P*P+H*H+P-3*H)/4.0;
  G4double Ab = (Pb*Pb+H*H+Pb-3*H)/4.0;
  U = std::max(U-A,0.0);
  if (U <= 0.0) return 0.0;

  G4double g = (6.0/pi2)*aFragment.GetA()*G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  G4double gb = (6.0/pi2)*GetRestA()*G4PreCompoundParameters::GetAddress()->GetLevelDensity();

  G4double Pf = P;
  G4double Hf = H;
  G4double Nf = N-1.0;
  for (G4int i = 1; i < GetA(); i++)
    {
      Pf *= (P-i);
      Hf *= (H-i);
      Nf *= (N-1-i);
    }

  G4double X = std::max(Up - Ab + GetBeta(),0.0);
  G4double Y = std::max(Up - Ab - Low, 0.0);

  G4double Probability = GetSpinFactor()/(pi*hbarc*hbarc*hbarc) * GetReducedMass() * GetAlpha() *
    r0 * r0 * std::pow(GetRestA(),2.0/3.0)/std::pow(U,N-1) * (std::pow(gb,Nb)/std::pow(g,N)) * Pf * Hf * Nf * K(aFragment) *
    std::pow(Y,Nb) * (X/Nb - Y/(Nb+1));

  return Probability;
}
