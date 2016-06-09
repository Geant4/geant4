//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// by V. Lara
 
#include "G4PreCompoundFragment.hh"

G4PreCompoundFragment::
G4PreCompoundFragment(const G4PreCompoundFragment &right) :
  G4VPreCompoundFragment(right)
{}


G4PreCompoundFragment::
G4PreCompoundFragment(const G4double anA,
		      const G4double aZ, 
		      G4VCoulombBarrier* aCoulombBarrier,
		      const G4String & aName):
  G4VPreCompoundFragment(anA,aZ,aCoulombBarrier,aName)
{}



G4PreCompoundFragment::~G4PreCompoundFragment()
{
}


const G4PreCompoundFragment & G4PreCompoundFragment::
operator= (const G4PreCompoundFragment & right)
{
  if (&right != this) this->G4VPreCompoundFragment::operator=(right);
  return *this;
}

G4int G4PreCompoundFragment::operator==(const G4PreCompoundFragment & right) const
{
  return G4VPreCompoundFragment::operator==(right);
}

G4int G4PreCompoundFragment::operator!=(const G4PreCompoundFragment & right) const
{
  return G4VPreCompoundFragment::operator!=(right);
}


G4double G4PreCompoundFragment::
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
  
  theEmissionProbability = 
    IntegrateEmissionProbability(LowerLimit,UpperLimit,aFragment);
  
  return theEmissionProbability;
}

G4double G4PreCompoundFragment::
IntegrateEmissionProbability(const G4double & Low, const G4double & Up,
			     const G4Fragment & aFragment)
{
  static const G4int N = 10;
  // 10-Points Gauss-Legendre abcisas and weights
  static const G4double w[N] = {
    0.0666713443086881,
    0.149451349150581,
    0.219086362515982,
    0.269266719309996,
    0.295524224714753,
    0.295524224714753,
    0.269266719309996,
    0.219086362515982,
    0.149451349150581,
    0.0666713443086881
  };
  static const G4double x[N] = {
    -0.973906528517172,
    -0.865063366688985,
    -0.679409568299024,
    -0.433395394129247,
    -0.148874338981631,
    0.148874338981631,
    0.433395394129247,
    0.679409568299024,
    0.865063366688985,
    0.973906528517172
  };
  
  G4double Total = 0.0;


  for (G4int i = 0; i < N; i++) 
    {
      G4double KineticE = ((Up-Low)*x[i]+(Up+Low))/2.0;
      Total += w[i]*ProbabilityDistributionFunction(KineticE, aFragment);
    }
  Total  *= (Up-Low)/2.0;

  return Total;
}




G4double G4PreCompoundFragment::
GetKineticEnergy(const G4Fragment & aFragment) 
{
  G4double V = this->GetCoulombBarrier();
  G4double Tmax =  this->GetMaximalKineticEnergy() - V;
  
  G4double T(0.0);
  G4double NormalizedProbability(1.0);
  do 
    {
      T = V + G4UniformRand()*Tmax;
      NormalizedProbability = this->ProbabilityDistributionFunction(T,aFragment)/
	this->GetEmissionProbability();
      
    } 
  while (G4UniformRand() > NormalizedProbability);
  
  return T;
}
