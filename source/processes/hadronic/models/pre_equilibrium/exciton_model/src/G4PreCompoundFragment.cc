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
// $Id: G4PreCompoundFragment.cc 68028 2013-03-13 13:48:15Z gcosmo $
//
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modified:
// 06.09.2008 JMQ Also external choice has been added for:
//               - superimposed Coulomb barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko cleanup
//

#include "G4PreCompoundFragment.hh"

G4PreCompoundFragment::
G4PreCompoundFragment(const G4ParticleDefinition* part,
		      G4VCoulombBarrier* aCoulombBarrier)
  : G4VPreCompoundFragment(part,aCoulombBarrier)
{}

G4PreCompoundFragment::~G4PreCompoundFragment()
{}

G4double G4PreCompoundFragment::
CalcEmissionProbability(const G4Fragment & aFragment)
{
  //G4cout << theCoulombBarrier << "  " << GetMaximalKineticEnergy() << G4endl;
  // If  theCoulombBarrier effect is included in the emission probabilities
  //if (GetMaximalKineticEnergy() <= 0.0) 
  G4double limit = 0.0;
  if(OPTxs==0 ||  useSICB) { limit = theCoulombBarrier; }
  if (GetMaximalKineticEnergy() <= limit) 
    {
      theEmissionProbability = 0.0;
      return 0.0;
    }    
  // If  theCoulombBarrier effect is included in the emission probabilities
  //  G4double LowerLimit = 0.;
  // Coulomb barrier is the lower limit 
  // of integration over kinetic energy
  G4double LowerLimit = limit;

  // Excitation energy of nucleus after fragment emission is the upper 
  //limit of integration over kinetic energy
  G4double UpperLimit = GetMaximalKineticEnergy();
  
  theEmissionProbability = 
    IntegrateEmissionProbability(LowerLimit,UpperLimit,aFragment);
  /*
  G4cout << "## G4PreCompoundFragment::CalcEmisProb "
         << "Z= " << aFragment.GetZ_asInt() 
	 << " A= " << aFragment.GetA_asInt()
	 << " Elow= " << LowerLimit/MeV
	 << " Eup= " << UpperLimit/MeV
	 << " prob= " << theEmissionProbability
	 << G4endl;
  */
  return theEmissionProbability;
}

G4double G4PreCompoundFragment::
IntegrateEmissionProbability(G4double Low, G4double Up,
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

  for (G4int i = 0; i < N; ++i) 
    {
      G4double KineticE = 0.5*((Up-Low)*x[i]+(Up+Low));
      Total += w[i]*ProbabilityDistributionFunction(KineticE, aFragment);
    }
  Total  *= 0.5*(Up-Low);
  return Total;
}

G4double G4PreCompoundFragment::
GetKineticEnergy(const G4Fragment & aFragment) 
{
  //let's keep this way for consistency with CalcEmissionProbability method
  G4double V = 0.0;
  if(OPTxs==0 || useSICB) { V = theCoulombBarrier; }

  G4double Tmax = GetMaximalKineticEnergy();  
  if(Tmax < V) { return 0.0; }
  G4double T(0.0);
  G4double Probability(1.0);
  G4double maxProbability = GetEmissionProbability();
  do 
    {
      T = V + G4UniformRand()*(Tmax-V);
      Probability = ProbabilityDistributionFunction(T,aFragment);      
    }   while (maxProbability*G4UniformRand() > Probability);  
  return T;
}
