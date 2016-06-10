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
// $Id: G4PreCompoundFragment.cc 94281 2015-11-10 15:00:15Z gcosmo $
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
#include "G4KalbachCrossSection.hh"
#include "G4ChatterjeeCrossSection.hh"
#include "Randomize.hh"

const G4int NPOINTSGL = 10;
// 10-Points Gauss-Legendre abcisas and weights
const G4double w[NPOINTSGL] = {
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
const G4double x[NPOINTSGL] = {
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

G4PreCompoundFragment::
G4PreCompoundFragment(const G4ParticleDefinition* part,
		      G4VCoulombBarrier* aCoulombBarrier)
  : G4VPreCompoundFragment(part,aCoulombBarrier)
{
  muu = probmax = 0.0;
  index = 0;
  if(1 == theZ) { index = theA; }
  else { index = theA + 1; }
  for(G4int i=0; i<12; ++i) { probability[i] = 0.0; }
}

G4PreCompoundFragment::~G4PreCompoundFragment()
{}

G4double G4PreCompoundFragment::
CalcEmissionProbability(const G4Fragment & aFragment)
{
  //G4cout << theCoulombBarrier << "  " << GetMaximalKineticEnergy() << G4endl;
  // If  theCoulombBarrier effect is included in the emission probabilities
  // Coulomb barrier is the lower limit of integration over kinetic energy

  G4double LowerLimit = 0.0;
  if(OPTxs==0 ||  useSICB) { LowerLimit = theCoulombBarrier; }

  if (theMaxKinEnergy <= LowerLimit) 
    {
      theEmissionProbability = 0.0;
      return 0.0;
    }    

  // compute power once
  if(OPTxs <= 2) { 
    muu =  G4ChatterjeeCrossSection::ComputePowerParameter(theResA, index);
  } else {
    muu = G4KalbachCrossSection::ComputePowerParameter(theResA, index);
  }
  
  theEmissionProbability = 
    IntegrateEmissionProbability(LowerLimit,theMaxKinEnergy,aFragment);
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
IntegrateEmissionProbability(G4double low, G4double up,
			     const G4Fragment & aFragment)
{  
  G4double sum = 0.0;
  G4double del = (up - low)*0.5;
  G4double avr = (up + low)*0.5;

  G4double e, y;
  probmax = 0.0;

  for (G4int i=0; i<NPOINTSGL; ++i) {
    e = del*x[i] + avr;
    y = ProbabilityDistributionFunction(e, aFragment);
    probability[i] = y;
    probmax = std::max(probmax, y);
    sum += w[i]*y;
  }
  return sum*del;
}

G4double G4PreCompoundFragment::CrossSection(G4double ekin) const
{
  G4double res;
  if(OPTxs == 0) { 
    res = GetOpt0(ekin);
  } else if(OPTxs <= 2) { 
    res = G4ChatterjeeCrossSection::ComputeCrossSection(ekin, theResA13, muu, 
							index, theZ, 
							theResZ, theResA); 
  } else { 
    res = G4KalbachCrossSection::ComputeCrossSection(ekin, theResA13, muu,
						     index, theZ, theA, 
						     theResZ, theResA);
  }
  return res;
}  

// *********************** OPT=0 : Dostrovski's cross section  ***************
G4double G4PreCompoundFragment::GetOpt0(G4double ekin) const
{
  G4double r0 = theParameters->Getr0()*theResA13;
  // cross section is now given in mb (r0 is in mm) for the sake of consistency
  //with the rest of the options
  return 1.e+25*CLHEP::pi*r0*r0*theResA13*GetAlpha()*(1.+GetBeta()/ekin);
}

G4double G4PreCompoundFragment::SampleKineticEnergy(const G4Fragment& fragment) 
{
  //let's keep this way for consistency with CalcEmissionProbability method
  G4double limit = 0.0;
  if(useSICB) { limit = theCoulombBarrier; }

  if(theMaxKinEnergy <= limit) { return 0.0; }

  G4double delta = theMaxKinEnergy - limit;
  static const G4double toler = 1.25;
  probmax *= toler;
  G4double prob, T(0.0);
  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  G4int i;
  for(i=0; i<100; ++i) {
    T = limit + delta*rndm->flat();
    prob = ProbabilityDistributionFunction(T, fragment);
    /*
    if(prob > probmax) { 
      G4cout << "G4PreCompoundFragment WARNING: prob= " << prob 
	     << " probmax= " << probmax << G4endl;
      G4cout << "i= " << i << " Z= " << theZ << " A= " << theA 
	     << " resZ= " << theResZ << " resA= " << theResA << "\n"
	     << " T= " << T << " Tmax= " << theMaxKinEnergy 
	     << " Tmin= " << limit
	     << G4endl;
      for(G4int i=0; i<N; ++i) { G4cout << " " << probability[i]; }
      G4cout << G4endl; 
    }
    */
    // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    if(probmax*rndm->flat() <= prob) { break; }
  }
  //G4cout << "G4PreCompoundFragment: i= " << i << "  T(MeV)= " << T 
  //<< " Emin(MeV)= " << limit << " Emax= " << theMaxKinEnergy << G4endl;
  return T;
}

