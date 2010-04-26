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
// $Id: G4FermiConfiguration.cc,v 1.13 2010-04-26 11:14:28 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// J. M. Quesada (12 October 2009) new implementation of Gamma function in configuration weight 
// J. M. Quesada (09 March 2010) Kappa is set to 6.

#include "G4FermiConfiguration.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include <set>

// Kappa = V/V_0 it is used in calculation of Coulomb energy
// Kappa is adimensional
// JMQ 090310 according to model developer (A. Botvina) no theoretical constraint for kappa below 10
// kappa values  larger than 1 seem to provide better results. 6 is a good choice  
// const G4double G4FermiConfiguration::Kappa = 1.0;
const G4double G4FermiConfiguration::Kappa = 6.0;

// r0 is the nuclear radius
const G4double G4FermiConfiguration::r0 = 1.3*fermi;


G4double G4FermiConfiguration::CoulombBarrier(void)
{
  //  Calculates Coulomb Barrier (MeV) for given channel with K fragments.
  const G4double Coef = (3./5.)*(elm_coupling/r0)*std::pow(1./(1.+Kappa), 1./3.);

  G4double SumA = 0;
  G4double SumZ = 0;
  G4double CoulombEnergy = 0.;
  for (std::vector<const G4VFermiFragment*>::iterator i = Configuration.begin(); 
       i != Configuration.end(); i++) 
    {
      G4double z = static_cast<G4double>((*i)->GetZ());
      G4double a = static_cast<G4double>((*i)->GetA());
      CoulombEnergy += (z*z) / std::pow(a, 1./3.);
      SumA += a;
      SumZ += z;
    }
  CoulombEnergy -= SumZ*SumZ/std::pow(SumA, 1./3.);
  return -Coef * CoulombEnergy;
}




G4double G4FermiConfiguration::DecayProbability(const G4int A, const G4double TotalE)
  // Decay probability  for a given channel with K fragments
{
  // A: Atomic Weight
  // TotalE: Total energy of nucleus


  G4double KineticEnergy = TotalE; // MeV
  G4double ProdMass = 1.0;
  G4double SumMass = 0.0;
  G4double S_n = 1.0;
  std::set<G4int>   combSet;
  std::multiset<G4int> combmSet;
  
  for (std::vector<const G4VFermiFragment*>::iterator i = Configuration.begin(); 
       i != Configuration.end(); i++) 
    {
      G4int a = (*i)->GetA();
      combSet.insert(a);
      combmSet.insert(a);
      G4double m = (*i)->GetFragmentMass();
      ProdMass *= m;
      SumMass += m;
      // Spin factor S_n
      S_n *= (*i)->GetPolarization();
      KineticEnergy -= m + (*i)->GetExcitationEnergy();
    }

  // Check that there is enough energy to produce K fragments
  if (KineticEnergy <= 0.0) return 0.0;
  if ((KineticEnergy -= this->CoulombBarrier()) <= 0.0) return 0.0;


  G4double MassFactor = ProdMass/SumMass;
  MassFactor *= std::sqrt(MassFactor);  
  
  // Number of fragments
  G4int K = Configuration.size();

  // This is the constant (doesn't depend on nucleus) part
  const G4double ConstCoeff = std::pow(r0/hbarc,3.0)*Kappa*std::sqrt(2.0/pi)/3.0;
  G4double Coeff = std::pow(ConstCoeff*A,K-1);

  //JMQ 111009 Bug fixed: gamma function for odd K was wrong  by a factor 2
  // new implementation explicitely according to standard properties of Gamma function
  // Calculation of 1/Gamma(3(k-1)/2)
  G4double Gamma = 1.0;
  //  G4double arg = 3.0*(K-1)/2.0 - 1.0;
  //  while (arg > 1.1) 
  //    {
  //      Gamma *= arg; 
  //      arg--;
  //    }
  //  if ((K-1)%2 == 1) Gamma *= std::sqrt(pi);

  if ((K-1)%2 != 1)

    {
      G4double arg = 3.0*(K-1)/2.0 - 1.0;
      while (arg > 1.1) 
	{
	  Gamma *= arg; 
	  arg--;
	}
    }
  else   { 
    G4double 	n= 3.0*K/2.0-2.0;
    G4double arg2=2*n-1;
    while (arg2>1.1)
      {
	Gamma*=arg2;
	arg2-=2;
      }
    Gamma=Gamma/std::pow(2.,n)*std::sqrt(pi);
  }
  // end of new implementation of Gamma function
  
  
  // Permutation Factor G_n
  G4double G_n = 1.0;
  for (std::set<G4int>::iterator s = combSet.begin(); s != combSet.end(); ++s)
    {
      for (G4int ni = combmSet.count(*s); ni > 1; ni--) G_n *= ni;
    }

  G4double Weight = Coeff * MassFactor * (S_n / G_n) / Gamma;
  Weight *= std::pow(KineticEnergy,3.0*(K-1)/2.0)/KineticEnergy;

  return Weight; 
}


G4FragmentVector * G4FermiConfiguration::GetFragments(const G4Fragment & theNucleus)
{
 
  G4FermiPhaseSpaceDecay thePhaseSpace;

  // Calculate Momenta of K fragments
  G4double M = theNucleus.GetMomentum().m();
  std::vector<G4double> m;
  m.reserve(Configuration.size());
  std::vector<const G4VFermiFragment*>::iterator i;
  for (i = Configuration.begin(); i != Configuration.end(); ++i)
    {
      m.push_back( (*i)->GetTotalEnergy() );
    }
  std::vector<G4LorentzVector*>* MomentumComponents = 
    thePhaseSpace.Decay(M,m);
  
  G4FragmentVector * theResult = new G4FragmentVector;

  G4ThreeVector boostVector = theNucleus.GetMomentum().boostVector();  

  // Go back to the Lab Frame
  for (i = Configuration.begin(); i != Configuration.end(); ++i)
    {
#ifdef G4NO_ISO_VECDIST
      std::vector<const G4VFermiFragment*>::difference_type n = 0;
      std::distance(Configuration.begin(), i, n);
      G4LorentzVector FourMomentum(*(MomentumComponents->operator[](n)));
#else
      G4LorentzVector FourMomentum(*(MomentumComponents->
				     operator[](std::distance(Configuration.begin(),i))));
#endif
    
      // Lorentz boost
      FourMomentum.boost(boostVector);
      
      G4FragmentVector * fragment = (*i)->GetFragment(FourMomentum);
 
  
      for (G4FragmentVector::reverse_iterator ri = fragment->rbegin();
           ri != fragment->rend(); ++ri)
        {
          theResult->push_back(*ri);
        }
      delete fragment;
    }
  
  if (!MomentumComponents->empty())
    {
      std::for_each(MomentumComponents->begin(),MomentumComponents->end(),
		      DeleteFragment());
    }
  
  delete MomentumComponents;
  
  return theResult;
}

    
G4ParticleMomentum G4FermiConfiguration::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = std::sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ParticleMomentum Vector(Magnitude*std::cos(Phi)*SinTheta,
			    Magnitude*std::sin(Phi)*SinTheta,
			    Magnitude*CosTheta);
  return Vector;
}
