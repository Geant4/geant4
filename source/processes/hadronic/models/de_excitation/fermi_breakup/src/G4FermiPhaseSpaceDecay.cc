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
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4FermiPhaseSpaceDecay.hh"
#include "G4HadronicException.hh"
#include "Randomize.hh"

#include <algorithm>
#include <numeric>
#include <functional>

G4FermiPhaseSpaceDecay::G4FermiPhaseSpaceDecay()
{
}

G4FermiPhaseSpaceDecay::~G4FermiPhaseSpaceDecay()
{
}

std::vector<G4LorentzVector*> *
G4FermiPhaseSpaceDecay::KopylovNBodyDecay(const G4double M, const std::vector<G4double>& m) const
  // Calculates momentum for N fragments (Kopylov's method of sampling is used)
{
  G4int N = m.size();

  std::vector<G4LorentzVector*>* P = new std::vector<G4LorentzVector*>;
  P->insert(P->begin(), N, static_cast<G4LorentzVector*>(0));

  G4double mtot = std::accumulate( m.begin(), m.end(), 0.0);
  G4double mu = mtot;
  G4double PFragMagCM = 0.0;
  G4double Mass = M;
  G4double T = M-mtot;
  G4LorentzVector PFragCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PFragLab(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestLab(0.0,0.0,0.0,Mass);

  for (G4int k = N-1; k > 0; k--)
    {
      mu -= m[k];
      if (k>1) T *= BetaKopylov(k);
      else T = 0.0;
      
      G4double RestMass = mu + T;

      PFragMagCM = PtwoBody(Mass,m[k],RestMass);
      if (PFragMagCM < 0) 
        {
          throw G4HadronicException(__FILE__, __LINE__, "G4FermiPhaseSpaceDecay::KopylovNBodyDecay: Error sampling fragments momenta!!");
        }
      

      // Create a unit vector with a random direction isotropically distributed
      G4ParticleMomentum RandVector(IsotropicVector(PFragMagCM));

      PFragCM.setVect(RandVector);
      PFragCM.setE(std::sqrt(RandVector.mag2()+m[k]*m[k]));

      PRestCM.setVect(-RandVector);
      PRestCM.setE(std::sqrt(RandVector.mag2()+RestMass*RestMass));


      G4ThreeVector BoostV = PRestLab.boostVector();

      PFragLab = PFragCM;
      PFragLab.boost(BoostV);
      PRestLab = PRestCM;
      PRestLab.boost(BoostV);

      P->operator[](k) = new G4LorentzVector(PFragLab);
      
      Mass = RestMass;
    }

    P->operator[](0) = new G4LorentzVector(PRestLab);

  return P;

}



std::vector<G4LorentzVector*> *
G4FermiPhaseSpaceDecay::NBodyDecay(const G4double M,  const std::vector<G4double>& m) const
{
  // Number of fragments
  G4int N = m.size();
  G4int i, j;
  // Total Daughters Mass
  G4double mtot = std::accumulate( m.begin(), m.end(), 0.0);
  G4double Emax = M - mtot + m[0];
  G4double Emin = 0.0;
  G4double Wmax = 1.0;
  for (i = 1; i < N; i++)
    {
      Emax += m[i];
      Emin += m[i-1];
      Wmax *= this->PtwoBody(Emax, Emin, m[i]);
    }

  G4int ntries = 0;
  G4double weight = 1.0;
  std::vector<G4double> p(N);
  do
    {
      // Sample uniform random numbers in increasing order
      std::vector<G4double> r;
      r.reserve(N);
      r.push_back(0.0);
      for (i = 1; i < N-1; i++) r.push_back(G4UniformRand());
      r.push_back(1.0);
      std::sort(r.begin(),r.end(), std::less<G4double>());
      
      // Calculate virtual masses
      std::vector<G4double> vm(N);
      vm[0] = 0.0;
      std::partial_sum(m.begin(), m.end(), vm.begin());
      std::transform(r.begin(), r.end(), r.begin(), std::bind2nd(std::multiplies<G4double>(), M-mtot));
      std::transform(r.begin(), r.end(), vm.begin(), vm.begin(), std::plus<G4double>());
      r.clear();

      // Calcualte daughter momenta
      weight = 1.0;
      for (j = 0; j < N-1; j++)
	{
	  p[j] = PtwoBody(vm[j+1],vm[j],m[j+1]); 
	  if (p[j] < 0.0)
	    {
	      G4cerr << "G4FermiPhaseSpaceDecay::Decay: Daughter momentum less than zero\n";
	      weight = 0.0;
	      break;
	    }
	  else
	    {
	      weight *= p[j]; 
	    }
	}
      p[N-1] = PtwoBody(vm[N-2], m[N-2], m[N-1]);

           
      if (ntries++ > 1000000)
	{
	  throw G4HadronicException(__FILE__, __LINE__, "G4FermiPhaseSpaceDecay::Decay: Cannot determine decay kinematics");
	}
    }
  while ( weight < G4UniformRand()*Wmax );

  std::vector<G4LorentzVector*> * P = new std::vector<G4LorentzVector*>;  
  P->insert(P->begin(),N, static_cast<G4LorentzVector*>(0));

  G4ParticleMomentum a3P = this->IsotropicVector(p[0]);

  P->operator[](0) = new G4LorentzVector( a3P, std::sqrt(a3P.mag2()+m[0]*m[0]) );  
  P->operator[](1) = new G4LorentzVector(-a3P, std::sqrt(a3P.mag2()+m[1]*m[1]) );
  for (i = 2; i < N; i++)
    {
      a3P = this->IsotropicVector(p[i-1]);
      P->operator[](i) = new G4LorentzVector(a3P, std::sqrt(a3P.mag2() + m[i]*m[i]));
      G4ThreeVector Beta = (-1.0)*P->operator[](i)->boostVector();
      // boost already created particles
      for (j = 0; j < i; j++)
	{
	  P->operator[](j)->boost(Beta);
	}
    }
  
  return P;
}




std::vector<G4LorentzVector*> *
G4FermiPhaseSpaceDecay::
TwoBodyDecay(const G4double M, const std::vector<G4double>& m) const
{
  G4double m0 = m.front();
  G4double m1 = m.back();
  G4double psqr = this->PtwoBody(M,m0,m1);
  G4ParticleMomentum p = this->IsotropicVector(std::sqrt(psqr));

  G4LorentzVector * P41 = new G4LorentzVector;
  P41->setVect(p);
  P41->setE(std::sqrt(psqr+m0*m0));

  G4LorentzVector * P42 = new G4LorentzVector;
  P42->setVect(-p);
  P42->setE(std::sqrt(psqr+m1*m1));

  std::vector<G4LorentzVector*> * result = new std::vector<G4LorentzVector*>;
  result->push_back(P41);
  result->push_back(P42);
  return result;
}




G4ParticleMomentum G4FermiPhaseSpaceDecay::IsotropicVector(const G4double Magnitude) const
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


