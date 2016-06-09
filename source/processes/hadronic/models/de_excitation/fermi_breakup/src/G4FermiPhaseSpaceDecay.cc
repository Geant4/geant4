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
// $Id: G4ExcitationHandler.hh,v 1.13 2010-11-17 16:20:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Phase space decay for the Fermi BreakUp model
// by V. Lara
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko: 
//          - IsotropicVector is inlined
//          - Momentum computation return zero or positive value
//          - DumpProblem method is added providing more information
//          - Reduced usage of exotic std functions  

#include <numeric>

#include "G4FermiPhaseSpaceDecay.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"

G4FermiPhaseSpaceDecay::G4FermiPhaseSpaceDecay()
{
  g4pow = G4Pow::GetInstance();
}

G4FermiPhaseSpaceDecay::~G4FermiPhaseSpaceDecay()
{
}

std::vector<G4LorentzVector*> *
G4FermiPhaseSpaceDecay::KopylovNBodyDecay(const G4double M, 
					  const std::vector<G4double>& mr) const
  // Calculates momentum for N fragments (Kopylov's method of sampling is used)
{
  size_t N = mr.size();

  std::vector<G4LorentzVector*>* P = new std::vector<G4LorentzVector*>(N, 0);

  G4double mtot = std::accumulate( mr.begin(), mr.end(), 0.0);
  G4double mu = mtot;
  G4double PFragMagCM = 0.0;
  G4double Mass = M;
  G4double T = Mass-mtot;
  G4LorentzVector PFragCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestLab(0.0,0.0,0.0,Mass);

  for (size_t k = N-1; k>0; --k)
    {
      mu -= mr[k];
      if (k>1) { T *= BetaKopylov(k); }
      else     { T = 0.0; }
      
      G4double RestMass = mu + T;

      PFragMagCM = PtwoBody(Mass,mr[k],RestMass);
      
      // Create a unit vector with a random direction isotropically distributed
      G4ThreeVector RandVector(IsotropicVector(PFragMagCM));

      PFragCM.setVect(RandVector);
      PFragCM.setE(std::sqrt(PFragMagCM*PFragMagCM + mr[k]*mr[k]));

      PRestCM.setVect(-RandVector);
      PRestCM.setE(std::sqrt(PFragMagCM*PFragMagCM + RestMass*RestMass));


      G4ThreeVector BoostV = PRestLab.boostVector();

      PFragCM.boost(BoostV);
      PRestCM.boost(BoostV);
      PRestLab = PRestCM;

      (*P)[k] = new G4LorentzVector(PFragCM);
      
      Mass = RestMass;
    }

  (*P)[0] = new G4LorentzVector(PRestLab);

  return P;
}


std::vector<G4LorentzVector*> *
G4FermiPhaseSpaceDecay::NBodyDecay(G4double M,  const std::vector<G4double>& mr) const
{
  // Number of fragments
  size_t N = mr.size();
  size_t i, j;
  // Total Daughters Mass
  G4double mtot = std::accumulate( mr.begin(), mr.end(), 0.0);
  G4double Emax = M - mtot + mr[0];
  G4double Emin = 0.0;
  G4double Wmax = 1.0;
  for (i = 1; i < N; i++)
    {
      Emax += mr[i];
      Emin += mr[i-1];
      Wmax *= PtwoBody(Emax, Emin, mr[i]);
    }

  G4int ntries = 0;
  G4double weight = 1.0;
  std::vector<G4double> p(N, 0.0);
  std::vector<G4double> r(N,0.0);
  std::vector<G4double> vm(N, 0.0);
  r[N-1] = 1.0;
 
  do
    {
      // Sample uniform random numbers in increasing order
      for (i = 1; i < N-1; i++) { r[i] = G4UniformRand(); }
      std::sort(r.begin(),r.end(), std::less<G4double>());
      
      // Calculate virtual masses
      std::partial_sum(mr.begin(), mr.end(), vm.begin());
      std::transform(r.begin(), r.end(), r.begin(), 
		     std::bind2nd(std::multiplies<G4double>(), M-mtot));
      std::transform(r.begin(), r.end(), vm.begin(), vm.begin(), std::plus<G4double>());

      // Calcualte daughter momenta
      weight = 1.0;
      for (j = 0; j < N-1; j++)
	{
	  p[j] = PtwoBody(vm[j+1],vm[j],mr[j+1]); 
	  weight *= p[j]; 
	}
      p[N-1] = PtwoBody(vm[N-2],mr[N-2],mr[N-1]);

           
      if (ntries++ > 1000000)
	{
	  throw G4HadronicException(__FILE__, __LINE__, "Failed to decay");
	}
    }
  while ( weight < G4UniformRand()*Wmax );

  std::vector<G4LorentzVector*> * P = new std::vector<G4LorentzVector*>(N, 0);  

  G4ThreeVector a3P = IsotropicVector(p[0]);

  (*P)[0] = new G4LorentzVector( a3P, std::sqrt(a3P.mag2()+mr[0]*mr[0]) );
  (*P)[1] = new G4LorentzVector(-a3P, std::sqrt(a3P.mag2()+mr[1]*mr[1]) );
  for (i = 2; i < N; i++)
    {
      a3P = IsotropicVector(p[i-1]);
      (*P)[i] = new G4LorentzVector(a3P, std::sqrt(a3P.mag2() + mr[i]*mr[i]));
      G4ThreeVector Beta = -((*P)[i]->boostVector());
      // boost already created particles
      for (j = 0; j < i; j++)
	{
	  (*P)[j]->boost(Beta);
	}
    }
  
  return P;
}

std::vector<G4LorentzVector*> *
G4FermiPhaseSpaceDecay::TwoBodyDecay(G4double M, 
				     const std::vector<G4double>& mass) const
{
  G4double m0 = mass.front();
  G4double m1 = mass.back();
  G4double mom = PtwoBody(M,m0,m1);
  G4ThreeVector p = IsotropicVector(mom);

  G4LorentzVector * P41 = new G4LorentzVector;
  P41->setVect(p);
  P41->setE(std::sqrt(mom*mom + m0*m0));

  G4LorentzVector * P42 = new G4LorentzVector;
  P42->setVect(-p);
  P42->setE(std::sqrt(mom*mom + m1*m1));

  std::vector<G4LorentzVector*> * result = new std::vector<G4LorentzVector*>;
  result->push_back(P41);
  result->push_back(P42);
  return result;
}

void 
G4FermiPhaseSpaceDecay::DumpProblem(G4double E, G4double P1, G4double P2, 
				    G4double P) const
{
  G4cout << "G4FermiPhaseSpaceDecay:  problem of decay of M(GeV)= " << E/GeV 
	 << " on M1(GeV)= " << P1/GeV << " and  M2(GeV)= " << P2/GeV
	 << " P(MeV)= " << P/MeV << " < 0" << G4endl;
  // exception only if the problem is numerically significant
  if(P < -CLHEP::eV) {
    throw G4HadronicException(__FILE__, __LINE__,"Error in decay kinematics");
  }
}


