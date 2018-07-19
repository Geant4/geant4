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

#include "G4FermiPhaseSpaceDecay.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4HadronicException.hh"

G4FermiPhaseSpaceDecay::G4FermiPhaseSpaceDecay()
{
  g4calc = G4Pow::GetInstance();
}

G4FermiPhaseSpaceDecay::~G4FermiPhaseSpaceDecay()
{}

std::vector<G4LorentzVector*>*
G4FermiPhaseSpaceDecay::KopylovNBodyDecay(G4double M, 
					  const std::vector<G4double>& mr) const
  // Calculates momentum for N fragments (Kopylov's method of sampling is used)
{
  size_t N = mr.size();

  std::vector<G4LorentzVector*>* P = 
    new std::vector<G4LorentzVector*>(N, 0);

  G4double mtot = 0.0;
  for(size_t k=0; k<N; ++k) { mtot += mr[k]; }
  G4double mu = mtot;
  G4double PFragMagCM = 0.0;
  G4double Mass = M;
  G4double T = Mass-mtot;
  G4LorentzVector PFragCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestLab(0.0,0.0,0.0,Mass);

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  for (size_t k = N-1; k>0; --k)
    {
      mu -= mr[k];
      if (k>1) { T *= BetaKopylov(k, rndmEngine); }
      else     { T = 0.0; }
      
      G4double RestMass = mu + T;

      PFragMagCM = PtwoBody(Mass,mr[k],RestMass);
      
      // Create a unit vector with a random direction isotropically distributed
      G4ThreeVector RandVector = PFragMagCM*G4RandomDirection();

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


