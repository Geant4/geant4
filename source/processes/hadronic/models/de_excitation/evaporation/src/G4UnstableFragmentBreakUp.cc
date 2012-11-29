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
// $Id$
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4UnstableFragmentBreakUp
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 11 May 2010
//
//Modifications: 
//      
// -------------------------------------------------------------------
//

#include "G4UnstableFragmentBreakUp.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"

G4int G4UnstableFragmentBreakUp::Zfr[6] = {0};
G4int G4UnstableFragmentBreakUp::Afr[6] = {0};
G4double G4UnstableFragmentBreakUp::masses[6] = {0.0};

G4UnstableFragmentBreakUp::G4UnstableFragmentBreakUp()
  :verbose(0)
{ 
  fNistManager = G4NistManager::Instance();
  if(0 == Afr[0]) {
    G4int z[6] = {0, 1, 1, 1, 2, 2};
    G4int a[6] = {1, 1, 2, 3, 3, 4};
    for(G4int i=0; i<6; ++i) {
      Zfr[i] = z[i];
      Afr[i] = a[i];
      masses[i] = G4NucleiProperties::GetNuclearMass(a[i], z[i]);
    }
  }
}

G4UnstableFragmentBreakUp::~G4UnstableFragmentBreakUp()
{}

G4Fragment* G4UnstableFragmentBreakUp::EmittedFragment(G4Fragment*)
{
  return 0;
}

G4FragmentVector* G4UnstableFragmentBreakUp::BreakUpFragment(G4Fragment* nucleus)
{
  //G4cout << "G4UnstableFragmentBreakUp::BreakUpFragment" << G4endl;
  G4FragmentVector * theResult = new G4FragmentVector();

  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  G4int Amax = A;
  G4LorentzVector lv = nucleus->GetMomentum();
  G4double time = nucleus->GetCreationTime();

  G4double deltaE, mass, mass1(0.0), mass2(0.0);
  G4int i, index;

  // Starts loop over evaporated particles, loop is limited by number
  // of nucleons
  for(G4int ia=0; ia<Amax; ++ia) {
 
    mass = lv.mag();
    deltaE = 0.0;
    index  = -1;
    for(i=0; i<6; ++i) {
      G4int Zres = Z - Zfr[i];
      G4int Ares = A - Afr[i];
      if(Zres >= 0 && Ares >= Zres && Ares > 0) {
	G4double m1 = G4NucleiProperties::GetNuclearMass(Ares, Zres);
	G4double de = mass - m1 - masses[i];
        if(de > deltaE) {
	  mass1 = m1;
          mass2 = masses[i];
          deltaE= de;
          index = i;
	}
      }
    }

    // no decay channels
    if(index < 0) { break; }

    // compute energy of light fragment
    G4double e2 = 0.5*((mass - mass1)*(mass + mass1) + mass2*mass2)/mass;
    if(e2 < mass2) { break; }

    // sample decay
    G4ThreeVector bst  = lv.boostVector();

    G4double cosTheta = 1. - 2. * G4UniformRand(); 
    G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
    G4double phi = twopi * G4UniformRand();
    G4double mom = std::sqrt((e2 - mass2)*(e2 + mass2));
    G4LorentzVector mom2(mom * sinTheta * std::cos(phi),
			 mom * sinTheta * std::sin(phi),
			 mom * cosTheta,
			 e2);
    mom2.boost(bst);  
    G4Fragment* fr = new G4Fragment(Afr[index], Zfr[index], mom2);
    fr->SetCreationTime(time);
    theResult->push_back(fr);

    // residual
    lv -= mom2;
    Z  -= Zfr[index];
    A  -= Afr[index];
  }

  // updated fragment
  if( theResult->size() > 0) {
    nucleus->SetZandA_asInt(Z, A);
    nucleus->SetMomentum(lv);
  }

  return theResult;
}

G4FragmentVector* G4UnstableFragmentBreakUp::BreakUp(const G4Fragment&)
{
  return 0;
}

G4double G4UnstableFragmentBreakUp::GetEmissionProbability(G4Fragment*)
{
  return 0.0;
}
