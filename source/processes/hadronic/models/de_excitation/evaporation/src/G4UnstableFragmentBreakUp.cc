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
// $Id: G4UnstableFragmentBreakUp.cc 99812 2016-10-06 08:49:27Z gcosmo $
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
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"

const G4int G4UnstableFragmentBreakUp::Zfr[] = {0, 1, 1, 1, 2, 2};
const G4int G4UnstableFragmentBreakUp::Afr[] = {1, 1, 2, 3, 3, 4};

G4UnstableFragmentBreakUp::G4UnstableFragmentBreakUp()
{ 
  fLevelData = G4NuclearLevelData::GetInstance();
  for(G4int i=0; i<6; ++i) {
    masses[i] = G4NucleiProperties::GetNuclearMass(Afr[i], Zfr[i]);
  }
}

G4UnstableFragmentBreakUp::~G4UnstableFragmentBreakUp()
{}

G4bool G4UnstableFragmentBreakUp::BreakUpChain(G4FragmentVector* results, 
					       G4Fragment* nucleus)
{
  //G4cout << "G4UnstableFragmentBreakUp::EmittedFragment" << G4endl;
  G4Fragment* frag = nullptr;
  
  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();

  // if the isotope is in the database it is not exotic
  // so, cannot be handled by this class
  if(fLevelData->GetLevelManager(Z, A)) { return false; }

  G4LorentzVector lv = nucleus->GetMomentum();
  G4double time = nucleus->GetCreationTime();

  G4double mass1(0.0), mass2(0.0);
 
  G4double mass = lv.mag();
  G4bool   done = false;

  G4int Amax = A;
  G4int Zres = Z;
  G4int Ares = A;

  for(G4int k=0; k<Amax; ++k) {
    // look for the decay channel with normal masses
    // without Coulomb barrier and paring corrections

    // 1 - recoil, 2 - emitted light ion 

    //G4cout << "Unstable decay #" << k << " Z= " << Z << " A= " << A << G4endl;

    G4double ekin = 0.0;
    G4int index = -1;
    for(G4int i=0; i<6; ++i) {
      Zres = Z - Zfr[i];
      Ares = A - Afr[i];
      for(G4int j=0; j<6; ++j) {
	if(Zres == Zfr[j] && Ares == Afr[j]) {
	  /*
	  G4cout << "i= " << i << " j= " << j << " Zres= " << Zres
		 << " Ares= " << Ares << " dm= " << mass - masses[i] - masses[j]
		 << G4endl;
	  */ 
	  if(mass >= masses[i] + masses[j]) {
 	    mass2 = masses[i];
	    mass1 = masses[j];
	    ekin  = 0.0;
	    done = true;
	    break;
	  }
	}
      }
      if(done) {
	index = i;
        break;
      }
      if(Zres >= 0 && Ares >= Zres && Ares > 0) {
	G4double m = G4NucleiProperties::GetNuclearMass(Ares, Zres);
	G4double e = mass - m - masses[i];
	if(e > ekin) {
	  mass2 = masses[i];
	  mass1 = m;
	  ekin  = e;
	  index = i;
	  if(fLevelData->GetLevelManager(Z, A)) {
	    done = true;
            break;
	  }	
	}
      }
    }

    // no decay channel - assume that primary mass is biased
    // only energy will be conserved
    if(index < 0) {
      for(G4int i=0; i<6; ++i) {
	Zres = Z - Zfr[i];
	Ares = A - Afr[i];
	if(Zres >= 0 && Ares >= Zres && Ares > 0) {
          
	  G4double m = proton_mass_c2*Zres + neutron_mass_c2*(Ares - Zres);
	  G4double e = mass - m - masses[i];
	  if(e > ekin) {
	    mass2 = masses[i];
	    mass1 = m;
	    ekin  = e;
	    index = i;
	  }
	}
      }
    }
    if(index < 0) {

      // further decay impossible
      // if no one decay sampled do not update primary 
      if(0 == k) { return false; }
      else       { break; }
    }

    // useful to left max excitation for the residual
    mass1 += ekin;

    // compute energy of light fragment
    G4double e2 = 0.5*((mass - mass1)*(mass + mass1) + mass2*mass2)/mass;
    e2 = std::max(e2, mass2);
    G4double mom = std::sqrt((e2 - mass2)*(e2 + mass2));

    // sample decay
    G4ThreeVector bst  = lv.boostVector();
    G4ThreeVector v = G4RandomDirection();
    G4LorentzVector mom2 = G4LorentzVector(v*mom, e2);
    mom2.boost(bst);  
    frag = new G4Fragment(Afr[index], Zfr[index], mom2);
    frag->SetCreationTime(time);
    results->push_back(frag);

    // residual
    lv -= mom2;
    Z  -= Zfr[index];
    A  -= Afr[index];
    
    mass = lv.mag();

    if(done) { break; }
  }

  // updated primary
  nucleus->SetZandA_asInt(Z, A);
  nucleus->SetMomentum(lv);
  return true;
}

G4double G4UnstableFragmentBreakUp::GetEmissionProbability(G4Fragment*)
{
  return 0.0;
}
