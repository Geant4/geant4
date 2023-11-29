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
#include "Randomize.hh"
#include "G4PhysicsModelCatalog.hh"

const G4int G4UnstableFragmentBreakUp::Zfr[] = {0, 1, 1, 1, 2, 2};
const G4int G4UnstableFragmentBreakUp::Afr[] = {1, 1, 2, 3, 3, 4};

G4UnstableFragmentBreakUp::G4UnstableFragmentBreakUp() : fVerbose(1), fSecID(-1)
{ 
  fLevelData = G4NuclearLevelData::GetInstance();
  for(G4int i=0; i<6; ++i) {
    masses[i] = G4NucleiProperties::GetNuclearMass(Afr[i], Zfr[i]);
  }
  fSecID = G4PhysicsModelCatalog::GetModelID("model_G4UnstableFragmentBreakUp");
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

  G4LorentzVector lv = nucleus->GetMomentum();
  G4double time = nucleus->GetCreationTime();

  G4double mass1(0.0), mass2(0.0);
 
  // look for the decay channel with normal masses
  // without Coulomb barrier and paring corrections
  // 1 - recoil, 2 - emitted light ion 
  if(fVerbose > 1) {
    G4cout << "#Unstable decay " << " Z= " << Z << " A= " << A 
	   << " Eex(MeV)= " << nucleus->GetExcitationEnergy() << G4endl;
  }
  const G4double tolerance = 10*CLHEP::eV;
  const G4double dmlimit   = 0.005*CLHEP::MeV;
  G4double mass = lv.mag();
  G4double exca = -1000.0;
  G4bool isChannel = false;
  G4int idx = -1;
  for(G4int i=0; i<6; ++i) {
    G4int Zres = Z - Zfr[i];
    G4int Ares = A - Afr[i];
    if(Zres >= 0 && Ares >= Zres && Ares >= Afr[i]) {
      if(Ares <= 4) {
	for(G4int j=0; j<6; ++j) {
	  if(Zres == Zfr[j] && Ares == Afr[j]) {
	    G4double delm = mass - masses[i] - masses[j];
	    /*
	    G4cout << "i=" << i << " j=" << j << " Zres=" << Zres
		   << " Ares=" << Ares << " delm=" << delm << G4endl;
	    */
	    if(delm > exca) {
	      mass2 = masses[i]; // emitted
	      mass1 = masses[j]; // recoil
              exca  = delm; 
	      idx = i;
              if(delm > 0.0) { 
		isChannel = true;
		break;
	      }
	    }
	  }
	}
      }
      if(isChannel) { break; }
      // no simple channel
      G4double mres = G4NucleiProperties::GetNuclearMass(Ares, Zres);
      G4double e = mass - mres - masses[i];
      // select excited state
      //const G4LevelManager* lman = fLevelData->GetLevelManager(Zres, Ares);
      /*
      G4cout << "i=" << i << " Zres=" << Zres
	     << " Ares=" << Ares << " delm=" << e << G4endl;
      */
      if(e >= exca) {
	mass2 = masses[i];
	mass1 = (Ares > 4 && e > 0.0) ? mres + e*G4UniformRand() : mres;
        exca  = e;
	idx = i;
	if(e > 0.0) {
	  isChannel = true;
	  break;
	}
      } 
    }
  }

  G4double massmin = mass1 + mass2;
  if(fVerbose > 1) {
    G4cout << "isChannel:" << isChannel << " idx=" << idx << " Zfr=" << Zfr[idx]
	   << " Arf=" << Afr[idx] << " delm=" << mass - massmin << G4endl;
  }
  if(!isChannel || mass < massmin) { 
    if(mass + dmlimit < massmin) { return false; }
    if(fVerbose > 1) {
      G4cout << "#Unstable decay correction: Z= " << Z << " A= " << A 
             << " idx= " << idx
             << " deltaM(MeV)= " << mass - massmin
             << G4endl;
    }      
    mass = massmin;
    G4double e = std::max(lv.e(), mass + tolerance);
    G4double mom = std::sqrt((e - mass)*(e + mass));
    G4ThreeVector dir = lv.vect().unit();
    lv.set(dir*mom, e);
  }

  // compute energy of light fragment
  G4double e2 = 0.5*((mass - mass1)*(mass + mass1) + mass2*mass2)/mass;
  e2 = std::max(e2, mass2);
  G4double mom = std::sqrt((e2 - mass2)*(e2 + mass2));

  // sample decay
  G4ThreeVector bst  = lv.boostVector();
  G4ThreeVector v = G4RandomDirection();
  G4LorentzVector mom2 = G4LorentzVector(v*mom, e2);
  mom2.boost(bst);  
  frag = new G4Fragment(Afr[idx], Zfr[idx], mom2);
  frag->SetCreationTime(time);
  frag->SetCreatorModelID(fSecID);
  results->push_back(frag);

  // residual
  lv -= mom2;
  Z  -= Zfr[idx];
  A  -= Afr[idx];
    
  nucleus->SetZandA_asInt(Z, A);
  nucleus->SetMomentum(lv);
  nucleus->SetCreatorModelID(fSecID);
  return true;
}

G4double G4UnstableFragmentBreakUp::GetEmissionProbability(G4Fragment*)
{
  return 0.0;
}
