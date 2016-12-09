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
// $Id: G4GammaTransition.cc 85659 2014-11-03 10:59:10Z vnivanch $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4GammaTransition
//
//      Author V.Ivanchenko 27 February 2015
//
// -------------------------------------------------------------------

#include "G4GammaTransition.hh"
#include "G4AtomicShells.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4GammaTransition::G4GammaTransition() 
  : polarFlag(false), fDirection(0.,0.,0.)
{}

G4GammaTransition::~G4GammaTransition() 
{}
  
G4Fragment* 
G4GammaTransition::SampleTransition(G4Fragment* nucleus,
				    G4double newExcEnergy,
				    G4double mpRatio,
				    G4int  JP1,
				    G4int  JP2,
				    G4int  MP,
				    size_t shell,
				    G4bool isDiscrete,
				    G4bool isGamma,
				    G4bool isLongLived)
{
  G4Fragment* result = nullptr;
  G4double bond_energy = 0.;

  if (!isGamma) { 
    G4int Z = nucleus->GetZ_asInt();
    if(Z <= 100) {
      G4int idx = (G4int)shell;
      idx = std::min(idx, G4AtomicShells::GetNumberOfShells(Z)-1);
      bond_energy = G4AtomicShells::GetBindingEnergy(Z, idx);
    }
  }
  G4double etrans = nucleus->GetExcitationEnergy() - newExcEnergy 
    - bond_energy;
  //  G4cout << "G4GammaTransition::GenerateGamma - Etrans(MeV)= " 
  //	 << etrans << "  Eexnew= " << newExcEnergy 
  //	 << " Ebond= " << bond_energy << G4endl; 
  if(etrans <= 0.0) { 
    etrans += bond_energy;
    bond_energy = 0.0;
  }
  
  // Do complete Lorentz computation 
  G4LorentzVector lv = nucleus->GetMomentum();
  G4double mass = nucleus->GetGroundStateMass() + newExcEnergy;

  // select secondary
  G4ParticleDefinition* part;

  if(isGamma) { part =  G4Gamma::Gamma(); }
  else {
    part = G4Electron::Electron();
    G4int ne = std::max(nucleus->GetNumberOfElectrons() - 1, 0);
    nucleus->SetNumberOfElectrons(ne);
    lv += G4LorentzVector(0.0,0.0,0.0,
                          CLHEP::electron_mass_c2 - bond_energy);
  }
  /*
  G4cout << "New GammaTransition: polarFlag: " << polarFlag 
	 << " isDiscrete: " << isDiscrete << " isGamma: " << isGamma
	 << " isLongLived: " << isLongLived << G4endl;
  */
  if(polarFlag && isDiscrete && isGamma && !isLongLived) {
    SampleDirection(nucleus, mpRatio, JP1, JP2, MP);
  } else {
    fDirection = G4RandomDirection();
  }

  G4double emass = part->GetPDGMass();

  // 2-body decay in rest frame
  G4double ecm       = lv.mag();
  G4ThreeVector bst  = lv.boostVector();

  //G4cout << "Ecm= " << ecm << " mass= " << mass << " emass= " << emass 
  //	 << "  isLongLived: " << isLongLived << G4endl;

  G4double energy = 0.5*((ecm - mass)*(ecm + mass) + emass*emass)/ecm;
  energy = std::max(energy, emass);

  G4double mom = std::sqrt((energy - emass)*(energy + emass));

  G4LorentzVector res4mom(mom * fDirection.x(),
			  mom * fDirection.y(),
			  mom * fDirection.z(), energy);

  // Lab system transform for short lived level
  if(!isLongLived) { 
    res4mom.boost(bst); 
    lv -= res4mom;
  } else {  
    // In exceptional case sample decay at rest at not correct position 
    // of stopping ion, 4-momentum balance is breaked but gamma energy
    // is correct
    lv -= res4mom;
    G4double E = lv.e();
    G4double P2= (E - mass)*(E + mass);
    G4ThreeVector v = lv.vect().unit();
    G4double p = 0.0;
    if(P2 > 0.0) { p = std::sqrt(P2); } 
    else { E = mass; }
    lv.set(v.x()*p, v.y()*p, v.z()*p, E);  
  }

  // modified primary fragment 
  nucleus->SetMomentum(lv);

  // gamma or e- are produced
  result = new G4Fragment(res4mom, part);

  //  G4cout << " DeltaE= " << e0 - lv.e() - res4mom.e() << G4endl;
  //G4cout << "G4GammaTransition::GenerateGamma : " << thePhoton << G4endl;
  //G4cout << "       Left nucleus: " << aNucleus << G4endl;
  return result;
}

void G4GammaTransition::SampleDirection(G4Fragment* nuc, G4double ratio, 
					G4int twoJ1, G4int twoJ2, G4int mp)
{
  //G4cout << "G4GammaTransition::SampleDirection" << G4endl;
  G4NuclearPolarization* np = nuc->GetNuclearPolarization();

  // initial state is non-polarized - create polarization
  if(!np) {
    np = new G4NuclearPolarization();
    nuc->SetNuclearPolarization(np);
    fDirection = G4RandomDirection();
    fPolTrans.UpdatePolarizationToFinalState(fDirection.z(), fDirection.phi(), nuc);
    return;
  }

  // PhotonEvaporation dataset:
  // The multipolarity number with 1,2,3,4,5,6,7 representing E1,M1,E2,M2,E3,M3
  // monopole transition and 100*Nx+Ny representing multipolarity transition with
  // Ny and Ny taking the value 1,2,3,4,5,6,7 referring to E1,M1,E2,M2,E3,M3,..
  // For example a M1+E2 transition would be written 203.
  //
  // In M1+E2, M1 is the primary transition (L) and E2 is the secondary (L')
  // So mp = 203 means L0 = 1, Lp = 2 while mp = 2 means L0 = 1, Lp = 0
  G4int L0 = 0, Lp = 0;
  if(mp >= 100) {
    L0 = mp / 200;
    Lp = (mp - 200*L0)/2;

  } else {
    L0 = mp / 2;
    //if(ratio != 0.0) {
    //  G4cout << "Warning: Got ratio = " << ratio 
    //         << " when 0 was expected... Setting to zero." << G4endl;
    ratio = 0.;
  }

  fPolTrans.SetGammaTransitionData(twoJ1, twoJ2, L0, ratio, Lp);

  G4double cosTheta = fPolTrans.GenerateGammaCosTheta(np->GetPolarization());
  G4double phi = fPolTrans.GenerateGammaPhi(cosTheta, np->GetPolarization());
  fPolTrans.UpdatePolarizationToFinalState(cosTheta, phi, nuc);

  G4double sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));
  fDirection.set(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
}
