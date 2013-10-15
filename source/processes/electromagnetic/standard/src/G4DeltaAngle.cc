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
// $Id: G4DeltaAngle.cc 68380 2013-03-22 18:39:29Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4DeltaAngle
//
// Author:        Vladimir Ivantcheko
// 
// Creation date: 23 August 2013
//
// Modifications: 
//
// Class Description: 
//
// Delta-electron Angular Distribution Generation 
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#include "G4DeltaAngle.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4AtomicShells.hh"
#include "G4Log.hh"

using namespace std;

G4DeltaAngle::G4DeltaAngle(const G4String&)
  : G4VEmAngularDistribution("deltaVI")
{
  fElectron = G4Electron::Electron();
  nprob = 26;
  prob.resize(nprob,0.0);
}    

G4DeltaAngle::~G4DeltaAngle() 
{}

G4ThreeVector& 
G4DeltaAngle::SampleDirection(const G4DynamicParticle* dp,
			      G4double kinEnergyFinal, G4int Z, 
			      const G4Material*)
{
  G4int nShells = G4AtomicShells::GetNumberOfShells(Z);
  if(nShells> nprob) {
    nprob = nShells;
    prob.resize(nprob,0.0);
  }
  G4int idx;
  G4double sum = 0.0;
  for(idx=0; idx<nShells; ++idx) {
    sum += G4AtomicShells::GetNumberOfElectrons(Z, idx)
      /G4AtomicShells::GetBindingEnergy(Z, idx);
    prob[idx] = sum;
  }
  sum *= G4UniformRand();
  for(idx=0; idx<nShells; ++idx) {
    if(sum <= prob[idx]) { break; }
  }
  G4double bindingEnergy = G4AtomicShells::GetBindingEnergy(Z, idx);
  G4double mass = dp->GetParticleDefinition()->GetPDGMass();

  G4ThreeVector bst(0.0,0.0,0.0);
  G4double cost, en, mom;

  do {
  
    // the atomic electron
    G4double x = -G4Log(G4UniformRand());
    G4double eKinEnergy = bindingEnergy*x;
    G4double ePotEnergy = bindingEnergy*(1.0 + x);
    G4double e = kinEnergyFinal + ePotEnergy + electron_mass_c2;

    G4double totEnergy = dp->GetTotalEnergy();
    G4double totMomentum = dp->GetTotalMomentum();
    if(dp->GetParticleDefinition() == fElectron) {
      totEnergy += ePotEnergy;
      totMomentum = sqrt((totEnergy + electron_mass_c2)
			 *(totEnergy - electron_mass_c2));
    }
 
    G4double eTotMomentum = sqrt(eKinEnergy*(eKinEnergy + 2*electron_mass_c2));
    G4double phi = G4UniformRand()*twopi;
    G4double costet = 2*G4UniformRand() - 1;
    G4double sintet = sqrt((1 - costet)*(1 + costet));
 
    G4LorentzVector lv0(eTotMomentum*sintet*cos(phi),
			eTotMomentum*sintet*sin(phi),
			eTotMomentum*costet + totMomentum,
			eKinEnergy + electron_mass_c2 + totEnergy);
    bst = lv0.boostVector();

    G4double m0  = lv0.mag();
    G4double bet = lv0.beta();
    G4double gam = lv0.gamma();

    en = 0.5*(m0*m0 - mass*mass + electron_mass_c2*electron_mass_c2)/m0;
    mom = sqrt((en + electron_mass_c2)*(en - electron_mass_c2));

    cost= (e/gam - en)/(mom*bet);

    //G4cout << "e= " << e << " gam= " << gam << " en= " << en 
    //	   << " mom= " << mom << "  beta= " << bet << " cost= " << cost 
    //	   << G4endl; 

  } while(std::fabs(cost) > 1.0);

  G4double sint = sqrt((1 - cost)*(1 + cost));
  G4double phi  = twopi*G4UniformRand(); 

  G4LorentzVector lv1(sint*std::cos(phi)*mom, sint*std::sin(phi)*mom,
		      mom*cost, en);
  lv1.boost(bst);

  fLocalDirection.set(lv1.x(), lv1.y(), lv1.z());
  fLocalDirection = fLocalDirection.unit();
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

void G4DeltaAngle::PrintGeneratorInformation() const
{} 
