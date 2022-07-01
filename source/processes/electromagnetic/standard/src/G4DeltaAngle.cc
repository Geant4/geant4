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
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"

using namespace std;

G4DeltaAngle::G4DeltaAngle(const G4String&)
  : G4VEmAngularDistribution("deltaVI")
{
  fElectron = G4Electron::Electron();
  nprob = 26;
  fShellIdx = -1;
  prob.resize(nprob,0.0);
}    

G4DeltaAngle::~G4DeltaAngle() = default;

G4ThreeVector& 
G4DeltaAngle::SampleDirectionForShell(const G4DynamicParticle* dp,
                              G4double kinEnergyFinal, G4int Z, G4int idx, 
                              const G4Material* mat)
{
  fShellIdx = idx;
  return SampleDirection(dp, kinEnergyFinal,Z, mat);
}

G4ThreeVector& 
G4DeltaAngle::SampleDirection(const G4DynamicParticle* dp,
                              G4double kinEnergyFinal, G4int Z, 
                              const G4Material*)
{
  G4int nShells = G4AtomicShells::GetNumberOfShells(Z);
  G4int idx = fShellIdx;

  // if idx is not properly defined sample shell index
  if(idx < 0 || idx >= nShells) {
    if(nShells> nprob) {
      nprob = nShells;
      prob.resize(nprob,0.0);
    }
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
  }
  G4double bindingEnergy = G4AtomicShells::GetBindingEnergy(Z, idx);
  G4double cost;
  /*
  G4cout << "E(keV)= " << kinEnergyFinal/keV 
         << " Ebind(keV)= " << bindingEnergy
         << " idx= " << idx << " nShells= " << nShells << G4endl;
  */
  G4int n = 0;
  G4bool isOK = false;
  static const G4int nmax = 100;
  do {
    ++n;
    // the atomic electron
    G4double x = -G4Log(G4UniformRand());
    G4double eKinEnergy = bindingEnergy*x;
    G4double ePotEnergy = bindingEnergy*(1.0 + x);
    G4double e = kinEnergyFinal + ePotEnergy + electron_mass_c2;
    G4double p = sqrt((e + electron_mass_c2)*(e - electron_mass_c2));

    G4double totEnergy = dp->GetTotalEnergy();
    G4double totMomentum = dp->GetTotalMomentum();
    if(dp->GetParticleDefinition() == fElectron) {
      totEnergy += ePotEnergy;
      totMomentum = sqrt((totEnergy + electron_mass_c2)
                         *(totEnergy - electron_mass_c2));
    }
 
    G4double eTotEnergy = eKinEnergy + electron_mass_c2;
    G4double eTotMomentum = sqrt(eKinEnergy*(eTotEnergy + electron_mass_c2));
    G4double costet = 2*G4UniformRand() - 1;
    G4double sintet = sqrt((1 - costet)*(1 + costet));

    cost = 1.0;
    if(n >= nmax) { 
      /*
      G4ExceptionDescription ed;
      ed << "### G4DeltaAngle Warning: " << n 
         << " iterations - stop the loop with cost= 1.0 " 
         << " for " << dp->GetDefinition()->GetParticleName() << "\n" 
         << " Ekin(MeV)= " << dp->GetKineticEnergy()/MeV 
         << " Efinal(MeV)= " << kinEnergyFinal/MeV 
         << " Ebinding(MeV)= " << bindingEnergy/MeV; 
      G4Exception("G4DeltaAngle::SampleDirection","em0044",
                  JustWarning, ed,"");
      */
      if(0.0 ==  bindingEnergy) { isOK = true; }
      bindingEnergy = 0.0; 
    } 

    G4double x0 = p*(totMomentum + eTotMomentum*costet);
    /*
    G4cout << " x0= " << x0 << " p= " << p 
           << "  ptot= " << totMomentum << " pe= " <<  eTotMomentum
           << " e= " << e << " totMom= " <<  totMomentum
           << G4endl;
    */
    if(x0 > 0.0) {
      G4double x1 = p*eTotMomentum*sintet;
      G4double x2 = totEnergy*(eTotEnergy - e) - e*eTotEnergy 
        - totMomentum*eTotMomentum*costet + electron_mass_c2*electron_mass_c2;
      G4double y = -x2/x0;
      if(std::abs(y) <= 1.0) { 
        cost = -(x2 + x1*sqrt(1. - y*y))/x0; 
        if(std::abs(cost) <= 1.0) { isOK = true; }
        else { cost = 1.0; }
      }

      /*
      G4cout << " Ekin(MeV)= " << dp->GetKineticEnergy() 
             << " e1(keV)= " <<  eKinEnergy/keV 
             << " e2(keV)= " << (e - electron_mass_c2)/keV 
             << " 1-cost= " << 1-cost 
             << " x0= " << x0 << " x1= " << x1 << " x2= " << x2 
             << G4endl; 
      */
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(!isOK);

  G4double sint = sqrt((1 - cost)*(1 + cost));
  G4double phi  = twopi*G4UniformRand(); 

  fLocalDirection.set(sint*cos(phi), sint*sin(phi), cost);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}
