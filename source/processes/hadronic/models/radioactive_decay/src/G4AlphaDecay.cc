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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4AlphaDecay.cc                                                   //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   20 October 2014                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4AlphaDecay.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>

G4AlphaDecay::G4AlphaDecay(const G4ParticleDefinition* theParentNucleus,
                           const G4double& branch, const G4double& Qvalue,
                           const G4double& excitationE,
                           const G4Ions::G4FloatLevelBase& flb)
 : G4NuclearDecay("alpha decay", Alpha, excitationE, flb), transitionQ(Qvalue)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent 
  SetBR(branch);

  SetNumberOfDaughters(2);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4int daughterZ = theParentNucleus->GetAtomicNumber() - 2;
  G4int daughterA = theParentNucleus->GetAtomicMass() - 4;
  SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, excitationE, flb) );
  SetDaughter(1, "alpha");
}


G4AlphaDecay::~G4AlphaDecay()
{}


G4DecayProducts* G4AlphaDecay::DecayIt(G4double)
{
  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Fill G4MT_daughters with alpha and residual nucleus (stored by SetDaughter)  
  CheckAndFillDaughters();

  G4double alphaMass = G4MT_daughters[1]->GetPDGMass();
  // Excitation energy included in PDG mass
  G4double nucleusMass = G4MT_daughters[0]->GetPDGMass();

  // Q value was calculated from atomic masses.
  // Use it to get correct alpha energy.
  G4double cmMomentum = std::sqrt(transitionQ*(transitionQ + 2.*alphaMass)*
                                (transitionQ + 2.*nucleusMass)*
                                (transitionQ + 2.*alphaMass + 2.*nucleusMass) )/
                                (transitionQ + alphaMass + nucleusMass)/2.; 

  // Set up final state
  // parentParticle is set at rest here because boost with correct momentum 
  // is done later
  G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt(1.0 - costheta*costheta);
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                          costheta);

  G4double KE = std::sqrt(cmMomentum*cmMomentum + alphaMass*alphaMass)
              - alphaMass;
  G4DynamicParticle* daughterparticle =
    new G4DynamicParticle(G4MT_daughters[1], direction, KE, alphaMass);
  products->PushProducts(daughterparticle);

  KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;
  daughterparticle =
    new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);
  products->PushProducts(daughterparticle);

  // Energy conservation check
  // For alpha decays, do final energy check against reaction Q value
  // which is well-measured using atomic mass differences.  Nuclear masses
  // should not be used since they are not usually directly measured and we
  // always decay atoms and not fully stripped nuclei.
  /*
  G4int nProd = products->entries();
  G4DynamicParticle* temp = 0;
  G4double Esum = 0.0;
  for (G4int i = 0; i < nProd; i++) {
    temp = products->operator[](i);
    Esum += temp->GetKineticEnergy();
  }
  G4double eCons = (transitionQ - Esum)/keV;
  if (eCons > 1.e-07) G4cout << " Alpha decay check: Ediff (keV) = " << eCons << G4endl;
  */
  return products;
}


void G4AlphaDecay::DumpNuclearInfo()
{
  G4cout << " G4AlphaDecay for parent nucleus " << GetParentName() << G4endl;
  G4cout << " decays to " << GetDaughterName(0) << " + " << GetDaughterName(1)
         << " with branching ratio " << GetBR() << "% and Q value "
         << transitionQ << G4endl;
}

