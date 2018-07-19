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
//  File:   G4BetaMinusDecay.cc                                               //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   25 October 2014                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4BetaMinusDecay.hh"
#include "G4BetaDecayCorrections.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>

G4BetaMinusDecay::G4BetaMinusDecay(const G4ParticleDefinition* theParentNucleus,
                                   const G4double& branch, const G4double& e0,
                                   const G4double& excitationE,
                                   const G4Ions::G4FloatLevelBase& flb,
                                   const G4BetaDecayType& betaType)
 : G4NuclearDecay("beta- decay", BetaMinus, excitationE, flb), endpointEnergy(e0)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent 
  SetBR(branch);

  SetNumberOfDaughters(3);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4int daughterZ = theParentNucleus->GetAtomicNumber() + 1;
  G4int daughterA = theParentNucleus->GetAtomicMass();
  SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, excitationE, flb) );
  SetDaughter(1, "e-");
  SetDaughter(2, "anti_nu_e");

  SetUpBetaSpectrumSampler(daughterZ, daughterA, betaType);
}


G4BetaMinusDecay::~G4BetaMinusDecay()
{
  delete spectrumSampler;
}


G4DecayProducts* G4BetaMinusDecay::DecayIt(G4double)
{
  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Fill G4MT_daughters with e-, nu and residual nucleus (stored by SetDaughter)  
  CheckAndFillDaughters();

  G4double parentMass = G4MT_parent->GetPDGMass();
  G4double eMass = G4MT_daughters[1]->GetPDGMass();
  G4double nucleusMass = G4MT_daughters[0]->GetPDGMass();
  // Set up final state
  // parentParticle is set at rest here because boost with correct momentum 
  // is done later
  G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  if (spectrumSampler) {
    // Electron, neutrino and daughter nucleus energies
    G4double eKE = endpointEnergy*spectrumSampler->shoot(G4Random::getTheEngine() );
    G4double eMomentum = std::sqrt(eKE*(eKE + 2.*eMass) );

    G4double cosThetaENu = 2.*G4UniformRand() - 1.;
    G4double eTE = eMass + eKE;
    G4double nuEnergy = ((endpointEnergy - eKE)*(parentMass + nucleusMass - eTE)
            - eMomentum*eMomentum)/(parentMass - eTE + eMomentum*cosThetaENu)/2.;

    // Electron 4-vector, isotropic angular distribution
    G4double cosTheta = 2.*G4UniformRand() - 1.0;
    G4double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);

    G4double phi = twopi*G4UniformRand()*rad;
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);

    G4ParticleMomentum eDirection(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
    G4DynamicParticle* dynamicElectron
      = new G4DynamicParticle(G4MT_daughters[1], eDirection*eMomentum);
    products->PushProducts(dynamicElectron);

    // Neutrino 4-vector
    G4double sinThetaENu = std::sqrt(1.0 - cosThetaENu*cosThetaENu);
    phi = twopi*G4UniformRand()*rad;
    G4double sinPhiNu = std::sin(phi);
    G4double cosPhiNu = std::cos(phi);

    G4ParticleMomentum nuDirection;
    nuDirection.setX(sinThetaENu*cosPhiNu*cosTheta*cosPhi -
                     sinThetaENu*sinPhiNu*sinPhi + cosThetaENu*sinTheta*cosPhi);
    nuDirection.setY(sinThetaENu*cosPhiNu*cosTheta*sinPhi +
                     sinThetaENu*sinPhiNu*cosPhi + cosThetaENu*sinTheta*sinPhi);
    nuDirection.setZ(-sinThetaENu*cosPhiNu*sinTheta + cosThetaENu*cosTheta);

    G4DynamicParticle* dynamicNeutrino
      = new G4DynamicParticle(G4MT_daughters[2], nuDirection*nuEnergy);
    products->PushProducts(dynamicNeutrino);

    // Daughter nucleus 4-vector
    // p_D = - p_e - p_nu
    G4DynamicParticle* dynamicDaughter =
      new G4DynamicParticle(G4MT_daughters[0],
                            -eDirection*eMomentum - nuDirection*nuEnergy);
    products->PushProducts(dynamicDaughter);

  } else {
    // electron energy below threshold -> no decay
    G4DynamicParticle* noDecay =
      new G4DynamicParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
    products->PushProducts(noDecay);
  }

  // Check energy conservation against Q value, not nuclear masses
  /*
  G4int nProd = products->entries();
  G4DynamicParticle* temp = 0;
  G4double Esum = 0.0;
  for (G4int i = 0; i < nProd; i++) {
    temp = products->operator[](i);
    // G4cout << temp->GetParticleDefinition()->GetParticleName() << " has " 
    //        << temp->GetTotalEnergy()/keV << " keV " << G4endl;
    Esum += temp->GetKineticEnergy();
  }
  G4double eCons = (endpointEnergy - Esum)/keV;
  if (std::abs(eCons) > 0.001) G4cout << " Beta- check: eCons = " << eCons << G4endl;
  */
  return products;
}


void
G4BetaMinusDecay::SetUpBetaSpectrumSampler(const G4int& daughterZ,
                                           const G4int& daughterA,
                                           const G4BetaDecayType& betaType)
{
  G4double e0 = endpointEnergy/CLHEP::electron_mass_c2;
  G4BetaDecayCorrections corrections(daughterZ, daughterA);
  spectrumSampler = 0;

  if (e0 > 0) {
    // Array to store spectrum pdf
    G4int npti = 100;
    G4double* pdf = new G4double[npti];

    G4double e;  // Total electron energy in units of electron mass
    G4double p;  // Electron momentum in units of electron mass
    G4double f;  // Spectral shape function
    for (G4int ptn = 0; ptn < npti; ptn++) {
      // Calculate simple phase space
      e = 1. + e0*(G4double(ptn) + 0.5)/G4double(npti);
      p = std::sqrt(e*e - 1.);
      f = p*e*(e0 - e + 1.)*(e0 - e + 1.);

      // Apply Fermi factor to get allowed shape
      f *= corrections.FermiFunction(e);

      // Apply shape factor for forbidden transitions
      f *= corrections.ShapeFactor(betaType, p, e0-e+1.);
      pdf[ptn] = f;
    }
    spectrumSampler = new G4RandGeneral(pdf, npti);
    delete[] pdf;
  }
}


void G4BetaMinusDecay::DumpNuclearInfo()
{
  G4cout << " G4BetaMinusDecay for parent nucleus " << GetParentName() << G4endl;
  G4cout << " decays to " << GetDaughterName(0) << " , " << GetDaughterName(1) 
         << " and " << GetDaughterName(2) << " with branching ratio " << GetBR()
         << "% and endpoint energy " << endpointEnergy/keV << " keV " << G4endl;
}

