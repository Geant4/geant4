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
//  File:   G4SFDecay.cc                                                      //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   27 March 2019                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4SFDecay.hh"
#include "G4fissionEvent.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


G4SFDecay::G4SFDecay(const G4ParticleDefinition* theParentNucleus,
                     const G4double& branch, const G4double& Qvalue,
                     const G4double& excitationE, 
                     const G4Ions::G4FloatLevelBase& flb)
 : G4NuclearDecay("SF decay", SpFission, excitationE, flb), transitionQ(Qvalue)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent 
  SetBR(branch);

  parentZ = theParentNucleus->GetAtomicNumber();
  parentA = theParentNucleus->GetAtomicMass(); 

//  SetNumberOfDaughters(0);  doesn't seem to work
//  Set to 1 for now
//  When we later install proper fragment generation, set this to 2
//  and modify Z and A in DecayIt() 
  SetNumberOfDaughters(1);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  SetDaughter(0, theIonTable->GetIon(parentZ, parentA, excitationE, flb) );
}


G4SFDecay::~G4SFDecay()
{}


G4DecayProducts* G4SFDecay::DecayIt(G4double)
{
  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Set up final state
  // parentParticle is set at rest here because boost with correct momentum 
  // is done later
  G4LorentzVector atRest(G4MT_parent->GetPDGMass(),
                         G4ThreeVector(0.,0.,0.) );
  G4DynamicParticle parentParticle(G4MT_parent, atRest);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  // Generate LLNL code for parent nucleus
  G4int code = 1000*G4MT_parent->GetAtomicNumber() +
                    G4MT_parent->GetAtomicMass();
  
  // Let G4fissionEvent do the decay
  // Argument 2 : time, passed through to set time of secondary
  // Argument 3 : nubar, when set to -1, selects spontaneous fission
  // Argument 4 : energy used in Watt spectrum 
  G4fissionEvent aFission(code, 10.0, -1, 0.);
  G4int nNeut = aFission.getNeutronNu();
  G4int nGam = aFission.getPhotonNu();

  G4DynamicParticle* dynPart = 0; 
  G4ThreeVector direction(0.,0.,0.);
  G4double KE = 0.;   // energy from generator is in MeV 

  if (nNeut > 0) {
    // Put neutrons on the stack
    for (G4int i = 0; i < nNeut; i++) {
      KE = aFission.getNeutronEnergy(i);
      direction.setX(aFission.getNeutronDircosu(i) );
      direction.setY(aFission.getNeutronDircosv(i) );
      direction.setZ(aFission.getNeutronDircosw(i) );

      dynPart = new G4DynamicParticle(G4Neutron::Neutron(), direction, KE);
    
      // Not needed - time comes from parent track
      // dynPart->SetProperTime(aFission.getNeutronAge(i) );
  
      products->PushProducts(dynPart);
    }

    // Put gammas on the stack
    for (G4int i = 0; i < nGam; i++) {
      KE = aFission.getPhotonEnergy(i);
      direction.setX(aFission.getPhotonDircosu(i) );
      direction.setY(aFission.getPhotonDircosv(i) );
      direction.setZ(aFission.getPhotonDircosw(i) );

      dynPart = new G4DynamicParticle(G4Gamma::Gamma(), direction, KE);

      // Not needed - time comes from parent track
      // dynPart->SetProperTime(aFission.getPhotonAge(i) );

      products->PushProducts(dynPart);
      // No residual nucleus in this model
    }

  } else {
    // No data for this isotope, return parent nucleus for now
//    G4cout << " No fission data for this isotope " << G4endl;
    G4DynamicParticle* parent =
      new G4DynamicParticle(G4MT_parent, G4ThreeVector(0.,0.,0.) );
    products->PushProducts(parent);
  } 

  // Energy conservation check not valid in this model

  return products;
}


void G4SFDecay::DumpNuclearInfo()
{
  G4cout << " G4SFDecay for parent nucleus " << GetParentName() << G4endl;
  G4cout << " decays to neutrons and gammas, with branching ratio " << GetBR()
         << "% and Q value " << transitionQ << G4endl;
}

