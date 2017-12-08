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
//  File:   G4ECDecay.cc                                                      //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   25 November 2014                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4ECDecay.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShells.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4ECDecay::G4ECDecay(const G4ParticleDefinition* theParentNucleus,
                     const G4double& branch, const G4double& Qvalue,
                     const G4double& excitationE,
                     const G4Ions::G4FloatLevelBase& flb,
                     const G4RadioactiveDecayMode& mode)
 : G4NuclearDecay("electron capture", mode, excitationE, flb), transitionQ(Qvalue),
   applyARM(true)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent 
  SetBR(branch);

  SetNumberOfDaughters(2);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4int daughterZ = theParentNucleus->GetAtomicNumber() - 1;
  G4int daughterA = theParentNucleus->GetAtomicMass(); 
  SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, excitationE, flb) );
  SetDaughter(1, "nu_e");
}


G4ECDecay::~G4ECDecay()
{}


G4DecayProducts* G4ECDecay::DecayIt(G4double)
{
  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Fill G4MT_daughters with alpha and residual nucleus (stored by SetDaughter)  
  CheckAndFillDaughters();

  // Get shell number of captured electron
  G4int shellIndex = -1;
  switch (theMode)
    {
    case KshellEC:
      shellIndex = 0;
      break;
    case LshellEC: //default PL1/PL=0.995 PL2/PL=0.005 from Fe55 Allowed transition
      if (G4UniformRand() < 0.995) shellIndex =1;
      else shellIndex =2;
      break;
    case MshellEC: //default PM1/PM=0.995 PM2/PM=0.005 from Fe55 Allowed transition
      if (G4UniformRand() < 0.995) shellIndex =4;
      else shellIndex =5;
      break;
    default:
      G4Exception("G4ECDecay::DecayIt()", "HAD_RDM_009",
                  FatalException, "Invalid electron shell selected");
    }

  // Initialize decay products with parent nucleus at rest
  G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);
  G4double eBind = 0.0;

  // G4LossTableManager must already be initialized with G4UAtomicDeexcitation
  // This is currently done in G4RadioactiveDecay::BuildPhysicsTable
  G4VAtomDeexcitation* atomDeex = 
          G4LossTableManager::Instance()->AtomDeexcitation();
  std::vector<G4DynamicParticle*> armProducts;

  if (applyARM) {
    if (atomDeex) {
      G4int aZ = G4MT_daughters[0]->GetAtomicNumber();
      G4int nShells = G4AtomicShells::GetNumberOfShells(aZ);
      if (shellIndex >= nShells) shellIndex = nShells;
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIndex);
      const G4AtomicShell* shell = atomDeex->GetAtomicShell(aZ, as);
      eBind = shell->BindingEnergy(); 
      if (atomDeex->IsFluoActive() && aZ > 5 && aZ < 100) {
        // Do atomic relaxation
          // VI, SI
          // Allows fixing of Bugzilla 1727
          //const G4double deexLimit = 0.1*keV;
          G4double deexLimit = 0.1*keV;
          if (G4EmParameters::Instance()->DeexcitationIgnoreCut())  deexLimit =0.;
          //
        atomDeex->GenerateParticles(&armProducts, shell, aZ, deexLimit, deexLimit);
      }

      G4double productEnergy = 0.;
      for (G4int i = 0; i < G4int(armProducts.size()); i++)
        productEnergy += armProducts[i]->GetKineticEnergy();

      G4double deficit = shell->BindingEnergy() - productEnergy;
      if (deficit > 0.0) {
        // Add a dummy electron to make up extra energy
        G4double cosTh = 1.-2.*G4UniformRand();
        G4double sinTh = std::sqrt(1.- cosTh*cosTh);
        G4double phi = twopi*G4UniformRand();

        G4ThreeVector electronDirection(sinTh*std::sin(phi),
                                        sinTh*std::cos(phi), cosTh);
        G4DynamicParticle* extra =
          new G4DynamicParticle(G4Electron::Electron(), electronDirection,
                                deficit);
        armProducts.push_back(extra);
      }
    } // atomDeex
  }  // applyARM

  G4double daughterMass = G4MT_daughters[0]->GetPDGMass();

  // CM momentum using Q value corrected for binding energy of captured electron
  G4double Q = transitionQ - eBind; 
  G4double cmMomentum = Q*(Q + 2.*daughterMass)/(Q + daughterMass)/2.;

  G4double costheta = 2.*G4UniformRand() - 1.0;
  G4double sintheta = std::sqrt(1.0 - costheta*costheta);
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                          costheta);
  G4double KE = cmMomentum;
  G4DynamicParticle* daughterParticle =
    new G4DynamicParticle(G4MT_daughters[1], direction, KE, 0.0);
  products->PushProducts(daughterParticle);

  KE = std::sqrt(cmMomentum*cmMomentum + daughterMass*daughterMass) - daughterMass;
  daughterParticle =
    new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, daughterMass);
  products->PushProducts(daughterParticle);

  G4int nArm = armProducts.size();
  if (nArm > 0) {
    G4ThreeVector bst = daughterParticle->Get4Momentum().boostVector();
    for (G4int i = 0; i < nArm; ++i) {
      G4DynamicParticle* dp = armProducts[i];
      G4LorentzVector lv = dp->Get4Momentum().boost(bst);
      dp->Set4Momentum(lv);
      products->PushProducts(dp);
    }
  }

  // Energy conservation check
  /*
  G4int newSize = products->entries();
  G4DynamicParticle* temp = 0;
  G4double KEsum = 0.0;
  for (G4int i = 0; i < newSize; i++) {
    temp = products->operator[](i);
    KEsum += temp->GetKineticEnergy();
  }

  G4double eCons = (transitionQ - KEsum)/keV; 
  G4cout << " EC check: Ediff (keV) = " << eCons << G4endl; 
  */
  return products;
}


void G4ECDecay::DumpNuclearInfo()
{
  G4cout << " G4ECDecay of parent nucleus " << GetParentName() << " from ";
  if (theMode == 3) {
    G4cout << "K shell";
  } else if (theMode == 4) {
    G4cout << "L shell";
  } else if (theMode == 5) {
    G4cout << "M shell";
  }
  G4cout << G4endl;
  G4cout << " to " << GetDaughterName(0) << " + " << GetDaughterName(1)
         << " with branching ratio " << GetBR() << "% and Q value "
         << transitionQ << G4endl;
}

