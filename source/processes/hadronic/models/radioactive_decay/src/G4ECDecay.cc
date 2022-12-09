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
  DefineSubshellProbabilities(daughterZ, daughterZ);
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
  G4double ran;
  switch (theMode)
    {
    case KshellEC:
      shellIndex = 0;
      break;
    case LshellEC: // PL1+PL2+PL3=1
      ran=G4UniformRand();
      if (ran <= PL1) shellIndex =1;
      else if (ran<= (PL1+PL2)) shellIndex =2;
      else shellIndex =3;
      break;
    case MshellEC:  // PM1+PM2+PM3=1
      ran=G4UniformRand();
      if (ran < PM1) shellIndex =4;
      else if (ran< (PM1+PM2)) shellIndex =5;
      else shellIndex = 6;
      break;
    case NshellEC: // PN1+PN2+PN3=1
      ran=G4UniformRand();
      if (ran < PN1) shellIndex = 9;
      else if (ran<= (PN1+PN2)) shellIndex =2;
      else shellIndex =10;
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
    if (nullptr != atomDeex) {
      G4int aZ = G4MT_daughters[0]->GetAtomicNumber();
      G4int nShells = G4AtomicShells::GetNumberOfShells(aZ);
      if (shellIndex >= nShells) shellIndex = nShells;
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIndex);
      const G4AtomicShell* shell = atomDeex->GetAtomicShell(aZ, as);
      eBind = shell->BindingEnergy(); 
      if (atomDeex->IsFluoActive() && aZ > 5 && aZ < 105) {
        // Do atomic relaxation
	// VI, SI
	// Allows fixing of Bugzilla 1727
	//const G4double deexLimit = 0.1*keV;
	G4double deexLimit = 0.1*keV;
	if (G4EmParameters::Instance()->DeexcitationIgnoreCut()) deexLimit = 0.;

        atomDeex->GenerateParticles(&armProducts, shell, aZ, deexLimit, deexLimit);
      }

      G4double productEnergy = 0.;
      for (std::size_t i = 0; i < armProducts.size(); ++i) {
        productEnergy += armProducts[i]->GetKineticEnergy();
      }
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

  // Negative transitionQ values for some rare nuclides require this
  // Absolute values in these cases are small 
  if (Q < 0.0) Q = 0.0;

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

  std::size_t nArm = armProducts.size();
  if (nArm > 0) {
    G4ThreeVector bst = daughterParticle->Get4Momentum().boostVector();
    for (std::size_t i = 0; i < nArm; ++i) {
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
  else if (theMode == 6) {
    G4cout << "N shell";
  }
  G4cout << G4endl;
  G4cout << " to " << GetDaughterName(0) << " + " << GetDaughterName(1)
         << " with branching ratio " << GetBR() << "% and Q value "
         << transitionQ << G4endl;
}
void G4ECDecay::DefineSubshellProbabilities(G4int Z, G4int )
{ //Implementation for the case of allowed transitions
  //PL1+PL2=1. , PM1+PM2=1., PN1+PN2=1.
	PL1 = 1./(1+PL2overPL1[Z-1]);
	PL2 = PL1*PL2overPL1[Z-1];
	PM1 = 1./(1+PM2overPM1[Z-1]);
	PM2 = PM1*PM2overPM1[Z-1];
	PN1 = 1./(1+PN2overPN1[Z-1]);
	PN2 = PN1*PN2overPN1[Z-1];
}
////////////////////////////////////////////////////////////////////////////
// Table of subshell ratio probability PL2/PL1  in function of Z
// PL2/PL1 = (fL2/gL1)^2
// with gL1 and fL2 the bound electron radial wave amplitudes taken from
//            Bambynek et al., Rev. Modern Physics, vol. 49, 1971, table IX
// For Z=18 interpolation  between Z=17 and Z=19 to avoid a jump in PL2/Pl1
////////////////////////////////////////////////////////////////////////////
const G4double G4ECDecay::PL2overPL1[100] = {
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   1.8722e-04,
2.6438e-04,   3.5456e-04,   4.5790e-04,   6.1588e-04,   7.8944e-04,   9.8530e-04,   1.2030e-03,
1.4361e-03,   1.6886e-03,   1.9609e-03,   2.2641e-03,   2.5674e-03,   2.9019e-03,   3.2577e-03,
3.6338e-03,   4.0310e-03,   4.4541e-03,   4.8943e-03,   5.3620e-03,   5.8523e-03,   6.3650e-03,
6.9061e-03,   7.4607e-03,   8.0398e-03,   8.6417e-03,   9.2665e-03,   9.9150e-03,   1.0588e-02,
1.1284e-02,   1.2004e-02,   1.2744e-02,   1.3518e-02,   1.4312e-02,   1.5136e-02,   1.5981e-02,
1.6857e-02,   1.7764e-02,   1.8696e-02,   1.9682e-02,   2.0642e-02,   2.1661e-02,   2.2708e-02,
2.3788e-02,   2.4896e-02,   2.6036e-02,   2.7217e-02,   2.8409e-02,   2.9646e-02,   3.0917e-02,
3.2220e-02,   3.3561e-02,   3.4937e-02,   3.6353e-02,   3.7805e-02,   3.9297e-02,   4.0826e-02,
4.2399e-02,   4.4010e-02,   4.5668e-02,   4.7368e-02,   4.9115e-02,   5.0896e-02,   5.2744e-02,
5.4625e-02,   5.6565e-02,   5.8547e-02,   6.0593e-02,   6.2690e-02,   6.4844e-02,   6.7068e-02,
6.9336e-02,   7.1667e-02,   7.4075e-02,   7.6544e-02,   7.9085e-02,   8.1688e-02,   8.4371e-02,
8.7135e-02,   8.9995e-02,   9.2919e-02,   9.5949e-02,   9.9036e-02,   1.0226e-01,   1.0555e-01,
1.0899e-01,   1.1249e-01,   1.1613e-01,   1.1989e-01,   1.2379e-01,   1.2780e-01,   1.3196e-01,
1.3627e-01,   1.4071e-01};
////////////////////////////////////////////////////////////////////////////
// Table of subshell ratio probability PM2/PM1  in function of Z
// PM2/PM1 = (fM2/gM1)^2
// with gM1 and fM2 the bound electron radial wave amplitudes taken from
//            Bambynek et al., Rev. Modern Physics, vol. 49, 1971, table IX
////////////////////////////////////////////////////////////////////////////
const G4double G4ECDecay::PM2overPM1[100] = {
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
1.0210e-03,   1.2641e-03,   1.5231e-03,   1.7990e-03,   2.1938e-03,   2.5863e-03,   2.9621e-03,
3.3637e-03,   3.7909e-03,   4.2049e-03,   4.7021e-03,   5.1791e-03,   5.6766e-03,   6.1952e-03,
6.7045e-03,   7.2997e-03,   7.9438e-03,   8.6271e-03,   9.3294e-03,   1.0058e-02,   1.0813e-02,
1.1594e-02,   1.2408e-02,   1.3244e-02,   1.4118e-02,   1.5023e-02,   1.5962e-02,   1.6919e-02,
1.7910e-02,   1.8934e-02,   1.9986e-02,   2.1072e-02,   2.2186e-02,   2.3336e-02,   2.4524e-02,
2.5750e-02,   2.7006e-02,   2.8302e-02,   2.9629e-02,   3.0994e-02,   3.2399e-02,   3.3845e-02,
3.5328e-02,   3.6852e-02,   3.8414e-02,   4.0025e-02,   4.1673e-02,   4.3368e-02,   4.5123e-02,
4.6909e-02,   4.8767e-02,   5.0662e-02,   5.2612e-02,   5.4612e-02,   5.6662e-02,   5.8773e-02,
6.0930e-02,   6.3141e-02,   6.5413e-02,   6.7752e-02,   7.0139e-02,   7.2603e-02,   7.5127e-02,
7.7721e-02,   8.0408e-02,   8.3128e-02,   8.5949e-02,   8.8843e-02,   9.1824e-02,   9.4888e-02,
9.8025e-02,   1.0130e-01,   1.0463e-01,   1.0806e-01,   1.1159e-01,   1.1526e-01,   1.1900e-01,
1.2290e-01,   1.2688e-01,   1.3101e-01,   1.3528e-01,   1.3969e-01,   1.4425e-01,   1.4896e-01,
1.5384e-01,   1.5887e-01};
////////////////////////////////////////////////////////////////////////////
// Table of subshell ratio probability PN2/PN1  in function of Z
// PN2/PN1 = (fN2/gN1)^2
// with gN1 and fN2 are the bound electron radial wave amplitude taken from
//            Bambynek et al., Rev. Modern Physics, vol. 49, 1971, table IX
// For Z=44 interpolation  between Z=43 and Z=45 to avoid a jump in PN2/PN1
////////////////////////////////////////////////////////////////////////////
const G4double G4ECDecay::PN2overPN1[100] = {
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,   0.0000e+00,
0.0000e+00,   9.6988e-03,   1.0797e-02,   1.1706e-02,   1.2603e-02,   1.3408e-02,   1.4352e-02,
1.5511e-02,   1.6579e-02,   1.7646e-02,   1.8731e-02,   1.9886e-02,   2.1069e-02,   2.2359e-02,
2.3710e-02,   2.5058e-02,   2.6438e-02,   2.7843e-02,   2.9283e-02,   3.0762e-02,   3.2275e-02,
3.3843e-02,   3.5377e-02,   3.6886e-02,   3.8502e-02,   4.0159e-02,   4.1867e-02,   4.3617e-02,
4.5470e-02,   4.7247e-02,   4.9138e-02,   5.1065e-02,   5.3049e-02,   5.5085e-02,   5.7173e-02,
5.9366e-02,   6.1800e-02,   6.3945e-02,   6.6333e-02,   6.8785e-02,   7.1303e-02,   7.3801e-02,
7.6538e-02,   7.9276e-02,   8.2070e-02,   8.4959e-02,   8.7940e-02,   9.0990e-02,   9.4124e-02,
9.7337e-02,   1.0069e-01,   1.0410e-01,   1.0761e-01,   1.1122e-01,   1.1499e-01,   1.1882e-01,
1.2282e-01,   1.2709e-01,   1.3114e-01,   1.3549e-01,   1.4001e-01,   1.4465e-01,   1.4946e-01,
1.5443e-01,   1.5954e-01};

