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
// 080721 Create adjust_final_state method by T. Koi
// 080801 Residual reconstruction with theNDLDataA,Z (A, Z, and momentum are adjusted) by T. Koi
// 101110 Set lower limit for gamma energy(1keV) by T. Koi
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//

#include "G4ParticleHPFinalState.hh"

#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Gamma.hh"
#include "G4He3.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4Triton.hh"

G4ParticleHPFinalState::G4ParticleHPFinalState()
{
  theProjectile = G4Neutron::Neutron();
  theResult.Put(nullptr);
  fManager = G4ParticleHPManager::GetInstance();
  ionTable = G4IonTable::GetIonTable();
}

G4ParticleHPFinalState::~G4ParticleHPFinalState()
{
  delete theResult.Get();
}

void G4ParticleHPFinalState::adjust_final_state(G4LorentzVector init_4p_lab)
{
  G4double minimum_energy = 1 * keV;

  if (fManager->GetDoNotAdjustFinalState()) return;

  auto nSecondaries = (G4int)theResult.Get()->GetNumberOfSecondaries();

  G4int sum_Z = 0;
  G4int sum_A = 0;
  G4int max_SecZ = 0;
  G4int max_SecA = 0;
  G4int imaxA = -1;
  for (G4int i = 0; i < nSecondaries; ++i) {
    auto ptr = theResult.Get()->GetSecondary(i)->GetParticle()->GetDefinition(); 
    sum_Z += ptr->GetAtomicNumber();
    max_SecZ = std::max(max_SecZ, ptr->GetAtomicNumber());
    sum_A += ptr->GetAtomicMass();
    max_SecA = std::max(max_SecA, ptr->GetAtomicMass());
    if (ptr->GetAtomicMass() == max_SecA)
      imaxA = i;
#ifdef G4PHPDEBUG
    if (fManager->GetDEBUG())
      G4cout << "G4ParticleHPFinalState::adjust_final_stat SECO " << i << " "
             << ptr->GetParticleName() << G4endl;
#endif
  }

  G4ParticleDefinition* resi_pd = nullptr;

  G4int baseZNew = theBaseZ;
  G4int baseANew = theBaseA;
  if (theProjectile == G4Neutron::Neutron()) {
    baseANew++;
  }
  else if (theProjectile == G4Proton::Proton()) {
    baseZNew++;
    baseANew++;
  }
  else if (theProjectile == G4Deuteron::Deuteron()) {
    baseZNew++;
    baseANew += 2;
  }
  else if (theProjectile == G4Triton::Triton()) {
    baseZNew++;
    baseANew += 3;
  }
  else if (theProjectile == G4He3::He3()) {
    baseZNew += 2;
    baseANew += 3;
  }
  else if (theProjectile == G4Alpha::Alpha()) {
    baseZNew += 2;
    baseANew += 4;
  }

#ifdef G4PHPDEBUG
  if (fManager->GetDEBUG())
    G4cout << "G4ParticleHPFinalState::adjust_final_stat  BaseZ " << baseZNew << " BaseA "
           << baseANew << " sum_Z " << sum_Z << " sum_A " << sum_A << G4endl;
#endif

  G4bool needOneMoreSec = false;
  G4ParticleDefinition* oneMoreSec_pd = nullptr;
  if (baseZNew == sum_Z && baseANew == sum_A) {
    // All secondaries are already created;
    resi_pd = theResult.Get()->GetSecondary(imaxA)->GetParticle()->GetDefinition();
  }
  else {
    if (max_SecA >= baseANew - sum_A) {
      // Most heavy secondary is interpreted as residual
      resi_pd = theResult.Get()->GetSecondary(imaxA)->GetParticle()->GetDefinition();
      needOneMoreSec = true;
    }
    else {
      // creation of residual is required
      resi_pd = ionTable->GetIon(baseZNew - sum_Z, baseANew - sum_A, 0.0);
    }

    if (needOneMoreSec) {
      if (baseZNew == sum_Z && baseANew == sum_A) {
        // In this case, one neutron is added to secondaries
        if (baseANew - sum_A > 1)
          G4cout << "More than one neutron is required for the balance of baryon number!"
                 << G4endl;
        oneMoreSec_pd = G4Neutron::Neutron();
      }
      else {
#ifdef G4PHPDEBUG
        if (fManager->GetDEBUG())
          G4cout << this << "G4ParticleHPFinalState oneMoreSec_pd Z "
                 << baseZNew << " - " << sum_Z
                 << " A " << baseANew << " - " << sum_A << " projectile "
                 << theProjectile->GetParticleName() << G4endl;
#endif
        oneMoreSec_pd =
          G4IonTable::GetIonTable()->GetIon(baseZNew - sum_Z, baseANew - sum_A, 0.0);
        if (oneMoreSec_pd == nullptr) {
	  G4ExceptionDescription ed;
          ed << "G4ParticleHPFinalState oneMoreSec_pd Z=" << baseZNew << " - " << sum_Z
                 << " A=" << baseANew << " - " << sum_A << " projectile "
                 << theProjectile->GetParticleName();
          G4Exception("G4ParticleHPFinalState:adjust_final_state", "hadr01", JustWarning,
                      ed, "No adjustment will be done!");
          return;
        }
      }
    }

    if (resi_pd == nullptr) {
      // theNDLDataZ,A has the Z and A of used NDL file
      G4int ndlZNew = theNDLDataZ;
      G4int ndlANew = theNDLDataA;
      if (theProjectile == G4Neutron::Neutron()) {
        ndlANew++;
      }
      else if (theProjectile == G4Proton::Proton()) {
        ndlZNew++;
        ndlANew++;
      }
      else if (theProjectile == G4Deuteron::Deuteron()) {
        ndlZNew++;
        ndlANew += 2;
      }
      else if (theProjectile == G4Triton::Triton()) {
        ndlZNew++;
        ndlANew += 3;
      }
      else if (theProjectile == G4He3::He3()) {
        ndlZNew += 2;
        ndlANew += 3;
      }
      else if (theProjectile == G4Alpha::Alpha()) {
        ndlZNew += 2;
        ndlANew += 4;
      }
      // theNDLDataZ,A has the Z and A of used NDL file
      if (ndlZNew == sum_Z && ndlANew == sum_A) {
        G4int dif_Z = theNDLDataZ - theBaseZ;
        G4int dif_A = theNDLDataA - theBaseA;
        resi_pd = ionTable->GetIon(max_SecZ - dif_Z, max_SecA - dif_A, 0.0);
        if (resi_pd == nullptr) {
	  G4ExceptionDescription ed;
          ed << "resi_pd Z=" << max_SecZ << " - " << dif_Z << " A="
                 << max_SecA << " - " << dif_A << " projectile " 
                 << theProjectile->GetParticleName();
          G4Exception("G4ParticleHPFinalState:adjust_final_state", "hadr02", JustWarning,
                      ed, "No adjustment will be done!");
          return;
        }

        for (G4int i = 0; i < nSecondaries; ++i) {
          auto ptr = theResult.Get()->GetSecondary(i)->GetParticle(); 
          if (ptr->GetDefinition()->GetAtomicNumber() == max_SecZ &&
              ptr->GetDefinition()->GetAtomicMass() == max_SecA)
          {
            G4ThreeVector p = ptr->GetMomentum() * resi_pd->GetPDGMass()
                / ionTable->GetIon(max_SecZ, max_SecA, 0.0)->GetPDGMass();
            ptr->SetDefinition(resi_pd);
            ptr->SetMomentum(p);
          }
        }
      }
    }
  }

  G4LorentzVector secs_4p_lab(0.0);

  auto n_sec = (G4int)theResult.Get()->GetNumberOfSecondaries();
  G4double fast = 0;
  G4double slow = 1;
  G4int ifast = 0;
  G4int islow = 0;
  G4int ires = -1;

  for (G4int i = 0; i < n_sec; ++i) {
    auto ptr = theResult.Get()->GetSecondary(i)->GetParticle();
    secs_4p_lab += ptr->Get4Momentum();

    G4double beta = 0;
    if (ptr->GetDefinition() != G4Gamma::Gamma()) {
      beta = ptr->Get4Momentum().beta();
    }
    else {
      beta = 1.0;
    }

    if (ptr->GetDefinition() == resi_pd) ires = i;

    if (slow > beta && beta != 0) {
      slow = beta;
      islow = i;
    }

    if (fast <= beta) {
      if (fast != 1) {
        fast = beta;
        ifast = i;
      }
      else {
        // fast is already photon then check E
        G4double e = ptr->Get4Momentum().e();
        if (e > theResult.Get()->GetSecondary(ifast)->GetParticle()->Get4Momentum().e()) {
          // among photons, the highest E becomes the fastest
          ifast = i;
        }
      }
    }
  }

  G4LorentzVector dif_4p = init_4p_lab - secs_4p_lab;

  G4LorentzVector p4(0);
  if (ires == -1) {
    // Create and Add Residual Nucleus
    ires = nSecondaries;
    nSecondaries += 1;

    auto res = new G4DynamicParticle(resi_pd, dif_4p.v());
    theResult.Get()->AddSecondary(res, secID);

    p4 = res->Get4Momentum();
    if (slow > p4.beta()) {
      slow = p4.beta();
      islow = ires;
    }
    dif_4p = init_4p_lab - (secs_4p_lab + p4);
  }

  if (needOneMoreSec && oneMoreSec_pd != nullptr)
  //
  // fca: this is not a fix, this is a crash avoidance...
  // fca: the baryon number is still wrong, most probably because it
  // fca: should have been decreased, but since we could not create a particle
  // fca: we just do not add it
  //
  {
    nSecondaries += 1;
    auto one = new G4DynamicParticle(oneMoreSec_pd, dif_4p.v());
    theResult.Get()->AddSecondary(one, secID);
    p4 = one->Get4Momentum();
    if (slow > p4.beta()) {
      slow = p4.beta();
      islow = nSecondaries - 1;  // Because the first is 0th, so the last becomes "nSecondaries-1"
    }
    dif_4p = init_4p_lab - (secs_4p_lab + p4);
  }

  // Which is bigger dif_p or dif_e

  if (dif_4p.v().mag() < std::abs(dif_4p.e())) {
    // Adjust p
    if (minimum_energy < dif_4p.v().mag() && dif_4p.v().mag() < 1 * MeV) {
      nSecondaries += 1;
      theResult.Get()->AddSecondary(new G4DynamicParticle(G4Gamma::Gamma(), dif_4p.v()), secID);
    }
  }
  else {
    // dif_p > dif_e
    // at first momentum
    // Move residual momentum
    p4 = theResult.Get()->GetSecondary(ires)->GetParticle()->Get4Momentum();
    theResult.Get()->GetSecondary(ires)->GetParticle()->SetMomentum(p4.v() + dif_4p.v());
    dif_4p = init_4p_lab -
      (secs_4p_lab - p4 + theResult.Get()->GetSecondary(ires)->GetParticle()->Get4Momentum());
  }

  G4double dif_e = dif_4p.e() - (dif_4p.v()).mag();

  if (dif_e > 0) {
    // create 2 gamma

    nSecondaries += 2;
    G4double e1 = (dif_4p.e() - dif_4p.v().mag()) / 2;

    if (minimum_energy < e1) {
      G4ThreeVector dir = G4RandomDirection();
      theResult.Get()->AddSecondary(new G4DynamicParticle(G4Gamma::Gamma(), e1 * dir), secID);
      theResult.Get()->AddSecondary(new G4DynamicParticle(G4Gamma::Gamma(), -e1 * dir), secID);
    }
  }
  else  // dif_e < 0
  {
    // At first reduce KE of the fastest secondary;
    auto ptr = theResult.Get()->GetSecondary(ifast)->GetParticle();
    G4double ke0 = ptr->GetKineticEnergy();
    G4ThreeVector p0 = ptr->GetMomentum();
    G4ThreeVector dir = p0.unit();

    if (ke0 + dif_e > 0) {
      ptr->SetKineticEnergy(ke0 + dif_e);
      G4ThreeVector dp = p0 - theResult.Get()->GetSecondary(ifast)->GetParticle()->GetMomentum();
      G4ThreeVector p = theResult.Get()->GetSecondary(islow)->GetParticle()->GetMomentum();
      theResult.Get()->GetSecondary(islow)->GetParticle()->SetMomentum(p + dp);
    }
  }
}
