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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 12-April-06 Enable IC electron emissions T. Koi
// 26-January-07 Add G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION flag
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
// 101203 Bugzilla/Geant4 Problem 1155 Lack of residual in some case
// 110430 Temporary solution in the case of being MF6 final state in Capture reaction (MT102)
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko July-2023 converted back
//
#include "G4NeutronHPCaptureFS.hh"

#include "G4Fragment.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4Nucleus.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhotonEvaporation.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4ReactionProduct.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

G4NeutronHPCaptureFS::G4NeutronHPCaptureFS()
{
  secID = G4PhysicsModelCatalog::GetModelID("model_NeutronHPCapture");
  hasXsec = false;
  hasExactMF6 = false;
  targetMass = 0;
}

G4HadFinalState*
G4NeutronHPCaptureFS::ApplyYourself(const G4HadProjectile& theTrack)
{
  if (theResult.Get() == nullptr) theResult.Put(new G4HadFinalState);
  theResult.Get()->Clear();

  G4int i;

  // prepare neutron
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4HadProjectile* incidentParticle = &theTrack;
  G4ReactionProduct theNeutron(theTrack.GetDefinition());
  theNeutron.SetMomentum(incidentParticle->Get4Momentum().vect());
  theNeutron.SetKineticEnergy(eKinetic);

  // Prepare target
  G4ReactionProduct theTarget;
  G4Nucleus aNucleus;
  if (targetMass < 500 * MeV)
    targetMass = G4NucleiProperties::GetNuclearMass(theBaseA, theBaseZ)
      / CLHEP::neutron_mass_c2;
  G4ThreeVector neutronVelocity = theNeutron.GetMomentum()/ CLHEP::neutron_mass_c2;
  G4double temperature = theTrack.GetMaterial()->GetTemperature();
  theTarget = aNucleus.GetBiasedThermalNucleus(targetMass, neutronVelocity, temperature);
  theTarget.SetDefinitionAndUpdateE(ionTable->GetIon(theBaseZ, theBaseA, 0.0));

  // Put neutron in nucleus rest system
  theNeutron.Lorentz(theNeutron, theTarget);
  eKinetic = theNeutron.GetKineticEnergy();

  // Sample the photons
  G4ReactionProductVector* thePhotons = nullptr;
  if (HasFSData() && !fManager->GetUseOnlyPhotoEvaporation()) {
    // NDL has final state data
    if (hasExactMF6) {
      theMF6FinalState.SetTarget(theTarget);
      theMF6FinalState.SetProjectileRP(theNeutron);
      thePhotons = theMF6FinalState.Sample(eKinetic);
    }
    else {
      thePhotons = theFinalStatePhotons.GetPhotons(eKinetic);
    }
    if (thePhotons == nullptr) {
      throw G4HadronicException(__FILE__, __LINE__,
                                "Final state data for photon is not properly allocated");
    }
  }
  else {
    // NDL does not have final state data or forced to use PhotoEvaporation model
    G4ThreeVector aCMSMomentum = theNeutron.GetMomentum() + theTarget.GetMomentum();
    G4LorentzVector p4(aCMSMomentum, theTarget.GetTotalEnergy() + theNeutron.GetTotalEnergy());
    G4Fragment nucleus(theBaseA + 1, theBaseZ, p4);
    G4PhotonEvaporation photonEvaporation;
    // T. K. add
    photonEvaporation.SetICM(true);
    G4FragmentVector* products = photonEvaporation.BreakItUp(nucleus);
    thePhotons = new G4ReactionProductVector;
    for (auto it = products->cbegin(); it != products->cend(); ++it) {
      auto theOne = new G4ReactionProduct;
      // T. K. add
      if ((*it)->GetParticleDefinition() != nullptr)
        theOne->SetDefinition((*it)->GetParticleDefinition());
      else
        theOne->SetDefinition(G4Gamma::Gamma());  // this definiion will be over writen

      // T. K. comment out below line
      if ((*it)->GetMomentum().mag() > 10 * CLHEP::MeV)
        theOne->SetDefinition(ionTable->GetIon(theBaseZ, theBaseA + 1, 0));

      if ((*it)->GetExcitationEnergy() > 1.0e-2 * eV) {
        G4double ex = (*it)->GetExcitationEnergy();
        auto aPhoton = new G4ReactionProduct;
        aPhoton->SetDefinition(G4Gamma::Gamma());
        aPhoton->SetMomentum((*it)->GetMomentum().vect().unit() * ex);
        // aPhoton->SetTotalEnergy( ex ); //will be calculated from momentum
        thePhotons->push_back(aPhoton);
      }

      theOne->SetMomentum((*it)->GetMomentum().vect()
                          * ((*it)->GetMomentum().t() - (*it)->GetExcitationEnergy())
                          / (*it)->GetMomentum().t());
      thePhotons->push_back(theOne);
      delete *it;
    }
    delete products;
  }

  // Add them to the final state
  G4int nPhotons = (G4int)thePhotons->size();

  if (!fManager->GetDoNotAdjustFinalState()) {
    // Make at least one photon
    // 101203 TK
    if (nPhotons == 0) {
      auto theOne = new G4ReactionProduct;
      theOne->SetDefinition(G4Gamma::Gamma());
      G4ThreeVector direction = G4RandomDirection();
      theOne->SetMomentum(direction);
      thePhotons->push_back(theOne);
      ++nPhotons;  // 0 -> 1
    }
    // One photon case: energy set to Q-value
    // 101203 TK
    if (nPhotons == 1 &&
	(*thePhotons)[0]->GetDefinition()->GetBaryonNumber() == 0) {
      G4ThreeVector direction = (*thePhotons)[0]->GetMomentum().unit();

      G4double Q = ionTable->GetIonMass(theBaseZ, theBaseA, 0)
	+ CLHEP::neutron_mass_c2
	- ionTable->GetIonMass(theBaseZ, theBaseA + 1, 0);
      (*thePhotons)[0]->SetMomentum(Q * direction);
    }
  }

  // back to lab system
  for (i = 0; i < nPhotons; i++) {
    (*thePhotons)[i]->Lorentz(*((*thePhotons)[i]), -1*theTarget);
  }

  // Recoil, if only one gamma
  // if (1==nPhotons)
  if (nPhotons == 1 && thePhotons->operator[](0)->GetDefinition()->GetBaryonNumber() == 0) {
    auto theOne = new G4DynamicParticle;
    G4ParticleDefinition* aRecoil = ionTable->GetIon(theBaseZ, theBaseA + 1, 0);
    theOne->SetDefinition(aRecoil);
    // Now energy;
    // Can be done slightly better @
    G4ThreeVector aMomentum = theTrack.Get4Momentum().vect() + theTarget.GetMomentum()
                              - thePhotons->operator[](0)->GetMomentum();

    theOne->SetMomentum(aMomentum);
    theResult.Get()->AddSecondary(theOne, secID);
  }

  // Now fill in the gammas.
  for (i = 0; i < nPhotons; ++i) {
    // back to lab system
    auto theOne = new G4DynamicParticle;
    theOne->SetDefinition(thePhotons->operator[](i)->GetDefinition());
    theOne->SetMomentum(thePhotons->operator[](i)->GetMomentum());
    theResult.Get()->AddSecondary(theOne, secID);
    delete thePhotons->operator[](i);
  }
  delete thePhotons;

  // 101203TK
  G4bool residual = false;
  G4ParticleDefinition* aRecoil = ionTable->GetIon(theBaseZ, theBaseA + 1, 0);
  for (std::size_t j = 0; j != theResult.Get()->GetNumberOfSecondaries(); j++) {
    if (theResult.Get()->GetSecondary(j)->GetParticle()->GetDefinition() == aRecoil)
      residual = true;
  }

  if (!residual) {
    G4int nNonZero = 0;
    G4LorentzVector p_photons(0, 0, 0, 0);
    for (std::size_t j = 0; j != theResult.Get()->GetNumberOfSecondaries(); ++j) {
      p_photons += theResult.Get()->GetSecondary(j)->GetParticle()->Get4Momentum();
      // To many 0 momentum photons -> Check PhotonDist
      if (theResult.Get()->GetSecondary(j)->GetParticle()->Get4Momentum().e() > 0) nNonZero++;
    }

    // Can we include kinetic energy here?
    G4double deltaE = (theTrack.Get4Momentum().e() + theTarget.GetTotalEnergy())
                      - (p_photons.e() + aRecoil->GetPDGMass());

    // Add photons
    if (nPhotons - nNonZero > 0) {
      // G4cout << "TKDB G4NeutronHPCaptureFS::ApplyYourself we will create additional " <<
      // nPhotons - nNonZero << " photons" << G4endl;
      std::vector<G4double> vRand;
      vRand.push_back(0.0);
      for (G4int j = 0; j != nPhotons - nNonZero - 1; j++) {
        vRand.push_back(G4UniformRand());
      }
      vRand.push_back(1.0);
      std::sort(vRand.begin(), vRand.end());

      std::vector<G4double> vEPhoton;
      for (G4int j = 0; j < (G4int)vRand.size() - 1; j++) {
        vEPhoton.push_back(deltaE * (vRand[j + 1] - vRand[j]));
      }
      std::sort(vEPhoton.begin(), vEPhoton.end());

      for (G4int j = 0; j < nPhotons - nNonZero - 1; j++) {
        // Isotopic in LAB OK?
        //  Bug # 1745 DHW G4double theta = pi*G4UniformRand();
        G4ThreeVector tempVector = G4RandomDirection()*vEPhoton[j]; 

        p_photons += G4LorentzVector(tempVector, tempVector.mag());
        auto theOne = new G4DynamicParticle;
        theOne->SetDefinition(G4Gamma::Gamma());
        theOne->SetMomentum(tempVector);
        theResult.Get()->AddSecondary(theOne, secID);
      }

      //        Add last photon
      auto theOne = new G4DynamicParticle;
      theOne->SetDefinition(G4Gamma::Gamma());
      //        For better momentum conservation
      G4ThreeVector lastPhoton = -p_photons.vect().unit() * vEPhoton.back();
      p_photons += G4LorentzVector(lastPhoton, lastPhoton.mag());
      theOne->SetMomentum(lastPhoton);
      theResult.Get()->AddSecondary(theOne, secID);
    }

    // Add residual
    auto theOne = new G4DynamicParticle;
    G4ThreeVector aMomentum =
      theTrack.Get4Momentum().vect() + theTarget.GetMomentum() - p_photons.vect();
    theOne->SetDefinition(aRecoil);
    theOne->SetMomentum(aMomentum);
    theResult.Get()->AddSecondary(theOne, secID);
  }
  // 101203TK END

  // clean up the primary neutron
  theResult.Get()->SetStatusChange(stopAndKill);
  return theResult.Get();
}

void G4NeutronHPCaptureFS::Init(G4double AA, G4double ZZ, G4int M,
                                G4String& dirName, G4String&,
                                G4ParticleDefinition*)
{
  G4int Z = G4lrint(ZZ);
  G4int A = G4lrint(AA);
  // TK110430 BEGIN
  std::stringstream ss;
  ss << Z;
  G4String sZ;
  ss >> sZ;
  ss.clear();
  ss << A;
  G4String sA;
  ss >> sA;

  ss.clear();
  G4String sM;
  if (M > 0) {
    ss << "m";
    ss << M;
    ss >> sM;
    ss.clear();
  }

  G4String element_name = theNames.GetName(Z - 1);
  G4String filenameMF6 = dirName + "/FSMF6/" + sZ + "_" + sA + sM + "_" + element_name;

  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filenameMF6, theData);

  // TK110430 Only use MF6MT102 which has exactly same A and Z
  // Even _nat_ do not select and there is no _nat_ case in ENDF-VII.0
  if (theData.good()) {
    hasExactMF6 = true;
    theMF6FinalState.Init(theData);
    return;
  }
  // TK110430 END

  G4String tString = "/FS";
  G4bool dbool;
  G4ParticleHPDataUsed aFile =
    theNames.GetName(A, Z, M, dirName, tString, dbool);

  G4String filename = aFile.GetName();
  SetAZMs(A, Z, M, aFile);
  if (!dbool || (Z <= 2 && (theBaseZ != Z || theBaseA != A))) {
    hasAnyData = false;
    hasFSData = false;
    hasXsec = false;
    return;
  }
  theData.clear();
  G4ParticleHPManager::GetInstance()->GetDataStream(filename, theData);
  hasFSData = theFinalStatePhotons.InitMean(theData);
  if (hasFSData) {
    targetMass = theFinalStatePhotons.GetTargetMass();
    theFinalStatePhotons.InitAngular(theData);
    theFinalStatePhotons.InitEnergies(theData);
  }
}
