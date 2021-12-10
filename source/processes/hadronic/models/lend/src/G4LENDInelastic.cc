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
//  File:   G4LENDInelastic.cc                                                //
//  Date:   24 March 2020                                                     //
//  Author: Dennis Wright                                                     //
//                                                                            //
//  Description: model for inelastic scattering of neutrons, light ions and   //
//               gammas at energies of order 20 MeV and lower.                //
//               This model uses GIDI particle data which are stored mostly   //
//               in spectrum mode.  In this mode, spectra are reproduced      //
//               for each possible particle type which can result from a      //
//               given interaction.  Unlike Geant4, this is done without      //
//               consideration of event-by-event conservation rules.          //
//               Indeed, forcing such conservation on GIDI output products    //
//               introduces correlations and distortions in the resulting     //
//               spectra which are not present in the data.                   //
//                                                                            //
//               In order to use GIDI data within the Geant4 framework, a     //
//               minimal event-by-event baryon number conservation is         //
//               enforced which allows deviations of up to 1 GeV without      //
//               giving warnings.  Neither charge, nor energy, nor momentum   //
//               conservation is enforced.  Under this scheme, light          //
//               fragment (n, p, d, t, alpha) spectra are well reproduced     //
//               after a large number of events.  Charge and energy           //
//               conservation also approach their event-by-event values in    //
//               this limit.  The mass, charge and energy distributions of    //
//               large fragments, however, are not expected to reproduce the  //
//               data very well.  This is a result of forcing the crude       //
//               baryon number conservation and ensuring that the light       //
//               fragment spectra are correct.                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4LENDInelastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include <algorithm>
#include <random>
  
G4HadFinalState* G4LENDInelastic::ApplyYourself(const G4HadProjectile& aTrack,
                                                G4Nucleus& aTarg)
{
  G4ThreeVector projMom = aTrack.Get4Momentum().vect();
  G4double temp = aTrack.GetMaterial()->GetTemperature();

  G4int iZ = aTarg.GetZ_asInt();
  G4int iA = aTarg.GetA_asInt();
  G4int iM = 0;
  if (aTarg.GetIsotope() != nullptr) iM = aTarg.GetIsotope()->Getm();

  G4double ke = aTrack.GetKineticEnergy();

  G4HadFinalState* theResult = &theParticleChange;
  theResult->Clear();

  G4GIDI_target* aGIDITarget =
    get_target_from_map(lend_manager->GetNucleusEncoding(iZ, iA, iM) );
  if (aGIDITarget == nullptr) {
//    G4cout << " No target found " << G4endl; 
    theParticleChange.Clear();
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }
// return returnUnchanged(aTrack, theResult);

  // Get GIDI final state products for givern target and projectile
  G4int loop(0);
  G4int loopMax = 1000;
  std::vector<G4GIDI_Product>* products;
  do {
    products = aGIDITarget->getOthersFinalState(ke*MeV, temp, MyRNG, NULL);
    loop++;
  } while (products == nullptr && loop < loopMax);

  // G4LENDInelastic accepts all light fragments and gammas from GIDI (A < 5)
  // and removes any heavy fragments which cause large baryon number violation.
  // Charge and energy non-conservation still occur, but over a large number 
  // of events, this improves on average.

  if (loop > loopMax - 1) {
//    G4cout << " too many loops, return initial state " << G4endl;

    theParticleChange.Clear();
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;

//    if (aTrack.GetDefinition() == G4Proton::Proton() ||
//        aTrack.GetDefinition() == G4Neutron::Neutron() ) {
//       theResult = preco->ApplyYourself(aTrack, aTarg);
//     } else {
//       theResult = returnUnchanged(aTrack, theResult);
//     }

  } else {
    G4int iTotZ = iZ + aTrack.GetDefinition()->GetAtomicNumber();
    G4int iTotA = iA + aTrack.GetDefinition()->GetAtomicMass();

    // Loop over GIDI products and separate light from heavy fragments
    G4int GZtot(0);
    G4int GAtot(0);
    G4int productA(0);
    G4int productZ(0);
    std::vector<G4int> lightProductIndex;
    std::vector<G4int> heavyProductIndex;
    for (G4int i = 0; i < int( products->size() ); i++ ) {
      productA = (*products)[i].A;
      if (productA < 5) { 
        lightProductIndex.push_back(i);
        GZtot += (*products)[i].Z;
        GAtot += productA;
      } else {
        heavyProductIndex.push_back(i);
      }
    } 

    // Randomize order of heavies to correct somewhat for sampling bias
    // std::random_shuffle(heavyProductIndex.begin(), heavyProductIndex.end() );
    // std::cout << " Heavy product index before shuffle : " ;
    // for (G4int i = 0; i < int(heavyProductIndex.size() ); i++) std::cout << heavyProductIndex[i] << ", " ;
    // std::cout << std::endl;

    auto rng = std::default_random_engine {};
    std::shuffle(heavyProductIndex.begin(), heavyProductIndex.end(), rng);

    // std::cout << " Heavy product index after shuffle : " ;
    // for (G4int i = 0; i < int(heavyProductIndex.size() ); i++) std::cout << heavyProductIndex[i] << ", " ;
    // std::cout << std::endl;

    std::vector<G4int> savedHeavyIndex;    
    G4int itest(0);
    for (G4int i = 0; i < int(heavyProductIndex.size() ); i++) {
      itest = heavyProductIndex[i];
      productA = (*products)[itest].A;
      productZ = (*products)[itest].Z;
      if ((GAtot + productA <= iTotA) && (GZtot + productZ <= iTotZ) ) {
        savedHeavyIndex.push_back(itest);
        GZtot += productZ;
        GAtot += productA;
      }
    }

/*
    G4cout << " saved light products = ";
    for (G4int k = 0; k < int(lightProductIndex.size() ); k++ ) {
      itest = lightProductIndex[k];
      G4cout << "(" << (*products)[itest].Z << ", " << (*products)[itest].A << "),  ";
    }
    G4cout << G4endl;

    G4cout << " saved heavy products = ";
    for (G4int k = 0; k < int(savedHeavyIndex.size() ); k++ ) {
      itest = savedHeavyIndex[k];
      G4cout << "(" << (*products)[itest].Z << ", " << (*products)[itest].A << "),  ";
    }
    G4cout << G4endl;
*/
    // Now convert saved products to Geant4 particles
    // Note that, at least for heavy fragments, GIDI masses and Geant4 masses
    // have slightly different values.

    G4DynamicParticle* theSec = nullptr;
    G4ThreeVector Psum;
    for (G4int i = 0; i < int(lightProductIndex.size() ); i++) {
      itest = lightProductIndex[i];
      productZ = (*products)[itest].Z;
      productA = (*products)[itest].A;
      theSec = new G4DynamicParticle();
      if (productA == 1 && productZ == 0) {
        theSec->SetDefinition(G4Neutron::Neutron() );
      } else if (productA == 1 && productZ == 1) {
        theSec->SetDefinition(G4Proton::Proton() );
      } else if (productA == 2 && productZ == 1) {
        theSec->SetDefinition(G4Deuteron::Deuteron() );
      } else if (productA == 3 && productZ == 1) {
        theSec->SetDefinition(G4Triton::Triton() );
      } else if (productA == 4 && productZ == 2) {
        theSec->SetDefinition(G4Alpha::Alpha() );
      } else {
        theSec->SetDefinition(G4Gamma::Gamma() );
      }

      G4ThreeVector momentum((*products)[itest].px*MeV, 
                             (*products)[itest].py*MeV,
                             (*products)[itest].pz*MeV );
      Psum += momentum;
      theSec->SetMomentum(momentum);
//      theResult->AddSecondary(theSec);
      theParticleChange.AddSecondary(theSec, secID);
    }

    G4int productM(0);
    for (G4int i = 0; i < int(savedHeavyIndex.size() ); i++) {
      itest = savedHeavyIndex[i];
      productZ = (*products)[itest].Z;
      productA = (*products)[itest].A;
      productM = (*products)[itest].m;
      theSec = new G4DynamicParticle();
      theSec->SetDefinition(G4IonTable::GetIonTable()->GetIon(productZ,
                                                              productA,
                                                              productM) );
      G4ThreeVector momentum((*products)[itest].px*MeV,
                             (*products)[itest].py*MeV,
                             (*products)[itest].pz*MeV );
      Psum += momentum;
      theSec->SetMomentum(momentum);
//      theResult->AddSecondary(theSec);
      theParticleChange.AddSecondary(theSec, secID);
    }

    // Create heavy fragment if necessary to try to balance A, Z
    // Note: this step is only required to prevent warnings at the process level
    //       where "catastrophic" non-conservation tolerances are set to 1 GeV.
    //       The residual generated will not necessarily be the one that would 
    //       occur in the actual reaction.
    if (iTotA - GAtot > 1) {
      theSec = new G4DynamicParticle();
      if (iTotZ == GZtot) {
        // Special case when a nucleus of only neutrons is requested
        // Violate charge conservation and set Z = 1
        // G4cout << " Z = 1, A = "<< iTotA - GAtot << " created " << G4endl;
        theSec->SetDefinition(G4IonTable::GetIonTable()->GetIon(1, iTotA-GAtot, 0) );
      } else {
        theSec->SetDefinition(G4IonTable::GetIonTable()->GetIon(iTotZ-GZtot, iTotA-GAtot, 0) );
      }
      theSec->SetMomentum(projMom - Psum);
//      theResult->AddSecondary(theSec);
      theParticleChange.AddSecondary(theSec, secID);
    }
  } // loop OK

  delete products;
//  theResult->SetStatusChange( stopAndKill );
  theParticleChange.SetStatusChange( stopAndKill );
//  return theResult; 
  return &theParticleChange;
}
