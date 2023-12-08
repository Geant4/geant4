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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 100413 Fix bug in incidence energy by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
// June-2019 - E. Mendoza --> Part of the code trying to preserve the baryonic number has been
// deleted. One has to assume that it is not preserved when using ENDF-6 data and it caused
// problems.
//
// V. Ivanchenko, May-2023 Basic revision of particle HP classes
//

#include "G4ParticleHPEnAngCorrelation.hh"

#include "G4IonTable.hh"
#include "G4LorentzRotation.hh"
#include "G4LorentzVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ParticleHPManager.hh"
#include "G4Neutron.hh"

G4ParticleHPEnAngCorrelation::G4ParticleHPEnAngCorrelation(const G4ParticleDefinition* proj)
{
  theProjectile = (nullptr == proj) ? G4Neutron::Neutron() : proj;
  toBeCached val;
  fCache.Put( val );
}

G4ParticleHPEnAngCorrelation::~G4ParticleHPEnAngCorrelation()
{
  delete[] theProducts;
}

void G4ParticleHPEnAngCorrelation::Init(std::istream& aDataFile)
{
  inCharge = true;
  aDataFile >> targetMass >> frameFlag >> nProducts;
  //G4cout << "G4ParticleHPEnAngCorrelation::Init " << theProjectile->GetParticleName()
  //	 << " frameFlag=" << frameFlag << " N=" << nProducts << " Mass=" << targetMass << G4endl;
  theProducts = new G4ParticleHPProduct[nProducts];
  for (G4int i = 0; i < nProducts; ++i) {
    theProducts[i].Init(aDataFile, theProjectile);
  }
}

G4ReactionProduct* G4ParticleHPEnAngCorrelation::SampleOne(G4double anEnergy)
{
  auto result = new G4ReactionProduct;

  // do we have an appropriate distribution
  if (nProducts != 1)
    throw G4HadronicException(__FILE__, __LINE__, "More than one product in SampleOne");

  // get the result
  G4ReactionProductVector* temp = nullptr;
  G4int i = 0;

  G4int icounter = 0;
  G4int icounter_max = 1024;
  while (temp == nullptr)
  {
    ++icounter;
    if (icounter > icounter_max)
    {
      G4cout << "Loop-counter exceeded the threshold value at " 
             << __LINE__ << "th line of "
             << __FILE__ << "." << G4endl;
      break;
    }
    temp = theProducts[++i].Sample(anEnergy, 1);
  }
  if (nullptr != temp)
  {
    if (temp->size() == 1)
    {
      result = (*temp)[0];
    }
    else
    {
      for (auto const& ptr : *temp) { delete ptr; }
      throw G4HadronicException(__FILE__, __LINE__,
				"SampleOne: Yield not correct");
    }
  }
  delete temp;
  return result;
}

G4ReactionProductVector* G4ParticleHPEnAngCorrelation::Sample(G4double anEnergy)
{
  auto result = new G4ReactionProductVector;
  G4int i;
  G4ReactionProductVector* it;
  G4ReactionProduct theCMS;
  G4LorentzRotation toZ;
  /*
  G4cout << "G4ParticleHPEnAngCorrelation::Sample for E=" << anEnergy 
	 << " flag=" << frameFlag << " nProducts=" << nProducts <<G4endl;
  */
  if (frameFlag == 2 || frameFlag == 3)  // Added for particle HP
  {
    G4ThreeVector the3IncidentPart = fCache.Get().theProjectileRP->GetMomentum();
      // theProjectileRP has value in LAB
    G4double nEnergy = fCache.Get().theProjectileRP->GetTotalEnergy();
    G4ThreeVector the3Target = fCache.Get().theTarget->GetMomentum();
      // theTarget has value in LAB
    G4double tEnergy = fCache.Get().theTarget->GetTotalEnergy();
    G4double totE = nEnergy + tEnergy;
    G4ThreeVector the3CMS = the3Target + the3IncidentPart;
    theCMS.SetMomentum(the3CMS);
    G4double cmsMom = std::sqrt(the3CMS * the3CMS);
    G4double sqrts = std::sqrt((totE - cmsMom) * (totE + cmsMom));
    theCMS.SetMass(sqrts);
    theCMS.SetTotalEnergy(totE);
    G4ReactionProduct aIncidentPart;
    aIncidentPart.Lorentz(*fCache.Get().theProjectileRP, theCMS);
    // TKDB 100413
    // ENDF-6 Formats Manual ENDF-102
    // CHAPTER 6. FILE 6: PRODUCT ENERGY-ANGLE DISTRIBUTIONS
    // LCT Reference system for secondary energy and angle
    // incident energy is always given in the LAB system
    anEnergy = fCache.Get().theProjectileRP->GetKineticEnergy();
    // should be same argumment of "anEnergy"

    G4LorentzVector Ptmp(aIncidentPart.GetMomentum(),
                         aIncidentPart.GetTotalEnergy());

    toZ.rotateZ(-1 * Ptmp.phi());
    toZ.rotateY(-1 * Ptmp.theta());
  }
  fCache.Get().theTotalMeanEnergy = 0;
  G4LorentzRotation toLab(toZ.inverse());
  // toLab only change axis NOT to LAB system

  for (i = 0; i < nProducts; ++i) {
    G4int nPart = theProducts[i].GetMultiplicity(anEnergy);
    it = theProducts[i].Sample(anEnergy, nPart);
    G4double aMeanEnergy = theProducts[i].MeanEnergyOfThisInteraction();

    // 151120 TK Modified for solving reproducibility problem
    // This change may have side effect.
    if (aMeanEnergy >= 0) {
      fCache.Get().theTotalMeanEnergy += aMeanEnergy;
    }
    else {
      fCache.Get().theTotalMeanEnergy = anEnergy/nProducts + theProducts[i].GetQValue();
    }
    if (it != nullptr) {
      for (unsigned int ii = 0; ii < it->size(); ++ii) {
        G4LorentzVector pTmp1(it->operator[](ii)->GetMomentum(),
                              it->operator[](ii)->GetTotalEnergy());
        pTmp1 = toLab * pTmp1;
        it->operator[](ii)->SetMomentum(pTmp1.vect());
        it->operator[](ii)->SetTotalEnergy(pTmp1.e());

        if (frameFlag == 1)  // target rest //TK 100413 should be LAB?
        {
          it->operator[](ii)->Lorentz(
            *(it->operator[](ii)),
            -1. * (*fCache.Get().theTarget));  // TK 100413 Is this really need?
        }
        else if (frameFlag == 2)  // CMS
        {
#ifdef G4PHPDEBUG
          if (G4ParticleHPManager::GetInstance()->GetDEBUG())
            G4cout << "G4ParticleHPEnAngCorrelation: before Lorentz boost "
                   << it->at(ii)->GetKineticEnergy() << " "
                   << it->at(ii)->GetMomentum() << G4endl;
#endif
          it->operator[](ii)->Lorentz(*(it->operator[](ii)), -1. * theCMS);
#ifdef G4PHPDEBUG
          if (G4ParticleHPManager::GetInstance()->GetDEBUG())
            G4cout << "G4ParticleHPEnAngCorrelation: after Lorentz boost "
                   << it->at(ii)->GetKineticEnergy() << " "
                   << it->at(ii)->GetMomentum() << G4endl;
#endif
        }
        // TK120515 migrate frameFlag (MF6 LCT) = 3
        else if (frameFlag == 3)  // CMS A<=4 other LAB
        {
          if (theProducts[i].GetMassCode() > 2004.5)  // Alpha AWP 3.96713
          {
            // LAB
            it->operator[](ii)->Lorentz(
              *(it->operator[](ii)),
              -1. * (*fCache.Get().theTarget));  // TK 100413 Is this really need?
#ifdef G4PHPDEBUG
            if (G4ParticleHPManager::GetInstance()->GetDEBUG())
              G4cout << "G4ParticleHPEnAngCorrelation: after Lorentz boost "
                     << it->at(ii)->GetKineticEnergy() << " "
                     << it->at(ii)->GetMomentum() << G4endl;
#endif
          }
          else {
            // CMS
            it->operator[](ii)->Lorentz(*(it->operator[](ii)), -1. * theCMS);
#ifdef G4PHPDEBUG
            if (G4ParticleHPManager::GetInstance()->GetDEBUG())
              G4cout << "G4ParticleHPEnAngCorrelation: after Lorentz boost "
                     << it->at(ii)->GetKineticEnergy() << " "
                     << it->at(ii)->GetMomentum() << G4endl;
#endif
          }
        }
        else
        {
          throw G4HadronicException(__FILE__, __LINE__,
            "Sample: The frame of the final state is not specified");
        }
        result->push_back(it->operator[](ii));
      }
      delete it;
    }
  }
  return result;
}
