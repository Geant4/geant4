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
// 080718 As for secondary photons, if its mean value has a value of integer,
//        then a sampling of multiplicity that based on Poisson Distribution
//        is not carried out and the mean is used as a multiplicity.
//        modified by T. Koi.
// 080721 Using ClearHistories() methodl for limiting the sum of secondary energies
//        modified by T. Koi.
// 080901 bug fix of too many secnodaries production in nd reactinos by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, May-2023 Basic revision of particle HP classes
//
#include "G4ParticleHPProduct.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicException.hh"
#include "G4ParticleHPContEnergyAngular.hh"
#include "G4ParticleHPDiscreteTwoBody.hh"
#include "G4ParticleHPIsotropic.hh"
#include "G4ParticleHPLabAngularEnergy.hh"
#include "G4ParticleHPNBodyPhaseSpace.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Proton.hh"

G4ParticleHPProduct::G4ParticleHPProduct()
{
  toBeCached val;
  fCache.Put(val);

  if (G4ParticleHPManager::GetInstance()->GetPHCUsePoisson()) {
    theMultiplicityMethod = G4HPMultiPoisson;
  } else {
    theMultiplicityMethod = G4HPMultiBetweenInts;
  }
}

G4ParticleHPProduct::~G4ParticleHPProduct()
{
  delete theDist;
}

void G4ParticleHPProduct::Init(std::istream& aDataFile, const G4ParticleDefinition* projectile)
{
  aDataFile >> theMassCode >> theMass >> theIsomerFlag >> theDistLaw >> theGroundStateQValue
	    >> theActualStateQValue;
  theGroundStateQValue *= CLHEP::eV;
  theActualStateQValue *= CLHEP::eV;
  theYield.Init(aDataFile, CLHEP::eV);
  theYield.Hash();
  if (theDistLaw == 0) {
    // distribution not known, use E-independent, isotropic
    // angular distribution
    theDist = new G4ParticleHPIsotropic;
  }
  else if (theDistLaw == 1) {
    // Continuum energy-angular distribution
    theDist = new G4ParticleHPContEnergyAngular(projectile);
  }
  else if (theDistLaw == 2) {
    // Discrete 2-body scattering
    theDist = new G4ParticleHPDiscreteTwoBody;
  }
  else if (theDistLaw == 3) {
    // Isotropic emission
    theDist = new G4ParticleHPIsotropic;
  }
  else if (theDistLaw == 4) {
    // Discrete 2-body recoil modification not used for now.
    // theDist = new G4ParticleHPDiscreteTwoBody;
    // the above is only temporary;
    // recoils need to be addressed properly
  }
  //    else if(theDistLaw == 5)
  //    {
  // charged particles only, to be used in a later stage. @@@@
  //    }
  else if (theDistLaw == 6) {
    // N-Body phase space
    theDist = new G4ParticleHPNBodyPhaseSpace;
  }
  else if (theDistLaw == 7) {
    // Laboratory angular energy paraetrisation
    theDist = new G4ParticleHPLabAngularEnergy;
  }
  else {
        throw G4HadronicException(__FILE__, __LINE__,
                                  "distribution law unknown to G4ParticleHPProduct");
  }
  if (theDist != nullptr) {
    theDist->SetQValue(theActualStateQValue);
    theDist->Init(aDataFile);
  }
}

G4int G4ParticleHPProduct::GetMultiplicity(G4double anEnergy)
{
  if (theDist == nullptr) {
    fCache.Get().theCurrentMultiplicity = 0;
    return 0;
  }

  G4double mean = theYield.GetY(anEnergy);
  if (mean <= 0.) {
    fCache.Get().theCurrentMultiplicity = 0;
    return 0;
  }
  G4int multi = (theMultiplicityMethod == G4HPMultiPoisson) ?
    (G4int)G4Poisson(mean) : G4lrint(mean);

#ifdef G4VERBOSE
  if (G4ParticleHPManager::GetInstance()->GetDEBUG())
    G4cout << "G4ParticleHPProduct::GetMultiplicity code=" << theMassCode << " M=" << theMass
	   << " multi=" << multi << " mean=" << mean << G4endl;
#endif
  fCache.Get().theCurrentMultiplicity = multi;

  return multi;
}

G4ReactionProductVector* G4ParticleHPProduct::Sample(G4double anEnergy, G4int multi)
{
  if (theDist == nullptr) {
    return nullptr;
  }
  auto result = new G4ReactionProductVector;

  theDist->SetTarget(fCache.Get().theTarget);
  theDist->SetProjectileRP(fCache.Get().theProjectileRP);
  G4ReactionProduct* tmp;
  theDist->ClearHistories();

  for (G4int i = 0; i < multi; ++i) {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    if (tmp != nullptr) {
      result->push_back(tmp);
#ifdef G4VERBOSE
      if (G4ParticleHPManager::GetInstance()->GetDEBUG())
	G4cout << "multi=" << multi << " i=" << i << " G4ParticleHPProduct::Sample "
	       << tmp->GetDefinition()->GetParticleName() << " E=" << tmp->GetKineticEnergy()
	       << G4endl;
#endif
    }
  }
  if (multi == 0) {
    tmp = theDist->Sample(anEnergy, theMassCode, theMass);
    delete tmp;
  }
  return result;
}
