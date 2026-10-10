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
// J. Shen, Oct-2026 - Correlated sampling of inclusive final state data with residual-nucleus
// yields, preserving energy, charge, and baryon number.
//

#include "G4ParticleHPEnAngCorrelation.hh"

#include "G4Fragment.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4LorentzRotation.hh"
#include "G4LorentzVector.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhotonEvaporation.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

G4ParticleHPEnAngCorrelation::G4ParticleHPEnAngCorrelation(const G4ParticleDefinition* proj)
{
  theProjectile = (nullptr == proj) ? G4Neutron::Neutron() : proj;
  toBeCached val;
  fCache.Put(val);
}

G4ParticleHPEnAngCorrelation::~G4ParticleHPEnAngCorrelation()
{
  delete[] theProducts;
}

void G4ParticleHPEnAngCorrelation::Init(std::istream& aDataFile)
{
  inCharge = true;
  aDataFile >> targetMass >> frameFlag >> nProducts;
  const std::size_t psize = nProducts > 0 ? nProducts : 1;
  theProducts = new G4ParticleHPProduct[psize];
  for (std::size_t i = 0; i < psize; ++i)
  {
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
      G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of "
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
      for (auto const& ptr : *temp)
      {
        delete ptr;
      }
      throw G4HadronicException(__FILE__, __LINE__, "SampleOne: Yield not correct");
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

    G4LorentzVector Ptmp(aIncidentPart.GetMomentum(), aIncidentPart.GetTotalEnergy());

    toZ.rotateZ(-1 * Ptmp.phi());
    toZ.rotateY(-1 * Ptmp.theta());
  }
  fCache.Get().theTotalMeanEnergy = 0;
  G4LorentzRotation toLab(toZ.inverse());
  // toLab only change axis NOT to LAB system

  if (isInclusive && correlatedSampling)
  {
    G4ReactionProductVector* correlatedResult = SampleCorrelated(anEnergy, theCMS, toLab);
    if (correlatedResult != nullptr)
    {
      delete result;
      return correlatedResult;
    }
  }

  for (i = 0; i < nProducts; ++i)
  {
    G4int nPart = theProducts[i].GetMultiplicity(anEnergy);
    it = theProducts[i].Sample(anEnergy, nPart);
    G4double aMeanEnergy = theProducts[i].MeanEnergyOfThisInteraction();

    // 151120 TK Modified for solving reproducibility problem
    // This change may have side effect.
    if (aMeanEnergy >= 0)
    {
      fCache.Get().theTotalMeanEnergy += aMeanEnergy;
    }
    else
    {
      fCache.Get().theTotalMeanEnergy = anEnergy / nProducts + theProducts[i].GetQValue();
    }
    if (it != nullptr)
    {
      for (unsigned int ii = 0; ii < it->size(); ++ii)
      {
        ToLab(it->operator[](ii), i, theCMS, toLab);
        result->push_back(it->operator[](ii));
      }
      delete it;
    }
  }
  return result;
}

void G4ParticleHPEnAngCorrelation::ToLab(G4ReactionProduct* aProduct, G4int iProduct,
                                         const G4ReactionProduct& theCMS,
                                         const G4LorentzRotation& toLab)
{
  G4LorentzVector pTmp1(aProduct->GetMomentum(), aProduct->GetTotalEnergy());
  pTmp1 = toLab * pTmp1;
  aProduct->SetMomentum(pTmp1.vect());
  aProduct->SetTotalEnergy(pTmp1.e());

  if (frameFlag == 1)  // target rest //TK 100413 should be LAB?
  {
    aProduct->Lorentz(*aProduct,
                      -1. * (*fCache.Get().theTarget));  // TK 100413 Is this really need?
  }
  else if (frameFlag == 2)  // CMS
  {
#ifdef G4PHPDEBUG
    if (G4ParticleHPManager::GetInstance()->GetDEBUG())
      G4cout << "G4ParticleHPEnAngCorrelation: before Lorentz boost "
             << aProduct->GetKineticEnergy() << " " << aProduct->GetMomentum() << G4endl;
#endif
    aProduct->Lorentz(*aProduct, -1. * theCMS);
#ifdef G4PHPDEBUG
    if (G4ParticleHPManager::GetInstance()->GetDEBUG())
      G4cout << "G4ParticleHPEnAngCorrelation: after Lorentz boost " << aProduct->GetKineticEnergy()
             << " " << aProduct->GetMomentum() << G4endl;
#endif
  }
  // TK120515 migrate frameFlag (MF6 LCT) = 3
  else if (frameFlag == 3)  // CMS A<=4 other LAB
  {
    if (theProducts[iProduct].GetMassCode() > 2004.5)  // Alpha AWP 3.96713
    {
      // LAB
      aProduct->Lorentz(*aProduct,
                        -1. * (*fCache.Get().theTarget));  // TK 100413 Is this really need?
    }
    else
    {
      // CMS
      aProduct->Lorentz(*aProduct, -1. * theCMS);
    }
#ifdef G4PHPDEBUG
    if (G4ParticleHPManager::GetInstance()->GetDEBUG())
      G4cout << "G4ParticleHPEnAngCorrelation: after Lorentz boost " << aProduct->GetKineticEnergy()
             << " " << aProduct->GetMomentum() << G4endl;
#endif
  }
  else
  {
    throw G4HadronicException(__FILE__, __LINE__,
                              "Sample: The frame of the final state is not specified");
  }
}

G4ReactionProductVector*
G4ParticleHPEnAngCorrelation::SampleCorrelated(G4double anEnergy, const G4ReactionProduct& theCMS,
                                               const G4LorentzRotation& toLab)
{
  const G4ReactionProduct* target = fCache.Get().theTarget;
  const G4ReactionProduct* projectile = fCache.Get().theProjectileRP;
  if (target == nullptr || projectile == nullptr || target->GetDefinition() == nullptr)
  {
    return nullptr;
  }
  const G4int targZ = G4lrint(target->GetDefinition()->GetPDGCharge() / CLHEP::eplus);
  const G4int targA = target->GetDefinition()->GetBaryonNumber();
  const G4int projZ = G4lrint(theProjectile->GetPDGCharge() / CLHEP::eplus);
  const G4int projA = theProjectile->GetBaryonNumber();
  if (targA <= 0)
  {
    return nullptr;
  }

  // classify the products of the data by their mass code ZA
  struct Light
  {
      G4int index, Z, A;
      G4double yield;
  };
  struct Residual
  {
      G4int index;
      G4double yield;
  };
  std::vector<Light> lights;
  std::vector<Residual> residuals;
  G4double sumResidual = 0.0;
  for (G4int i = 0; i < nProducts; ++i)
  {
    const G4int za = G4lrint(theProducts[i].GetMassCode());
    if (za <= 0)
    {
      continue;
    }  // photons: replaced by the de-excitation of the residual
    const G4int Z = za / 1000;
    const G4int A = za % 1000;
    const G4double y = theProducts[i].GetMeanYield(anEnergy);
    if (y <= 0.0)
    {
      continue;
    }
    if ((za == 1 ||  // n
         za == 1001 ||  // p
         za == 1002 ||  // d
         za == 1003 ||  // t
         za == 2003 ||  // he3
         za == 2004  // alpha
         )
        && !(Z == projZ + targZ && A == projA + targA))
    {
      lights.push_back({i, Z, A, y});
    }
    else if (A > 4 || Z == projZ + targZ)
    {
      residuals.push_back({i, y});
      sumResidual += y;
    }
  }
  // the residual yields must be the branching ratios; otherwise these data are not of that kind
  if (residuals.empty() || sumResidual < 0.9 || sumResidual > 1.1)
  {
    return nullptr;
  }

  // initial state in the lab, with Geant4 nuclear masses (excluding shell electrons)
  const G4double targMass = G4NucleiProperties::GetNuclearMass(targA, targZ);
  const G4LorentzVector initial =
    G4LorentzVector(projectile->GetMomentum(), projectile->GetTotalEnergy())
    + G4LorentzVector(target->GetMomentum(),
                      std::sqrt(targMass * targMass + target->GetMomentum().mag2()));
  const G4double sqrts = initial.m();

  // exit channels: a residual and the light particles that lead to it. Several combinations can
  // give the same residual (e.g. alpha or 2n+2p); kinematically closed ones are dropped and the
  // open ones share the residual yield in proportion to the product of their mean yields (below
  // the multi-particle thresholds there is one open combination per residual).
  struct Channel
  {
      G4int resZ, resA;
      std::vector<G4int> counts;
      G4double weight;
  };
  std::vector<Channel> channels;
  for (auto const& [resIndex, resYield] : residuals)
  {
    const G4int resZA = G4lrint(theProducts[resIndex].GetMassCode());
    const G4int resZ = resZA / 1000;
    const G4int resA = resZA % 1000;
    const G4int dZ = projZ + targZ - resZ;
    const G4int dA = projA + targA - resA;
    if (dZ < 0 || dA < 0)
    {
      continue;
    }
    const G4double resMass = G4NucleiProperties::GetNuclearMass(resA, resZ);
    std::vector<std::vector<G4int>> combos;
    std::vector<G4double> weights;
    std::vector<G4int> counts(lights.size(), 0);
    // recursively enumerate all combinations of light particles that give the right residual
    std::function<void(std::size_t, G4int, G4int, G4int)> enumerate = [&](std::size_t k, G4int z,
                                                                          G4int a, G4int n) {
      if (z == 0 && a == 0)
      {
        // valid combination found; compute its weight and check if it is kinematically open
        G4double w = 1.0;
        G4double mass = resMass;
        for (std::size_t j = 0; j < lights.size(); ++j)
        {
          w *= std::pow(lights[j].yield, counts[j]);
          mass += counts[j] * G4NucleiProperties::GetNuclearMass(lights[j].A, lights[j].Z);
        }
        if (mass < sqrts)
        {
          combos.push_back(counts);
          weights.push_back(w);
        }
        return;
      }
      // no particle types left to try, or already 8 particles
      if (k == lights.size() || n >= 8)
      {
        return;
      }
      // Also search to the next light particle
      enumerate(k + 1, z, a, n);
      if (lights[k].Z <= z && lights[k].A <= a)
      {
        // if we can use this light particle, increment its count and continue the search
        ++counts[k];
        enumerate(k, z - lights[k].Z, a - lights[k].A, n + 1);
        --counts[k];
      }
    };  // end of enumerate

    enumerate(0, dZ, dA, 0);
    G4double wsum = 0.0;
    for (auto w : weights)
    {
      wsum += w;
    }
    for (std::size_t c = 0; c < combos.size(); ++c)
    {
      channels.push_back({resZ, resA, combos[c], resYield * weights[c] / wsum});
    }
  }

  // choose a channel and sample its light particles from their (channel-summed) energy-angle
  // distributions; draws that leave the residual below its ground state are rejected. A channel
  // that cannot be sampled (e.g. open only through the linear interpolation of the yields) is
  // removed and the choice is repeated among the others.
  auto result = new G4ReactionProductVector;
  G4LorentzVector residual;
  G4int resZ = 0;
  G4int resA = 0;
  G4bool accepted = false;
  constexpr G4int channel_choice_max = 10;
  constexpr G4int attempt_max = 1024;
  for (G4int choice = 0; choice < channel_choice_max && !accepted; ++choice)
  {
    G4double wsum =
      std::accumulate(channels.begin(), channels.end(), 0.0,
                      [](G4double sum, const Channel& ch) { return sum + ch.weight; });
    if (wsum <= 0.0)
    {
      break;
    }
    std::size_t iCh = channels.size() - 1;
    G4double rc = G4UniformRand() * wsum;
    for (std::size_t c = 0; c < channels.size(); ++c)
    {
      rc -= channels[c].weight;
      if (rc <= 0.0 && channels[c].weight > 0.0)
      {
        iCh = c;
        break;
      }
    }
    const Channel& ch = channels[iCh];
    resZ = ch.resZ;
    resA = ch.resA;
    const G4double resMass = G4NucleiProperties::GetNuclearMass(resA, resZ);
    for (G4int attempt = 0; attempt < attempt_max && !accepted; ++attempt)
    {
      for (auto p : *result)
      {
        delete p;
      }
      result->clear();
      G4LorentzVector emitted(0., 0., 0., 0.);
      G4bool complete = true;
      for (std::size_t j = 0; j < lights.size() && complete; ++j)
      {
        const G4int n = ch.counts[j];
        if (n == 0)
        {
          continue;
        }
        G4ReactionProductVector* it = theProducts[lights[j].index].Sample(anEnergy, n);
        if (it == nullptr || (G4int)it->size() != n)
        {
          complete = false;
        }
        if (it != nullptr)
        {
          for (auto p : *it)
          {
            ToLab(p, lights[j].index, theCMS, toLab);
            emitted += G4LorentzVector(p->GetMomentum(), p->GetTotalEnergy());
            result->push_back(p);
          }
          delete it;
        }
      }
      if (!complete)
      {
        continue;
      }
      residual = initial - emitted;
      accepted = residual.e() > 0.0 && residual.m2() > 0.0 && residual.m() - resMass > -CLHEP::keV;
    }
    if (!accepted)
    {
      channels[iCh].weight = 0.0;
    }
  }
  if (!accepted)
  {
    for (auto p : *result)
    {
      delete p;
    }
    delete result;
    return nullptr;
  }
  const G4double resMass = G4NucleiProperties::GetNuclearMass(resA, resZ);

  // de-excite the residual using photonEvaporation
  if (residual.m() < resMass)
  {
    residual.setE(std::sqrt(resMass * resMass + residual.vect().mag2()));
  }
  G4Fragment nucleus(resA, resZ, residual);
  static G4ThreadLocal G4PhotonEvaporation* photonEvaporation = nullptr;
  if (photonEvaporation == nullptr)
  {
    photonEvaporation = new G4PhotonEvaporation();
    photonEvaporation->Initialise();
  }
  if (nucleus.GetExcitationEnergy() > CLHEP::keV)
  {
    G4FragmentVector emittedFragments;
    photonEvaporation->BreakUpChain(&emittedFragments, &nucleus);
    for (auto frag : emittedFragments)
    {
      auto p = new G4ReactionProduct;
      const G4ParticleDefinition* def = frag->GetParticleDefinition();
      if (def == nullptr)
      {
        def = G4Gamma::Gamma();
      }
      p->SetDefinition(def);
      p->SetMomentum(frag->GetMomentum().vect());
      p->SetTotalEnergy(frag->GetMomentum().e());
      result->push_back(p);
      delete frag;
    }
  }
  auto res = new G4ReactionProduct;
  const G4double exRes = nucleus.GetExcitationEnergy();
  res->SetDefinition(exRes > CLHEP::keV ? G4IonTable::GetIonTable()->GetIon(
                       resZ, resA, exRes, G4Ions::FloatLevelBase(nucleus.GetFloatingLevelNumber()))
                                        : G4IonTable::GetIonTable()->GetIon(resZ, resA));
  res->SetMomentum(nucleus.GetMomentum().vect());
  res->SetTotalEnergy(nucleus.GetMomentum().e());
  result->push_back(res);

  // the residual excitation is emitted here: no photons for the binding energy downstream
  fCache.Get().theTotalMeanEnergy = DBL_MAX;
  return result;
}
