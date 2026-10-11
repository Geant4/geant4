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
// 070523 add neglecting doppler broadening on the fly. T. Koi
// 070613 fix memory leaking by T. Koi
// 071002 enable cross section dump by T. Koi
// 080428 change checking point of "neglecting doppler broadening" flag
//        from GetCrossSection to BuildPhysicsTable by T. Koi
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko July-2023 converted back
//
#include "G4NeutronHPCaptureData.hh"

#include "G4AutoLock.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4HadronicParameters.hh"
#include "G4Material.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleHPData.hh"
#include "G4ParticleHPIsoData.hh"
#include "G4ParticleHPVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4Isotope.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsTable.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"

G4bool G4NeutronHPCaptureData::fLock = true;
G4PhysicsTable* G4NeutronHPCaptureData::theCrossSections = nullptr;
std::vector<std::vector<G4PhysicsVector*>>* G4NeutronHPCaptureData::theIsoCrossSections = nullptr;

namespace
{
G4Mutex theHPCaptureData = G4MUTEX_INITIALIZER;
}

G4NeutronHPCaptureData::G4NeutronHPCaptureData() : G4VCrossSectionDataSet("NeutronHPCaptureXS")
{
  emax = 20. * CLHEP::MeV;
  fManager = G4ParticleHPManager::GetInstance();
}

G4NeutronHPCaptureData::~G4NeutronHPCaptureData()
{
  if (isFirst)
  {
    if (nullptr != theCrossSections) theCrossSections->clearAndDestroy();
    delete theCrossSections;
    theCrossSections = nullptr;
    ClearIsotopeTables();
  }
}

G4bool G4NeutronHPCaptureData::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int,
                                               const G4Element*, const G4Material*)
{
  return true;
}

G4double G4NeutronHPCaptureData::GetIsoCrossSection(const G4DynamicParticle* dp, G4int Z,
                                                    G4int A, const G4Isotope* iso,
                                                    const G4Element* element,
                                                    const G4Material* material)
{
  // The cache key includes the isotope: G4CrossSectionDataStore asks for every
  // isotope of an element at the same energy to sample the target isotope.
  if (dp->GetKineticEnergy() == ke_cache && element == element_cache && material == material_cache
      && iso == isotope_cache)
    return xs_cache;

  ke_cache = dp->GetKineticEnergy();
  element_cache = element;
  material_cache = material;
  isotope_cache = iso;

  // Isotope-wise cross section; fall back to the element cross section if the
  // isotope has no data of its own
  const G4PhysicsVector* isoVec = nullptr;
  if (nullptr != iso && nullptr != element && nullptr != theIsoCrossSections)
  {
    std::size_t idx = element->GetIndex();
    if (idx < theIsoCrossSections->size())
    {
      const auto& vecs = (*theIsoCrossSections)[idx];
      for (std::size_t j = 0; j < element->GetNumberOfIsotopes() && j < vecs.size(); ++j)
      {
        if (element->GetIsotope((G4int)j) == iso)
        {
          isoVec = vecs[j];
          break;
        }
      }
    }
  }

  G4double xs;
  if (nullptr != isoVec && dp->GetKineticEnergy() < emax)
  {
    G4double isoMass = G4NucleiProperties::GetNuclearMass(A, Z) / CLHEP::neutron_mass_c2;
    xs = DopplerAverage(dp, isoVec, isoMass, material->GetTemperature());
  }
  else
  {
    xs = GetCrossSection(dp, element, material->GetTemperature());
  }
  xs_cache = xs;
  return xs;
}

void G4NeutronHPCaptureData::ClearIsotopeTables()
{
  if (nullptr == theIsoCrossSections) return;
  for (auto& vecs : *theIsoCrossSections)
  {
    for (auto v : vecs)
    {
      delete v;
    }
  }
  delete theIsoCrossSections;
  theIsoCrossSections = nullptr;
}

void G4NeutronHPCaptureData::BuildIsotopeTables()
{
  // Load the capture cross section of every isotope of every element on its own
  // (abundance 100 %), so that the target isotope can be sampled correctly.
  ClearIsotopeTables();
  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  theIsoCrossSections = new std::vector<std::vector<G4PhysicsVector*>>(numberOfElements);

  G4String dirName = fManager->GetNeutronHPPath() + "/Capture";
  auto theElementTable = G4Element::GetElementTable();
  for (std::size_t i = 0; i < numberOfElements; ++i)
  {
    auto elm = (*theElementTable)[i];
    std::size_t nIso = elm->GetNumberOfIsotopes();
    auto& vecs = (*theIsoCrossSections)[i];
    vecs.assign(nIso, nullptr);
    if (nIso < 2) continue;  // single isotope: element cross section is the isotope one

    for (std::size_t j = 0; j < nIso; ++j)
    {
      const G4Isotope* iso = elm->GetIsotope((G4int)j);
      G4ParticleHPIsoData isoData;
      isoData.Init(iso->GetN(), iso->GetZ(), iso->Getm(), 100., dirName, "/CrossSection");
      G4ParticleHPVector* hpVec = isoData.MakeChannelData();
      if (nullptr == hpVec || hpVec->GetVectorLength() == 0) continue;

      G4int len = hpVec->GetVectorLength();
      auto v = new G4PhysicsFreeVector(len, hpVec->GetX(0), hpVec->GetX(len - 1));
      for (G4int k = 0; k < len; ++k)
      {
        v->PutValues(k, hpVec->GetX(k), hpVec->GetY(k));
      }
      vecs[j] = v;
#ifdef G4VERBOSE
      if (fManager->GetDEBUG())
      {
        G4cout << "G4NeutronHPCaptureData: isotope table " << elm->GetName() << " "
               << iso->GetName() << " " << len << " points" << G4endl;
      }
#endif
    }
  }
}

void G4NeutronHPCaptureData::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  // the choice of the first instance
  if (fLock)
  {
    G4AutoLock l(&theHPCaptureData);
    if (fLock)
    {
      isFirst = true;
      fLock = false;
    }
    l.unlock();
  }
  if (!isFirst)
  {
    return;
  }
  if (p.GetParticleName() != "neutron")
  {
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only neutron is allowed";
    G4Exception("G4NeutronHPCaptureData::BuildPhysicsTable(..)", "had012", FatalException, ed, "");
    return;
  }

  // initialisation for the first instance, others are locked
  G4AutoLock l(&theHPCaptureData);
  if (theCrossSections != nullptr)
  {
    theCrossSections->clearAndDestroy();
    delete theCrossSections;
  }
  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  theCrossSections = new G4PhysicsTable(numberOfElements);

  // make a PhysicsVector for each element
  auto theElementTable = G4Element::GetElementTable();
  for (std::size_t i = 0; i < numberOfElements; ++i)
  {
    auto elm = (*theElementTable)[i];
#ifdef G4VERBOSE
    if (fManager->GetDEBUG())
    {
      G4cout << "ElementIndex " << elm->GetIndex() << "  " << elm->GetName() << G4endl;
    }
#endif
    G4PhysicsVector* physVec =
      G4ParticleHPData::Instance(G4Neutron::Neutron())->MakePhysicsVector(elm, this);
    theCrossSections->push_back(physVec);
  }
  fManager->RegisterCaptureCrossSections(theCrossSections);
  BuildIsotopeTables();
  l.unlock();
}

void G4NeutronHPCaptureData::DumpPhysicsTable(const G4ParticleDefinition&)
{
#ifdef G4VERBOSE
  if (fManager->GetVerboseLevel() == 0) return;

  //
  // Dump element based cross section
  // range 10e-5 eV to 20 MeV
  // 10 point per decade
  // in barn
  //

  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "Capture Cross Section of Neutron HP" << G4endl;
  G4cout << "(Pointwise cross-section at 0 Kelvin.)" << G4endl;
  G4cout << G4endl;
  G4cout << "Name of Element" << G4endl;
  G4cout << "Energy[eV]  XS[barn]" << G4endl;
  G4cout << G4endl;

  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  auto theElementTable = G4Element::GetElementTable();

  for (std::size_t i = 0; i < numberOfElements; ++i)
  {
    G4cout << (*theElementTable)[i]->GetName() << G4endl;
    G4cout << *((*theCrossSections)(i)) << G4endl;
  }
#endif
}

G4double G4NeutronHPCaptureData::GetCrossSection(const G4DynamicParticle* aP, const G4Element* anE,
                                                 G4double aT)
{
  auto idx = (G4int)anE->GetIndex();

  // prepare neutron
  G4double eKinetic = aP->GetKineticEnergy();
  if (eKinetic >= emax)
  {
    return 0.0;
  }

  G4int theA = anE->GetN();
  G4int theZ = anE->GetZasInt();
  G4double eleMass = G4NucleiProperties::GetNuclearMass(theA, theZ) / CLHEP::neutron_mass_c2;

  return DopplerAverage(aP, (*theCrossSections)(idx), eleMass, aT);
}

G4double G4NeutronHPCaptureData::DopplerAverage(const G4DynamicParticle* aP,
                                                const G4PhysicsVector* xsVec,
                                                G4double targetMass, G4double aT)
{
  G4double result = 0;
  G4double eKinetic = aP->GetKineticEnergy();

  // NEGLECT_DOPPLER
  if (fManager->GetNeglectDoppler())
  {
    return xsVec->Value(eKinetic);
  }

  G4ReactionProduct theNeutron(aP->GetDefinition());
  theNeutron.SetMomentum(aP->GetMomentum());
  theNeutron.SetKineticEnergy(eKinetic);

  // prepare thermal nucleus
  G4Nucleus aNuc;

  G4ReactionProduct boosted;
  G4double aXsection;

  // MC integration loop
  G4int counter = 0;
  G4double buffer = 0;
  G4int size = G4int(std::max(10., aT / 60 * kelvin));
  G4ThreeVector neutronVelocity =
    1. / G4Neutron::Neutron()->GetPDGMass() * theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();

  while (counter == 0 || std::abs(buffer - result / std::max(1, counter)) > 0.03 * buffer)
  // Loop checking, 11.05.2015, T. Koi
  {
    if (counter != 0) buffer = result / counter;
    while (counter < size)  // Loop checking, 11.05.2015, T. Koi
    {
      ++counter;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus(targetMass, aT);
      boosted.Lorentz(theNeutron, aThermalNuc);
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = xsVec->Value(theEkin);
      // velocity correction, or luminosity factor...
      G4ThreeVector targetVelocity = 1. / aThermalNuc.GetMass() * aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity - neutronVelocity).mag() / neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  return result;
}

void G4NeutronHPCaptureData::CrossSectionDescription(std::ostream& outF) const
{
  outF << "High Precision cross data based on Evaluated Nuclear Data Files"
       << " (ENDF) for radiative capture reaction of neutrons below 20 MeV";
}
