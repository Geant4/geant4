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

#include "G4ElementTable.hh"
#include "G4HadronicParameters.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleHPData.hh"
#include "G4ParticleHPManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

G4bool G4NeutronHPCaptureData::fLock = true;
G4PhysicsTable* G4NeutronHPCaptureData::theCrossSections = nullptr;

namespace
{
  G4Mutex theHPCaptureData = G4MUTEX_INITIALIZER;
}

G4NeutronHPCaptureData::G4NeutronHPCaptureData()
  : G4VCrossSectionDataSet("NeutronHPCaptureXS")
{
  emax = 20.*CLHEP::MeV;
  fManager = G4ParticleHPManager::GetInstance();
}

G4NeutronHPCaptureData::~G4NeutronHPCaptureData()
{
  if (isFirst) {
    if (nullptr != theCrossSections)
      theCrossSections->clearAndDestroy();
    delete theCrossSections;
    theCrossSections = nullptr;
  }
}

G4bool G4NeutronHPCaptureData::IsIsoApplicable(const G4DynamicParticle*,
                                               G4int, G4int,
                                               const G4Element*,
                                               const G4Material*)
{
  return true;
}

G4double G4NeutronHPCaptureData::GetIsoCrossSection(const G4DynamicParticle* dp,
                                                    G4int /*Z*/, G4int /*A*/,
                                                    const G4Isotope* /*iso*/,
                                                    const G4Element* element,
                                                    const G4Material* material)
{
  if (dp->GetKineticEnergy() == ke_cache && element == element_cache &&
      material == material_cache)
    return xs_cache;

  ke_cache = dp->GetKineticEnergy();
  element_cache = element;
  material_cache = material;
  G4double xs = GetCrossSection(dp, element, material->GetTemperature());
  xs_cache = xs;
  return xs;
}

void G4NeutronHPCaptureData::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  // the choice of the first instance
  if (fLock) {
    G4AutoLock l(&theHPCaptureData);
    if (fLock) {
      isFirst = true;
      fLock = false;
    }
    l.unlock();
  }
  if (!isFirst) { return; }
  if (p.GetParticleName() != "neutron") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only neutron is allowed";
    G4Exception("G4NeutronHPCaptureData::BuildPhysicsTable(..)","had012",
                FatalException, ed, "");
    return; 
  }

  // initialisation for the first instance, others are locked
  G4AutoLock l(&theHPCaptureData);
  if (theCrossSections != nullptr) {
    theCrossSections->clearAndDestroy();
    delete theCrossSections;
  }
  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  theCrossSections = new G4PhysicsTable(numberOfElements);
  
  // make a PhysicsVector for each element
  auto theElementTable = G4Element::GetElementTable();
  for (std::size_t i = 0; i < numberOfElements; ++i) {
    auto elm = (*theElementTable)[i];
#ifdef G4VERBOSE
    if (fManager->GetDEBUG()) {
      G4cout << "ElementIndex " << elm->GetIndex() << "  " 
             << elm->GetName() << G4endl;
    }
#endif
    G4PhysicsVector* physVec =
      G4ParticleHPData::Instance(G4Neutron::Neutron())->MakePhysicsVector(elm, this);
    theCrossSections->push_back(physVec);
  }
  fManager->RegisterCaptureCrossSections(theCrossSections);
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

  for (std::size_t i = 0; i < numberOfElements; ++i) {
    G4cout << (*theElementTable)[i]->GetName() << G4endl;
    G4cout << *((*theCrossSections)(i)) << G4endl;
  }
#endif
}

G4double G4NeutronHPCaptureData::GetCrossSection(const G4DynamicParticle* aP,
                                                 const G4Element* anE,
                                                 G4double aT)
{
  G4double result = 0;
  auto idx = (G4int)anE->GetIndex();

  // prepare neutron
  G4double eKinetic = aP->GetKineticEnergy();
  if (eKinetic >= emax) { return 0.0; }

  // NEGLECT_DOPPLER
  if (fManager->GetNeglectDoppler()) {
    return (*((*theCrossSections)(idx))).Value(eKinetic);
  }

  G4ReactionProduct theNeutron(aP->GetDefinition());
  theNeutron.SetMomentum(aP->GetMomentum());
  theNeutron.SetKineticEnergy(eKinetic);

  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4int theA = anE->GetN();
  G4int theZ = anE->GetZasInt();
  G4double eleMass = G4NucleiProperties::GetNuclearMass(theA, theZ)
    / CLHEP::neutron_mass_c2;

  G4ReactionProduct boosted;
  G4double aXsection;

  // MC integration loop
  G4int counter = 0;
  G4double buffer = 0;
  G4int size = G4int(std::max(10., aT / 60 * kelvin));
  G4ThreeVector neutronVelocity =
    1. / G4Neutron::Neutron()->GetPDGMass() * theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();

  while (counter == 0 ||
         std::abs(buffer - result / std::max(1, counter)) > 0.03 * buffer)
  // Loop checking, 11.05.2015, T. Koi
  {
    if (counter != 0) buffer = result / counter;
    while (counter < size)  // Loop checking, 11.05.2015, T. Koi
    {
      ++counter;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus(eleMass, aT);
      boosted.Lorentz(theNeutron, aThermalNuc);
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*((*theCrossSections)(idx))).Value(theEkin);
      // velocity correction, or luminosity factor...
      G4ThreeVector targetVelocity = 1. / aThermalNuc.GetMass() * aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity - neutronVelocity).mag() / neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  /*
    // Checking impact of  G4NEUTRONHP_NEGLECT_DOPPLER
    G4cout << " result " << result << " "
           << (*((*theCrossSections)(index))).Value(eKinetic) << " "
           << (*((*theCrossSections)(index))).Value(eKinetic) /result << G4endl;
  */
  return result;
}

void G4NeutronHPCaptureData::CrossSectionDescription(std::ostream& outF) const
{
  outF << "High Precision cross data based on Evaluated Nuclear Data Files"
       << " (ENDF) for radiative capture reaction of neutrons below 20 MeV";
}
