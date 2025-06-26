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
// G4ParticleHPThermalScatteringData
//
// 15-Nov-06 First implementation is done by T. Koi (SLAC/SCCS)
// 070625 implement clearCurrentXSData to fix memory leaking by T. Koi
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// ---------------------------------------------------------------------

#include "G4ParticleHPThermalScatteringData.hh"

#include "G4ElementTable.hh"
#include "G4Neutron.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

#include <algorithm>
#include <list>

G4ParticleHPThermalScatteringNames* G4ParticleHPThermalScatteringData::names = nullptr;
std::vector<G4int>* G4ParticleHPThermalScatteringData::indexOfThermalElement = nullptr;
std::map<std::pair<const G4Material*, const G4Element*>, G4int>* G4ParticleHPThermalScatteringData::dic = nullptr;

G4ParticleHPThermalScatteringData::G4ParticleHPThermalScatteringData()
  : G4VCrossSectionDataSet("NeutronHPThermalScatteringData")
{
  // Upper limit of neutron energy
  emax = 4 * CLHEP::eV;
  SetMinKinEnergy(0 * MeV);
  SetMaxKinEnergy(emax);

  ke_cache = 0.0;
  xs_cache = 0.0;
  element_cache = nullptr;
  material_cache = nullptr;

  if (nullptr == names) {
    isInitializer = true;
    indexOfThermalElement = new std::vector<G4int>;
    names = new G4ParticleHPThermalScatteringNames();
    dic = new std::map<std::pair<const G4Material*, const G4Element*>, G4int>;

    G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();
    coherent = new std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>;
    incoherent = new std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>;
    inelastic = new std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>;

    hpmanager->RegisterThermalScatteringCoherentCrossSections(coherent);
    hpmanager->RegisterThermalScatteringIncoherentCrossSections(incoherent);
    hpmanager->RegisterThermalScatteringInelasticCrossSections(inelastic);
  }
}

G4ParticleHPThermalScatteringData::~G4ParticleHPThermalScatteringData()
{
  if (!isInitializer) return;  

  clearCurrentXSData(coherent);
  clearCurrentXSData(incoherent);
  clearCurrentXSData(inelastic);

  delete names;
  delete dic;
  delete indexOfThermalElement;
  names = nullptr;
  dic = nullptr;
  indexOfThermalElement = nullptr;
}

G4bool G4ParticleHPThermalScatteringData::IsIsoApplicable(const G4DynamicParticle* dp, G4int /*Z*/,
                                                          G4int /*A*/, const G4Element* element,
                                                          const G4Material* material)
{
  if (dp->GetKineticEnergy() > emax)
    return false;

  if (dic->find(std::pair<const G4Material*, const G4Element*>((G4Material*)nullptr, element))
        != dic->end()
      || dic->find(std::pair<const G4Material*, const G4Element*>(material, element)) != dic->end())
    return true;

  return false;
}

G4double G4ParticleHPThermalScatteringData::GetIsoCrossSection(const G4DynamicParticle* dp,
                                                               G4int /*Z*/, G4int /*A*/,
                                                               const G4Isotope* /*iso*/,
                                                               const G4Element* element,
                                                               const G4Material* material)
{
  ke_cache = dp->GetKineticEnergy();
  element_cache = element;
  material_cache = material;
  G4double xs = GetCrossSection(dp, element, material);
  xs_cache = xs;
  return xs;
}

void G4ParticleHPThermalScatteringData::clearCurrentXSData(std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* ptr)
{
  if (nullptr == ptr) return;
  for (auto it = ptr->begin(); it != ptr->end(); ++it) {
    auto p = it->second;
    if (nullptr != p) {
      for (auto itt = p->begin(); itt != p->end(); ++itt) {
	delete itt->second;
      }
      delete p;
    }
  }
  delete ptr;
}

G4bool G4ParticleHPThermalScatteringData::IsApplicable(const G4DynamicParticle* aP,
                                                       const G4Element* anEle)
{
  // Check energy
  if (aP->GetKineticEnergy() > emax)
    return false;

  // anEle is one of Thermal elements
  auto ie = (G4int)anEle->GetIndex();
  for (auto const& it : *indexOfThermalElement) {
    if (ie == it) return true;
  }
  return false;
}

void G4ParticleHPThermalScatteringData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{
  if (&aP != G4Neutron::Neutron()) {
    G4ExceptionDescription ed;
    ed << "Neutron thermal scattering cannot be applied to " << aP.GetParticleName();

    G4Exception("G4ParticleHPThermalScatteringData::BuildPhysicsTable","hp0001",
		FatalException, ed, " run stopped");
    return;
  }

  // Common initialisation for all threads
  G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();
  verbose = hpmanager->GetVerboseLevel();

  coherent = hpmanager->GetThermalScatteringCoherentCrossSections();
  incoherent = hpmanager->GetThermalScatteringIncoherentCrossSections();
  inelastic = hpmanager->GetThermalScatteringInelasticCrossSections();

  // The initialisation is performed only in the master thread
  if (!isInitializer)
    return;

  std::map<G4String, G4int> co_dic;

  // Searching Materials
  auto const theMaterialTable = G4Material::GetMaterialTable();
  std::size_t numberOfMaterials = G4Material::GetNumberOfMaterials();
  for (std::size_t i = 0; i < numberOfMaterials; ++i) {
    G4Material* material = (*theMaterialTable)[i];
    auto numberOfElements = (G4int)material->GetNumberOfElements();
    for (G4int j = 0; j < numberOfElements; ++j) {
      const G4Element* element = material->GetElement(j);
      if (names->IsThisThermalElement(material->GetName(), element->GetName())) {
        G4int ts_ID_of_this_geometry;
        G4String ts_ndl_name = names->GetTS_NDL_Name(material->GetName(), element->GetName());
        if (co_dic.find(ts_ndl_name) != co_dic.cend()) {
          ts_ID_of_this_geometry = co_dic.find(ts_ndl_name)->second;
        }
        else {
          ts_ID_of_this_geometry = (G4int)co_dic.size();
          co_dic.insert(std::pair<G4String, G4int>(ts_ndl_name, ts_ID_of_this_geometry));
        }

        dic->insert(std::pair<std::pair<G4Material*, const G4Element*>, G4int>(
          std::pair<G4Material*, const G4Element*>(material, element), ts_ID_of_this_geometry));
      }
    }
  }

  // Searching TS Elements
  auto const theElementTable = G4Element::GetElementTable();
  std::size_t numberOfElements = G4Element::GetNumberOfElements();

  for (std::size_t i = 0; i < numberOfElements; ++i) {
    const G4Element* element = (*theElementTable)[i];
    if (names->IsThisThermalElement(element->GetName())) {
      G4int ts_ID_of_this_geometry;
      const G4String ts_ndl_name = names->GetTS_NDL_Name(element->GetName());
      if (co_dic.find(ts_ndl_name) != co_dic.cend()) {
	ts_ID_of_this_geometry = co_dic.find(ts_ndl_name)->second;
      }
      else {
	ts_ID_of_this_geometry = (G4int)co_dic.size();
	co_dic.insert(std::pair<G4String, G4int>(ts_ndl_name, ts_ID_of_this_geometry));
      }

      dic->insert(std::pair<std::pair<const G4Material*, const G4Element*>, G4int>(
	  std::pair<const G4Material*, const G4Element*>((G4Material*)nullptr, element),
          ts_ID_of_this_geometry));
    }
  }

  if (0 < verbose) {
    G4cout << "##T## Neutron HP Thermal Scattering Data: Following material-element pairs and/or elements "
      "are registered for " << dic->size() << " materials." << G4endl;

    if (dic->empty()) return;
    
    for (const auto& it : *dic) {
      if (it.first.first != nullptr) {
	G4cout << "    Material " << it.first.first->GetName() << " - Element "
	       << it.first.second->GetName() << ", internal thermal scattering id " << it.second
	       << G4endl;
      }
      else {
	G4cout << "    Element " << it.first.second->GetName() << ",  internal thermal scattering id "
	       << it.second << G4endl;
      }
    }
    G4cout << G4endl;
  }
  
  // Read Cross Section Data files
  const G4String dirName = hpmanager->GetNeutronHPPath() + "/ThermalScattering";

  G4String ndl_filename;
  G4String full_name;

  for (const auto& it : co_dic) {
    ndl_filename = it.first;
    G4int ts_ID = it.second;

    // Coherent
    full_name = dirName + "/Coherent/CrossSection/" + ndl_filename;
    auto coh_amapTemp_EnergyCross = readData(full_name);
    coherent->insert(std::pair<G4int, std::map<G4double, G4ParticleHPVector*>*>(
	ts_ID, coh_amapTemp_EnergyCross));

    // Incoherent
    full_name = dirName + "/Incoherent/CrossSection/" + ndl_filename;
    auto incoh_amapTemp_EnergyCross = readData(full_name);
    incoherent->insert(std::pair<G4int, std::map<G4double, G4ParticleHPVector*>*>(
        ts_ID, incoh_amapTemp_EnergyCross));

    // Inelastic
    full_name = dirName + "/Inelastic/CrossSection/" + ndl_filename;
    auto inela_amapTemp_EnergyCross = readData(full_name);
    inelastic->insert(std::pair<G4int, std::map<G4double, G4ParticleHPVector*>*>(
        ts_ID, inela_amapTemp_EnergyCross));
  }
}

std::map<G4double, G4ParticleHPVector*>*
G4ParticleHPThermalScatteringData::readData(const G4String& full_name)
{
  auto aData = new std::map<G4double, G4ParticleHPVector*>;

  std::istringstream theChannel;
  G4ParticleHPManager::GetInstance()->GetDataStream(full_name, theChannel);

  G4int dummy;
  while (theChannel >> dummy)  // MF // Loop checking, 11.05.2015, T. Koi
  {
    theChannel >> dummy;  // MT
    G4double temp;
    theChannel >> temp;
    auto anEnergyCross = new G4ParticleHPVector;
    G4int nData;
    theChannel >> nData;
    anEnergyCross->Init(theChannel, nData, eV, barn);
    aData->insert(std::pair<G4double, G4ParticleHPVector*>(temp, anEnergyCross));
  }

  return aData;
}

void G4ParticleHPThermalScatteringData::DumpPhysicsTable(const G4ParticleDefinition&)
{}

G4double G4ParticleHPThermalScatteringData::GetCrossSection(const G4DynamicParticle* aP,
                                                            const G4Element* anE,
                                                            const G4Material* aM)
{
  G4double result = 0;
  G4int ts_id = getTS_ID(aM, anE);
  
  if (ts_id == -1) return result;

  G4double aT = aM->GetTemperature();

  auto u = coherent->find(ts_id);
  G4double Xcoh = (u != coherent->end()) ? GetX(aP, aT, u->second) : 0.0;
  auto v = incoherent->find(ts_id);
  G4double Xincoh = (v != incoherent->end()) ? GetX(aP, aT, v->second) : 0.0;
  auto w = inelastic->find(ts_id);
  G4double Xinela = (w != inelastic->end()) ? GetX(aP, aT, w->second) : 0.0;

  result = Xcoh + Xincoh + Xinela;

  return result;
}

G4double G4ParticleHPThermalScatteringData::GetInelasticCrossSection(const G4DynamicParticle* aP,
                                                                     const G4Element* anE,
                                                                     const G4Material* aM)
{
  G4double result = 0;
  G4int ts_id = getTS_ID(aM, anE);
  if (ts_id == -1) return result;
  G4double aT = aM->GetTemperature();
  auto ptr = inelastic->find(ts_id);
  if (ptr != inelastic->end()) {
    result = GetX(aP, aT, ptr->second);
  }
  return result;
}

G4double G4ParticleHPThermalScatteringData::GetCoherentCrossSection(const G4DynamicParticle* aP,
                                                                    const G4Element* anE,
                                                                    const G4Material* aM)
{
  G4double result = 0;
  G4int ts_id = getTS_ID(aM, anE);
  if (ts_id == -1) return result;
  G4double aT = aM->GetTemperature();
  auto u = coherent->find(ts_id);
  if (u != coherent->end()) { result = GetX(aP, aT, u->second); }
  return result;
}

G4double G4ParticleHPThermalScatteringData::GetIncoherentCrossSection(const G4DynamicParticle* aP,
                                                                      const G4Element* anE,
                                                                      const G4Material* aM)
{
  G4double result = 0;
  G4int ts_id = getTS_ID(aM, anE);
  if (ts_id == -1) return result;
  G4double aT = aM->GetTemperature();
  auto u = incoherent->find(ts_id);
  if (u != incoherent->end()) { result = GetX(aP, aT, u->second); }
  return result;
}

G4int G4ParticleHPThermalScatteringData::getTS_ID(const G4Material* material,
                                                  const G4Element* element)
{
  auto it = dic->find(std::pair<const G4Material*, const G4Element*>((G4Material*)nullptr, element));
  if (it != dic->end()) { return it->second; }
  auto jt = dic->find(std::pair<const G4Material*, const G4Element*>(material, element));
  if (jt != dic->end()) { return jt->second; }
  return -1;
}

G4double G4ParticleHPThermalScatteringData::
GetX(const G4DynamicParticle* aP, G4double aT,
     std::map<G4double, G4ParticleHPVector*>* amapTemp_EnergyCross)
{
  if (amapTemp_EnergyCross->empty())
    return 0.0;

  G4double eKinetic = aP->GetKineticEnergy();
  auto it_begin = amapTemp_EnergyCross->begin();
  G4double Tmin = it_begin->first;
  std::size_t n = amapTemp_EnergyCross->size();

  // special cases
  if (n == 1 || aT <= Tmin) {
    return it_begin->second->GetXsec(eKinetic);
  }

  // high temperature
  auto it_end = amapTemp_EnergyCross->end();
  --it_end;
  G4double Tmax = it_end->first;
  if (aT >= Tmax) {
    return it_end->second->GetXsec(eKinetic);
  }

  // linear interpolation between two temperature values
  ++it_begin;
  auto it = it_begin;
  G4double TH = Tmin;
  for (;;) {
    TH = it->first;
    if (aT <= TH || it == it_end) break;
    ++it;
  }

  G4double XH = it->second->GetXsec(eKinetic);
  --it;
  G4double TL = it->first;
  if (TH == TL) return XH;
  
  G4double XL = it->second->GetXsec(eKinetic);
  return (aT - TL) * (XH - XL) / (TH - TL) + XL;  
}

void G4ParticleHPThermalScatteringData::AddUserThermalScatteringFile(const G4String& nameG4Element,
                                                                     const G4String& filename)
{
  names->AddThermalElement(nameG4Element, filename);
}

void G4ParticleHPThermalScatteringData::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "High Precision cross data based on thermal scattering data in evaluated nuclear data "
             "libraries for neutrons below 5eV on specific materials\n";
}
