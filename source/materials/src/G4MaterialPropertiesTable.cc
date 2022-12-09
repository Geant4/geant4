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
//
//
////////////////////////////////////////////////////////////////////////
// G4MaterialPropertiesTable Implementation
////////////////////////////////////////////////////////////////////////
//
// File: G4MaterialPropertiesTable.cc
// Version:     1.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     2005-05-12 add SetGROUPVEL(), courtesy of
//              Horton-Smith (bug report #741), by P. Gumplinger
//              2002-11-05 add named material constants by P. Gumplinger
//              1999-11-05 Migration from G4RWTPtrHashDictionary to STL
//                         by John Allison
//              1997-03-26 by Peter Gumplinger
//              > cosmetics (only)
//
////////////////////////////////////////////////////////////////////////

#include "G4MaterialPropertiesTable.hh"

#include "globals.hh"
#include "G4Log.hh"
#include "G4OpticalMaterialProperties.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <cassert>

#ifdef G4MULTITHREADED
#  include "G4AutoLock.hh"
namespace
{
  G4Mutex materialPropertyTableMutex = G4MUTEX_INITIALIZER;
}
#endif

G4MaterialPropertiesTable::G4MaterialPropertiesTable()
{
  // elements of these 2 vectors must be in same order as
  // the corresponding enums in G4MaterialPropertiesIndex.hh
  fMatPropNames.emplace_back("RINDEX");
  fMatPropNames.emplace_back("REFLECTIVITY");
  fMatPropNames.emplace_back("REALRINDEX");
  fMatPropNames.emplace_back("IMAGINARYRINDEX");
  fMatPropNames.emplace_back("EFFICIENCY");
  fMatPropNames.emplace_back("TRANSMITTANCE");
  fMatPropNames.emplace_back("SPECULARLOBECONSTANT");
  fMatPropNames.emplace_back("SPECULARSPIKECONSTANT");
  fMatPropNames.emplace_back("BACKSCATTERCONSTANT");
  fMatPropNames.emplace_back("GROUPVEL");
  fMatPropNames.emplace_back("MIEHG");
  fMatPropNames.emplace_back("RAYLEIGH");
  fMatPropNames.emplace_back("WLSCOMPONENT");
  fMatPropNames.emplace_back("WLSABSLENGTH");
  fMatPropNames.emplace_back("WLSCOMPONENT2");
  fMatPropNames.emplace_back("WLSABSLENGTH2");
  fMatPropNames.emplace_back("ABSLENGTH");
  fMatPropNames.emplace_back("PROTONSCINTILLATIONYIELD");
  fMatPropNames.emplace_back("DEUTERONSCINTILLATIONYIELD");
  fMatPropNames.emplace_back("TRITONSCINTILLATIONYIELD");
  fMatPropNames.emplace_back("ALPHASCINTILLATIONYIELD");
  fMatPropNames.emplace_back("IONSCINTILLATIONYIELD");
  fMatPropNames.emplace_back("ELECTRONSCINTILLATIONYIELD");
  fMatPropNames.emplace_back("SCINTILLATIONCOMPONENT1");
  fMatPropNames.emplace_back("SCINTILLATIONCOMPONENT2");
  fMatPropNames.emplace_back("SCINTILLATIONCOMPONENT3");
  fMatPropNames.emplace_back("COATEDRINDEX");

  assert(fMatPropNames.size() == kNumberOfPropertyIndex);

  fMP.assign(kNumberOfPropertyIndex, nullptr);

  fMatConstPropNames.emplace_back("SURFACEROUGHNESS");
  fMatConstPropNames.emplace_back("ISOTHERMAL_COMPRESSIBILITY");
  fMatConstPropNames.emplace_back("RS_SCALE_FACTOR");
  fMatConstPropNames.emplace_back("WLSMEANNUMBERPHOTONS");
  fMatConstPropNames.emplace_back("WLSTIMECONSTANT");
  fMatConstPropNames.emplace_back("WLSMEANNUMBERPHOTONS2");
  fMatConstPropNames.emplace_back("WLSTIMECONSTANT2");
  fMatConstPropNames.emplace_back("MIEHG_FORWARD");
  fMatConstPropNames.emplace_back("MIEHG_BACKWARD");
  fMatConstPropNames.emplace_back("MIEHG_FORWARD_RATIO");
  fMatConstPropNames.emplace_back("SCINTILLATIONYIELD");
  fMatConstPropNames.emplace_back("RESOLUTIONSCALE");
  fMatConstPropNames.emplace_back("FERMIPOT");
  fMatConstPropNames.emplace_back("DIFFUSION");
  fMatConstPropNames.emplace_back("SPINFLIP");
  fMatConstPropNames.emplace_back("LOSS");
  fMatConstPropNames.emplace_back("LOSSCS");
  fMatConstPropNames.emplace_back("ABSCS");
  fMatConstPropNames.emplace_back("SCATCS");
  fMatConstPropNames.emplace_back("MR_NBTHETA");
  fMatConstPropNames.emplace_back("MR_NBE");
  fMatConstPropNames.emplace_back("MR_RRMS");
  fMatConstPropNames.emplace_back("MR_CORRLEN");
  fMatConstPropNames.emplace_back("MR_THETAMIN");
  fMatConstPropNames.emplace_back("MR_THETAMAX");
  fMatConstPropNames.emplace_back("MR_EMIN");
  fMatConstPropNames.emplace_back("MR_EMAX");
  fMatConstPropNames.emplace_back("MR_ANGNOTHETA");
  fMatConstPropNames.emplace_back("MR_ANGNOPHI");
  fMatConstPropNames.emplace_back("MR_ANGCUT");
  fMatConstPropNames.emplace_back("SCINTILLATIONTIMECONSTANT1");
  fMatConstPropNames.emplace_back("SCINTILLATIONTIMECONSTANT2");
  fMatConstPropNames.emplace_back("SCINTILLATIONTIMECONSTANT3");
  fMatConstPropNames.emplace_back("SCINTILLATIONRISETIME1");
  fMatConstPropNames.emplace_back("SCINTILLATIONRISETIME2");
  fMatConstPropNames.emplace_back("SCINTILLATIONRISETIME3");
  fMatConstPropNames.emplace_back("SCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("SCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("SCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("PROTONSCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("PROTONSCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("PROTONSCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("DEUTERONSCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("DEUTERONSCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("DEUTERONSCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("TRITONSCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("TRITONSCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("TRITONSCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("ALPHASCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("ALPHASCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("ALPHASCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("IONSCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("IONSCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("IONSCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("ELECTRONSCINTILLATIONYIELD1");
  fMatConstPropNames.emplace_back("ELECTRONSCINTILLATIONYIELD2");
  fMatConstPropNames.emplace_back("ELECTRONSCINTILLATIONYIELD3");
  fMatConstPropNames.emplace_back("COATEDTHICKNESS");
  fMatConstPropNames.emplace_back("COATEDFRUSTRATEDTRANSMISSION");

  assert(fMatConstPropNames.size() == kNumberOfConstPropertyIndex);

  fMCP.assign(kNumberOfConstPropertyIndex, { 0., false });
}

G4MaterialPropertiesTable::~G4MaterialPropertiesTable()
{
  for(auto prop : fMP)
  {
    delete(prop);
  }
}

G4int G4MaterialPropertiesTable::GetConstPropertyIndex(
  const G4String& key) const
{
  // Returns the constant material property index corresponding to a key

  std::size_t index = std::distance(
    fMatConstPropNames.cbegin(),
    std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key));
  if(index < fMatConstPropNames.size())
  {
    return (G4int)index;
  }

  G4ExceptionDescription ed;
  ed << "Constant Material Property Index for key " << key << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetConstPropertyIndex()", "mat200",
              FatalException, ed);
  return 0;
}

G4int G4MaterialPropertiesTable::GetPropertyIndex(const G4String& key) const
{
  // Returns the material property index corresponding to a key
  std::size_t index =
    std::distance(fMatPropNames.cbegin(),
                  std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key));
  if(index < fMatPropNames.size())
  {
    return (G4int)index;
  }
  G4ExceptionDescription ed;
  ed << "Material Property Index for key " << key << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetPropertyIndex()", "mat201",
              FatalException, ed);
  return 0;
}

G4double G4MaterialPropertiesTable::GetConstProperty(const G4int index) const
{
  // Returns the constant material property corresponding to an index
  // fatal exception if property not found

  if(index < (G4int) fMCP.size() && fMCP[index].second)
  {
    return fMCP[index].first;
  }
  G4ExceptionDescription ed;
  ed << "Constant Material Property Index " << index << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetConstProperty()", "mat202",
              FatalException, ed);
  return 0.;
}

G4double G4MaterialPropertiesTable::GetConstProperty(const G4String& key) const
{
  // Returns the constant material property corresponding to a key
  // fatal exception if property not found

  return GetConstProperty(GetConstPropertyIndex(key));
}

G4double G4MaterialPropertiesTable::GetConstProperty(const char* key) const
{
  return GetConstProperty(GetConstPropertyIndex(G4String(key)));
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4int index) const
{
  // Returns true if a const property corresponding to 'index' exists

  return index >= 0 && index < (G4int) fMCP.size() && fMCP[index].second;
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4String& key) const
{
  // Returns true if a const property 'key' exists
  std::size_t index = std::distance(
    fMatConstPropNames.cbegin(),
    std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key));
  if(index < fMatConstPropNames.size())
  {  // index is type std::size_t so >= 0
    return ConstPropertyExists((G4int)index);
  }
  return false;
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const char* key) const
{
  std::size_t index = std::distance(
    fMatConstPropNames.cbegin(),
    std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key));
  if(index < fMatConstPropNames.size())
  {  // index is type std::size_t so >= 0
    return ConstPropertyExists((G4int)index);
  }
  return false;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(
  const G4String& key) const
{
  // Returns a Material Property Vector corresponding to a key
  if(std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) !=
     fMatPropNames.cend())
  {
    const G4int index = GetPropertyIndex(G4String(key));
    return GetProperty(index);
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(
  const char* key) const
{
  if(std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) !=
     fMatPropNames.cend())
  {
    const G4int index = GetPropertyIndex(G4String(key));
    return GetProperty(index);
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(
  const G4int index) const
{
  // Returns a Material Property Vector corresponding to an index
  // returns nullptr if the property has not been defined by user
  if(index >= 0 && index < (G4int) fMP.size())
  {
    return fMP[index];
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::AddProperty(
  const G4String& key, const std::vector<G4double>& photonEnergies,
  const std::vector<G4double>& propertyValues, G4bool createNewKey,
  G4bool spline)
{
  if(photonEnergies.size() != propertyValues.size())
  {
    G4ExceptionDescription ed;
    ed << "AddProperty error!";
    G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat204",
                FatalException, ed);
  }

  // G4PhysicsVector assumes energies are in increasing order
  for (std::size_t i = 0; i < photonEnergies.size() - 1; ++i)
  {
    if(photonEnergies.at(i+1) < photonEnergies.at(i))
    {
      G4ExceptionDescription ed;
      ed << "Energies in material property table must be in increasing "
         << "order. Key: " << key << " Energy: " << photonEnergies.at(i+1);
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat215",
                  FatalException, ed);
    }
  }

  // if the key doesn't exist, add it if requested
  if(std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) ==
     fMatPropNames.cend())
  {
    if(createNewKey)
    {
      fMatPropNames.push_back(key);
      fMP.push_back(nullptr);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Attempting to create a new material property key " << key
         << " without setting\n"
         << "createNewKey parameter of AddProperty to true.";
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat205",
                  FatalException, ed);
    }
  }

  auto* mpv =
    new G4MaterialPropertyVector(photonEnergies, propertyValues, spline);
  mpv->SetVerboseLevel(1);
  if(spline)
  {
    mpv->FillSecondDerivatives();
  }
  G4int index = GetPropertyIndex(key);
  fMP[index]  = mpv;

  // if key is RINDEX, we calculate GROUPVEL -
  // contribution from Tao Lin (IHEP, the JUNO experiment)
  if(key == "RINDEX")
  {
    CalculateGROUPVEL();
  }

  return mpv;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::AddProperty(
  const char* key, G4double* photonEnergies, G4double* propertyValues,
  G4int numEntries, G4bool createNewKey, G4bool spline)
{

  // Provides a way of adding a property to the Material Properties
  // Table given a pair of arrays and a key
  G4String k(key);

  std::vector<G4double> energies(photonEnergies, photonEnergies + numEntries);
  std::vector<G4double> values(propertyValues, propertyValues + numEntries);
  return AddProperty(k, energies, values, createNewKey, spline);
}

void G4MaterialPropertiesTable::AddProperty(const G4String& key,
                                            G4MaterialPropertyVector* mpv,
                                            G4bool createNewKey)
{
  //  Provides a way of adding a property to the Material Properties
  //  Table given an G4MaterialPropertyVector Reference and a key

  // G4PhysicsVector assumes energies are in increasing order
  // An MPV with size==0 or 1 is also ok
  if (mpv->GetVectorLength() > 1) {
    for (std::size_t i = 0; i < mpv->GetVectorLength() - 1; ++i)
    {
      if(mpv->Energy(i+1) < mpv->Energy(i))
      {
        G4ExceptionDescription ed;
        ed << "Energies in material property vector must be in increasing "
           << "order. Key: " << key << " Energy: " << mpv->Energy(i+1);
        G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat216",
                    FatalException, ed);
      }
    }
  }

  // if the key doesn't exist, add it
  if(std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) ==
     fMatPropNames.cend())
  {
    if(createNewKey)
    {
      fMatPropNames.push_back(key);
      fMP.push_back(nullptr);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Attempting to create a new material property key " << key
         << " without setting\n"
         << "createNewKey parameter of AddProperty to true.";
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat206",
                  FatalException, ed);
    }
  }
  G4int index = GetPropertyIndex(key);
  fMP[index]  = mpv;

  // if key is RINDEX, we calculate GROUPVEL -
  // contribution from Tao Lin (IHEP, the JUNO experiment)
  if(key == "RINDEX")
  {
    CalculateGROUPVEL();
  }
}

void G4MaterialPropertiesTable::AddProperty(const char* key,
                                            G4MaterialPropertyVector* mpv,
                                            G4bool createNewKey)
{
  AddProperty(G4String(key), mpv, createNewKey);
}

void G4MaterialPropertiesTable::AddProperty(const G4String& key,
                                            const G4String& mat)
{
  // load a material property vector defined in Geant4 source
  G4MaterialPropertyVector* v =
    G4OpticalMaterialProperties::GetProperty(key, mat);
  AddProperty(key, v);
}

void G4MaterialPropertiesTable::AddConstProperty(const G4String& key,
                                                 G4double propertyValue,
                                                 G4bool createNewKey)
{
  // Provides a way of adding a constant property to the Material Properties
  // Table given a key
  if(std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key) ==
     fMatConstPropNames.cend())
  {
    if(createNewKey)
    {
      fMatConstPropNames.push_back(key);
      fMCP.emplace_back(0., true);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Attempting to create a new material constant property key " << key
         << " without setting"
         << " createNewKey parameter of AddProperty to true.";
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat207",
                  FatalException, ed);
    }
  }
  G4int index = GetConstPropertyIndex(key);

  fMCP[index] = std::pair<G4double, G4bool>{ propertyValue, true };
}

void G4MaterialPropertiesTable::AddConstProperty(const char* key,
                                                 G4double propertyValue,
                                                 G4bool createNewKey)
{
  // Provides a way of adding a constant property to the Material Properties
  // Table given a key
  AddConstProperty(G4String(key), propertyValue, createNewKey);
}

void G4MaterialPropertiesTable::RemoveConstProperty(const G4String& key)
{
  G4int index = GetConstPropertyIndex(key);
  if(index < (G4int) fMCP.size())
  {
    fMCP[index] = std::pair<G4double, G4bool>{ 0., false };
  }
}

void G4MaterialPropertiesTable::RemoveConstProperty(const char* key)
{
  RemoveConstProperty(G4String(key));
}

void G4MaterialPropertiesTable::RemoveProperty(const G4String& key)
{
  G4int index = GetPropertyIndex(key);
  delete fMP[index];
  fMP[index] = nullptr;
}

void G4MaterialPropertiesTable::RemoveProperty(const char* key)
{
  RemoveProperty(G4String(key));
}

void G4MaterialPropertiesTable::AddEntry(const G4String& key,
                                         G4double aPhotonEnergy,
                                         G4double aPropertyValue)
{
  // Allows to add an entry pair directly to the Material Property Vector
  // given a key.
  if(std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) ==
     fMatPropNames.cend())
  {
    G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat214",
                FatalException, "Material Property Vector not found.");
  }
  G4int index = GetPropertyIndex(key);

  G4MaterialPropertyVector* targetVector = fMP[index];
  if(targetVector != nullptr)
  {
    // do not allow duplicate energies
    for (std::size_t i = 0; i < targetVector->GetVectorLength(); ++i)
    {
      if(aPhotonEnergy == targetVector->Energy(i))
      {
        G4ExceptionDescription ed;
        ed << "Energy values in material property vector must be unique. "
           << "Key: " << key;
        G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat217",
                    FatalException, ed);
      }
    }

    targetVector->InsertValues(aPhotonEnergy, aPropertyValue);
  }
  else
  {
    G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat208",
                FatalException, "Material Property Vector not found.");
  }
  if(key == "RINDEX")
  {
    CalculateGROUPVEL();
  }
}

void G4MaterialPropertiesTable::AddEntry(const char* key,
                                         G4double aPhotonEnergy,
                                         G4double aPropertyValue)
{
  AddEntry(G4String(key), aPhotonEnergy, aPropertyValue);
}

void G4MaterialPropertiesTable::DumpTable() const
{
  // material properties
  G4int j = 0;
  for(const auto& prop : fMP)
  {
    if(prop != nullptr)
    {
      G4cout << j << ": " << fMatPropNames[j] << G4endl;
      prop->DumpValues();
    }
    ++j;
  }
  // material constant properties
  j = 0;
  for(const auto& cprop : fMCP)
  {
    if(cprop.second)
    {
      G4cout << j << ": " << fMatConstPropNames[j] << " " << cprop.first
             << G4endl;
    }
    ++j;
  }
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::CalculateGROUPVEL()
{
#ifdef G4MULTITHREADED
  G4AutoLock mptm(&materialPropertyTableMutex);
#endif

  // check if "GROUPVEL" already exists. If so, remove it.
  if(fMP[kGROUPVEL] != nullptr)
  {
    this->RemoveProperty("GROUPVEL");
  }

  // fetch RINDEX data, give up if unavailable
  G4MaterialPropertyVector* rindex = this->GetProperty(kRINDEX);
  if(rindex == nullptr)
  {
    return nullptr;
  }

  // RINDEX exists but has no entries, give up
  if(rindex->GetVectorLength() == 0)
  {
    return nullptr;
  }

  // add GROUPVEL vector
  auto* groupvel = new G4MaterialPropertyVector();
  groupvel->SetVerboseLevel(1);

  // fill GROUPVEL vector using RINDEX values
  // rindex built-in "iterator" was advanced to first entry above
  G4double E0 = rindex->Energy(0);
  G4double n0 = (*rindex)[0];

  if(E0 <= 0.)
  {
    G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat211",
                FatalException, "Optical Photon Energy <= 0");
  }

  if(rindex->GetVectorLength() >= 2)
  {
    // good, we have at least two entries in RINDEX
    // get next energy/value pair

    G4double E1 = rindex->Energy(1);
    G4double n1 = (*rindex)[1];

    if(E1 <= 0.)
    {
      G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat212",
                  FatalException, "Optical Photon Energy <= 0");
    }

    G4double vg;

    // add entry at first photon energy
    vg = c_light / (n0 + (n1 - n0) / G4Log(E1 / E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    if((vg < 0) || (vg > c_light / n0))
    {
      vg = c_light / n0;
    }

    groupvel->InsertValues(E0, vg);

    // add entries at midpoints between remaining photon energies
    for(std::size_t i = 2; i < rindex->GetVectorLength(); ++i)
    {
      vg = c_light / (0.5 * (n0 + n1) + (n1 - n0) / G4Log(E1 / E0));

      // allow only for 'normal dispersion' -> dn/d(logE) > 0
      if((vg < 0) || (vg > c_light / (0.5 * (n0 + n1))))
      {
        vg = c_light / (0.5 * (n0 + n1));
      }
      groupvel->InsertValues(0.5 * (E0 + E1), vg);

      // get next energy/value pair, or exit loop
      E0 = E1;
      n0 = n1;
      E1 = rindex->Energy(i);
      n1 = (*rindex)[i];

      if(E1 <= 0.)
      {
        G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat213",
                    FatalException, "Optical Photon Energy <= 0");
      }
    }

    // add entry at last photon energy
    vg = c_light / (n1 + (n1 - n0) / G4Log(E1 / E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    if((vg < 0) || (vg > c_light / n1))
    {
      vg = c_light / n1;
    }
    groupvel->InsertValues(E1, vg);
  }
  else  // only one entry in RINDEX -- weird!
  {
    groupvel->InsertValues(E0, c_light / n0);
  }

  this->AddProperty("GROUPVEL", groupvel);

  return groupvel;
}
