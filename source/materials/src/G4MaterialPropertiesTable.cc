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

#include "G4Log.hh"
#include "G4OpticalMaterialProperties.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"

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
  fMatPropNames.assign(kNumberOfPropertyIndex, "");
  fMatPropNames[kRINDEX] =                    "RINDEX";
  fMatPropNames[kREFLECTIVITY] =              "REFLECTIVITY";
  fMatPropNames[kREALRINDEX] =                "REALRINDEX";
  fMatPropNames[kIMAGINARYRINDEX] =           "IMAGINARYRINDEX";
  fMatPropNames[kEFFICIENCY] =                "EFFICIENCY";
  fMatPropNames[kTRANSMITTANCE] =             "TRANSMITTANCE";
  fMatPropNames[kSPECULARLOBECONSTANT] =      "SPECULARLOBECONSTANT";
  fMatPropNames[kSPECULARSPIKECONSTANT] =     "SPECULARSPIKECONSTANT";
  fMatPropNames[kBACKSCATTERCONSTANT] =       "BACKSCATTERCONSTANT";
  fMatPropNames[kGROUPVEL] =                  "GROUPVEL";
  fMatPropNames[kMIEHG] =                     "MIEHG";
  fMatPropNames[kRAYLEIGH] =                  "RAYLEIGH";
  fMatPropNames[kWLSCOMPONENT] =              "WLSCOMPONENT";
  fMatPropNames[kWLSABSLENGTH] =              "WLSABSLENGTH";
  fMatPropNames[kWLSCOMPONENT2] =             "WLSCOMPONENT2";
  fMatPropNames[kWLSABSLENGTH2] =             "WLSABSLENGTH2";
  fMatPropNames[kABSLENGTH] =                 "ABSLENGTH";
  fMatPropNames[kPROTONSCINTILLATIONYIELD] =  "PROTONSCINTILLATIONYIELD";
  fMatPropNames[kDEUTERONSCINTILLATIONYIELD] ="DEUTERONSCINTILLATIONYIELD";
  fMatPropNames[kTRITONSCINTILLATIONYIELD] =  "TRITONSCINTILLATIONYIELD";
  fMatPropNames[kALPHASCINTILLATIONYIELD] =   "ALPHASCINTILLATIONYIELD";
  fMatPropNames[kIONSCINTILLATIONYIELD] =     "IONSCINTILLATIONYIELD";
  fMatPropNames[kELECTRONSCINTILLATIONYIELD] ="ELECTRONSCINTILLATIONYIELD";
  fMatPropNames[kSCINTILLATIONCOMPONENT1] =   "SCINTILLATIONCOMPONENT1";
  fMatPropNames[kSCINTILLATIONCOMPONENT2] =   "SCINTILLATIONCOMPONENT2";
  fMatPropNames[kSCINTILLATIONCOMPONENT3] =   "SCINTILLATIONCOMPONENT3";
  fMatPropNames[kCOATEDRINDEX] =              "COATEDRINDEX";

  fMP.assign(kNumberOfPropertyIndex, nullptr);

  fMatConstPropNames.assign(kNumberOfConstPropertyIndex, "");
  fMatConstPropNames[kSURFACEROUGHNESS] =             "SURFACEROUGHNESS";
  fMatConstPropNames[kISOTHERMAL_COMPRESSIBILITY] =   "ISOTHERMAL_COMPRESSIBILITY";
  fMatConstPropNames[kRS_SCALE_FACTOR] =              "RS_SCALE_FACTOR";
  fMatConstPropNames[kWLSMEANNUMBERPHOTONS] =         "WLSMEANNUMBERPHOTONS";
  fMatConstPropNames[kWLSTIMECONSTANT] =              "WLSTIMECONSTANT";
  fMatConstPropNames[kWLSMEANNUMBERPHOTONS2] =        "WLSMEANNUMBERPHOTONS2";
  fMatConstPropNames[kWLSTIMECONSTANT2] =             "WLSTIMECONSTANT2";
  fMatConstPropNames[kMIEHG_FORWARD] =                "MIEHG_FORWARD";
  fMatConstPropNames[kMIEHG_BACKWARD] =               "MIEHG_BACKWARD";
  fMatConstPropNames[kMIEHG_FORWARD_RATIO] =          "MIEHG_FORWARD_RATIO";
  fMatConstPropNames[kSCINTILLATIONYIELD] =           "SCINTILLATIONYIELD";
  fMatConstPropNames[kRESOLUTIONSCALE] =              "RESOLUTIONSCALE";
  fMatConstPropNames[kFERMIPOT] =                     "FERMIPOT";
  fMatConstPropNames[kDIFFUSION] =                    "DIFFUSION";
  fMatConstPropNames[kSPINFLIP] =                     "SPINFLIP";
  fMatConstPropNames[kLOSS] =                         "LOSS";
  fMatConstPropNames[kLOSSCS] =                       "LOSSCS";
  fMatConstPropNames[kABSCS] =                        "ABSCS";
  fMatConstPropNames[kSCATCS] =                       "SCATCS";
  fMatConstPropNames[kMR_NBTHETA] =                   "MR_NBTHETA";
  fMatConstPropNames[kMR_NBE] =                       "MR_NBE";
  fMatConstPropNames[kMR_RRMS] =                      "MR_RRMS";
  fMatConstPropNames[kMR_CORRLEN] =                   "MR_CORRLEN";
  fMatConstPropNames[kMR_THETAMIN] =                  "MR_THETAMIN";
  fMatConstPropNames[kMR_THETAMAX] =                  "MR_THETAMAX";
  fMatConstPropNames[kMR_EMIN] =                      "MR_EMIN";
  fMatConstPropNames[kMR_EMAX] =                      "MR_EMAX";
  fMatConstPropNames[kMR_ANGNOTHETA] =                "MR_ANGNOTHETA";
  fMatConstPropNames[kMR_ANGNOPHI] =                  "MR_ANGNOPHI";
  fMatConstPropNames[kMR_ANGCUT] =                    "MR_ANGCUT";
  fMatConstPropNames[kSCINTILLATIONTIMECONSTANT1] =   "SCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kSCINTILLATIONTIMECONSTANT2] =   "SCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kSCINTILLATIONTIMECONSTANT3] =   "SCINTILLATIONTIMECONSTANT3";
  fMatConstPropNames[kSCINTILLATIONRISETIME1] =       "SCINTILLATIONRISETIME1";
  fMatConstPropNames[kSCINTILLATIONRISETIME2] =       "SCINTILLATIONRISETIME2";
  fMatConstPropNames[kSCINTILLATIONRISETIME3] =       "SCINTILLATIONRISETIME3";
  fMatConstPropNames[kSCINTILLATIONYIELD1] =          "SCINTILLATIONYIELD1";
  fMatConstPropNames[kSCINTILLATIONYIELD2] =          "SCINTILLATIONYIELD2";
  fMatConstPropNames[kSCINTILLATIONYIELD3] =          "SCINTILLATIONYIELD3";
  fMatConstPropNames[kPROTONSCINTILLATIONYIELD1] =    "PROTONSCINTILLATIONYIELD1";
  fMatConstPropNames[kPROTONSCINTILLATIONYIELD2] =    "PROTONSCINTILLATIONYIELD2";
  fMatConstPropNames[kPROTONSCINTILLATIONYIELD3] =    "PROTONSCINTILLATIONYIELD3";
  fMatConstPropNames[kDEUTERONSCINTILLATIONYIELD1] =  "DEUTERONSCINTILLATIONYIELD1";
  fMatConstPropNames[kDEUTERONSCINTILLATIONYIELD2] =  "DEUTERONSCINTILLATIONYIELD2";
  fMatConstPropNames[kDEUTERONSCINTILLATIONYIELD3] =  "DEUTERONSCINTILLATIONYIELD3";
  fMatConstPropNames[kTRITONSCINTILLATIONYIELD1] =    "TRITONSCINTILLATIONYIELD1";
  fMatConstPropNames[kTRITONSCINTILLATIONYIELD2] =    "TRITONSCINTILLATIONYIELD2";
  fMatConstPropNames[kTRITONSCINTILLATIONYIELD3] =    "TRITONSCINTILLATIONYIELD3";
  fMatConstPropNames[kALPHASCINTILLATIONYIELD1] =     "ALPHASCINTILLATIONYIELD1";
  fMatConstPropNames[kALPHASCINTILLATIONYIELD2] =     "ALPHASCINTILLATIONYIELD2";
  fMatConstPropNames[kALPHASCINTILLATIONYIELD3] =     "ALPHASCINTILLATIONYIELD3";
  fMatConstPropNames[kIONSCINTILLATIONYIELD1] =       "IONSCINTILLATIONYIELD1";
  fMatConstPropNames[kIONSCINTILLATIONYIELD2] =       "IONSCINTILLATIONYIELD2";
  fMatConstPropNames[kIONSCINTILLATIONYIELD3] =       "IONSCINTILLATIONYIELD3";
  fMatConstPropNames[kELECTRONSCINTILLATIONYIELD1] =  "ELECTRONSCINTILLATIONYIELD1";
  fMatConstPropNames[kELECTRONSCINTILLATIONYIELD2] =  "ELECTRONSCINTILLATIONYIELD2";
  fMatConstPropNames[kELECTRONSCINTILLATIONYIELD3] =  "ELECTRONSCINTILLATIONYIELD3";
  fMatConstPropNames[kCOATEDTHICKNESS] =              "COATEDTHICKNESS";
  fMatConstPropNames[kCOATEDFRUSTRATEDTRANSMISSION] = "COATEDFRUSTRATEDTRANSMISSION";
  fMatConstPropNames[kPROTONSCINTILLATIONTIMECONSTANT1] =   "PROTONSCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kPROTONSCINTILLATIONTIMECONSTANT2] =   "PROTONSCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kPROTONSCINTILLATIONTIMECONSTANT3] =   "PROTONSCINTILLATIONTIMECONSTANT3";
  fMatConstPropNames[kDEUTERONSCINTILLATIONTIMECONSTANT1] = "DEUTERONSCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kDEUTERONSCINTILLATIONTIMECONSTANT2] = "DEUTERONSCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kDEUTERONSCINTILLATIONTIMECONSTANT3] = "DEUTERONSCINTILLATIONTIMECONSTANT3";
  fMatConstPropNames[kTRITONSCINTILLATIONTIMECONSTANT1] =   "TRITONSCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kTRITONSCINTILLATIONTIMECONSTANT2] =   "TRITONSCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kTRITONSCINTILLATIONTIMECONSTANT3] =   "TRITONSCINTILLATIONTIMECONSTANT3";
  fMatConstPropNames[kALPHASCINTILLATIONTIMECONSTANT1] =    "ALPHASCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kALPHASCINTILLATIONTIMECONSTANT2] =    "ALPHASCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kALPHASCINTILLATIONTIMECONSTANT3] =    "ALPHASCINTILLATIONTIMECONSTANT3";
  fMatConstPropNames[kIONSCINTILLATIONTIMECONSTANT1] =      "IONSCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kIONSCINTILLATIONTIMECONSTANT2] =      "IONSCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kIONSCINTILLATIONTIMECONSTANT3] =      "IONSCINTILLATIONTIMECONSTANT3";
  fMatConstPropNames[kELECTRONSCINTILLATIONTIMECONSTANT1] = "ELECTRONSCINTILLATIONTIMECONSTANT1";
  fMatConstPropNames[kELECTRONSCINTILLATIONTIMECONSTANT2] = "ELECTRONSCINTILLATIONTIMECONSTANT2";
  fMatConstPropNames[kELECTRONSCINTILLATIONTIMECONSTANT3] = "ELECTRONSCINTILLATIONTIMECONSTANT3";
  fMCP.assign(kNumberOfConstPropertyIndex, {0., false});
}

G4MaterialPropertiesTable::~G4MaterialPropertiesTable()
{
  for (auto prop : fMP) {
    delete (prop);
  }
}

G4int G4MaterialPropertiesTable::GetConstPropertyIndex(const G4String& key) const
{
  // Returns the constant material property index corresponding to a key

  std::size_t index = std::distance(fMatConstPropNames.cbegin(),
    std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key));
  if (index < fMatConstPropNames.size()) {
    return (G4int)index;
  }

  G4ExceptionDescription ed;
  ed << "Constant Material Property Index for key " << key << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetConstPropertyIndex()", "mat200", FatalException, ed);
  return 0;
}

G4int G4MaterialPropertiesTable::GetPropertyIndex(const G4String& key) const
{
  // Returns the material property index corresponding to a key
  std::size_t index = std::distance(
    fMatPropNames.cbegin(), std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key));
  if (index < fMatPropNames.size()) {
    return (G4int)index;
  }
  G4ExceptionDescription ed;
  ed << "Material Property Index for key " << key << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetPropertyIndex()", "mat201", FatalException, ed);
  return 0;
}

G4double G4MaterialPropertiesTable::GetConstProperty(const G4int index) const
{
  // Returns the constant material property corresponding to an index
  // fatal exception if property not found

  if (index < (G4int)fMCP.size() && fMCP[index].second) {
    return fMCP[index].first;
  }
  G4ExceptionDescription ed;
  ed << "Constant Material Property " << fMatConstPropNames[index] << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetConstProperty()", "mat202", FatalException, ed);
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

  return index >= 0 && index < (G4int)fMCP.size() && fMCP[index].second;
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4String& key) const
{
  // Returns true if a const property 'key' exists
  std::size_t index = std::distance(fMatConstPropNames.cbegin(),
    std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key));
  if (index < fMatConstPropNames.size()) {  // index is type std::size_t so >= 0
    return ConstPropertyExists((G4int)index);
  }
  return false;
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const char* key) const
{
  std::size_t index = std::distance(fMatConstPropNames.cbegin(),
    std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key));
  if (index < fMatConstPropNames.size()) {  // index is type std::size_t so >= 0
    return ConstPropertyExists((G4int)index);
  }
  return false;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(const G4String& key) const
{
  // Returns a Material Property Vector corresponding to a key
  if (std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) != fMatPropNames.cend()) {
    const G4int index = GetPropertyIndex(G4String(key));
    return GetProperty(index);
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(const char* key) const
{
  if (std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) != fMatPropNames.cend()) {
    const G4int index = GetPropertyIndex(G4String(key));
    return GetProperty(index);
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(const G4int index) const
{
  // Returns a Material Property Vector corresponding to an index
  // returns nullptr if the property has not been defined by user
  if (index >= 0 && index < (G4int)fMP.size()) {
    return fMP[index];
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::AddProperty(const G4String& key,
  const std::vector<G4double>& photonEnergies, const std::vector<G4double>& propertyValues,
  G4bool createNewKey, G4bool spline)
{
  if (photonEnergies.size() != propertyValues.size()) {
    G4ExceptionDescription ed;
    ed << "AddProperty error. Number of property values must be equal to the number of\n"
       << "energy values. Property name: " << key;
    G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat204", FatalException, ed);
  }

  if (photonEnergies.size() == 1) {
    G4ExceptionDescription ed;
    ed << "AddProperty warning. A material property vector must have more than one value.\n"
       << "Unless you will later add an entry, this is an error.\n"
       << "Property name: " << key;
    G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat218", JustWarning, ed);
  }

  // G4PhysicsVector assumes energies are in increasing order
  for (std::size_t i = 0; i < photonEnergies.size() - 1; ++i) {
    if (photonEnergies.at(i + 1) < photonEnergies.at(i)) {
      G4ExceptionDescription ed;
      ed << "Energies in material property vector must be in increasing "
         << "order. Key: " << key << " Energy: " << photonEnergies.at(i + 1);
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat215", FatalException, ed);
    }
  }

  // if the key doesn't exist, add it if requested
  if (std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) == fMatPropNames.cend()) {
    if (createNewKey) {
      fMatPropNames.push_back(key);
      fMP.push_back(nullptr);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Attempting to create a new material property vector " << key << " without setting\n"
         << "createNewKey parameter of AddProperty to true.";
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat205", FatalException, ed);
    }
  }

  auto* mpv = new G4MaterialPropertyVector(photonEnergies, propertyValues, spline);
  mpv->SetVerboseLevel(1);
  if (spline) {
    mpv->FillSecondDerivatives();
  }
  G4int index = GetPropertyIndex(key);
  fMP[index] = mpv;

  // if key is RINDEX, we calculate GROUPVEL -
  // contribution from Tao Lin (IHEP, the JUNO experiment)
  if (key == "RINDEX") {
    CalculateGROUPVEL();
  }

  return mpv;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::AddProperty(const char* key,
  G4double* photonEnergies, G4double* propertyValues, G4int numEntries, G4bool createNewKey,
  G4bool spline)
{
  // Provides a way of adding a property to the Material Properties
  // Table given a pair of arrays and a key
  G4String k(key);

  std::vector<G4double> energies(photonEnergies, photonEnergies + numEntries);
  std::vector<G4double> values(propertyValues, propertyValues + numEntries);
  return AddProperty(k, energies, values, createNewKey, spline);
}

void G4MaterialPropertiesTable::AddProperty(
  const G4String& key, G4MaterialPropertyVector* mpv, G4bool createNewKey)
{
  //  Provides a way of adding a property to the Material Properties
  //  Table given an G4MaterialPropertyVector Reference and a key

  // G4PhysicsVector assumes energies are in increasing order
  if (mpv->GetVectorLength() > 1) {
    for (std::size_t i = 0; i < mpv->GetVectorLength() - 1; ++i) {
      if (mpv->Energy(i + 1) < mpv->Energy(i)) {
        G4ExceptionDescription ed;
        ed << "Energies in material property vector must be in increasing "
           << "order. Key: " << key << " Energy: " << mpv->Energy(i + 1);
        G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat216", FatalException, ed);
      }
    }
  }

  if (mpv->GetVectorLength() <= 1) {
    G4ExceptionDescription ed;
    ed << "AddProperty warning. A material property vector must have more than one value.\n"
       << "Unless you will later add an entry, this is an error.\n"
       << "Property name: " << key;
    G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat219", JustWarning, ed);
  }

  // if the key doesn't exist, add it
  if (std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) == fMatPropNames.cend()) {
    if (createNewKey) {
      fMatPropNames.push_back(key);
      fMP.push_back(nullptr);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Attempting to create a new material property key " << key << " without setting\n"
         << "createNewKey parameter of AddProperty to true.";
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat206", FatalException, ed);
    }
  }
  G4int index = GetPropertyIndex(key);
  fMP[index] = mpv;

  // if key is RINDEX, we calculate GROUPVEL -
  // contribution from Tao Lin (IHEP, the JUNO experiment)
  if (key == "RINDEX") {
    CalculateGROUPVEL();
  }
}

void G4MaterialPropertiesTable::AddProperty(
  const char* key, G4MaterialPropertyVector* mpv, G4bool createNewKey)
{
  AddProperty(G4String(key), mpv, createNewKey);
}

void G4MaterialPropertiesTable::AddProperty(const G4String& key, const G4String& mat)
{
  // load a material property vector defined in Geant4 source
  G4MaterialPropertyVector* v = G4OpticalMaterialProperties::GetProperty(key, mat);
  AddProperty(key, v);
}

void G4MaterialPropertiesTable::AddConstProperty(
  const G4String& key, G4double propertyValue, G4bool createNewKey)
{
  // Provides a way of adding a constant property to the Material Properties
  // Table given a key
  if (std::find(fMatConstPropNames.cbegin(), fMatConstPropNames.cend(), key) ==
      fMatConstPropNames.cend())
  {
    if (createNewKey) {
      fMatConstPropNames.push_back(key);
      fMCP.emplace_back(0., true);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Attempting to create a new material constant property key " << key
         << " without setting\n"
         << "createNewKey parameter of AddProperty to true.";
      G4Exception("G4MaterialPropertiesTable::AddProperty()", "mat207", FatalException, ed);
    }
  }
  G4int index = GetConstPropertyIndex(key);

  fMCP[index] = std::pair<G4double, G4bool>{propertyValue, true};
}

void G4MaterialPropertiesTable::AddConstProperty(
  const char* key, G4double propertyValue, G4bool createNewKey)
{
  // Provides a way of adding a constant property to the Material Properties
  // Table given a key
  AddConstProperty(G4String(key), propertyValue, createNewKey);
}

void G4MaterialPropertiesTable::RemoveConstProperty(const G4String& key)
{
  G4int index = GetConstPropertyIndex(key);
  if (index < (G4int)fMCP.size()) {
    fMCP[index] = std::pair<G4double, G4bool>{0., false};
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

void G4MaterialPropertiesTable::RemoveProperty(const char* key) { RemoveProperty(G4String(key)); }

void G4MaterialPropertiesTable::AddEntry(
  const G4String& key, G4double aPhotonEnergy, G4double aPropertyValue)
{
  // Allows to add an entry pair directly to the Material Property Vector
  // given a key.
  if (std::find(fMatPropNames.cbegin(), fMatPropNames.cend(), key) == fMatPropNames.cend()) {
    G4ExceptionDescription ed;
    ed << "Material Property Vector " << key << " not found.";
    G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat214", FatalException, ed);
  }
  G4int index = GetPropertyIndex(key);

  G4MaterialPropertyVector* targetVector = fMP[index];
  if (targetVector != nullptr) {
    // do not allow duplicate energies
    for (std::size_t i = 0; i < targetVector->GetVectorLength(); ++i) {
      if (aPhotonEnergy == targetVector->Energy(i)) {
        G4ExceptionDescription ed;
        ed << "Energy values in material property vector must be unique. "
           << "Key: " << key;
        G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat217", FatalException, ed);
      }
    }

    targetVector->InsertValues(aPhotonEnergy, aPropertyValue);
  }
  else {
    G4ExceptionDescription ed;
    ed << "Material Property Vector " << key << " not found.";
    G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat208", FatalException, ed);
  }
  if (key == "RINDEX") {
    CalculateGROUPVEL();
  }
}

void G4MaterialPropertiesTable::AddEntry(
  const char* key, G4double aPhotonEnergy, G4double aPropertyValue)
{
  AddEntry(G4String(key), aPhotonEnergy, aPropertyValue);
}

void G4MaterialPropertiesTable::DumpTable() const
{
  // material properties
  G4int j = 0;
  for (const auto& prop : fMP) {
    if (prop != nullptr) {
      G4cout << j << ": " << fMatPropNames[j] << G4endl;
      prop->DumpValues();
    }
    ++j;
  }
  // material constant properties
  j = 0;
  for (const auto& cprop : fMCP) {
    if (cprop.second) {
      G4cout << j << ": " << fMatConstPropNames[j] << " " << cprop.first << G4endl;
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
  if (fMP[kGROUPVEL] != nullptr) {
    this->RemoveProperty("GROUPVEL");
  }

  // fetch RINDEX data, give up if unavailable
  G4MaterialPropertyVector* rindex = this->GetProperty(kRINDEX);
  if (rindex == nullptr) {
    return nullptr;
  }

  // RINDEX exists but has no entries, give up
  if (rindex->GetVectorLength() == 0) {
    return nullptr;
  }

  // add GROUPVEL vector
  auto* groupvel = new G4MaterialPropertyVector();
  groupvel->SetVerboseLevel(1);

  // fill GROUPVEL vector using RINDEX values
  // rindex built-in "iterator" was advanced to first entry above
  G4double E0 = rindex->Energy(0);
  G4double n0 = (*rindex)[0];

  if (E0 <= 0.) {
    G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat211", FatalException,
      "Optical Photon Energy <= 0");
  }

  if (rindex->GetVectorLength() >= 2) {
    // good, we have at least two entries in RINDEX
    // get next energy/value pair

    G4double E1 = rindex->Energy(1);
    G4double n1 = (*rindex)[1];

    if (E1 <= 0.) {
      G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat212", FatalException,
        "Optical Photon Energy <= 0");
    }

    G4double vg;

    // add entry at first photon energy
    vg = c_light / (n0 + (n1 - n0) / G4Log(E1 / E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    if ((vg < 0) || (vg > c_light / n0)) {
      vg = c_light / n0;
    }

    groupvel->InsertValues(E0, vg);

    // add entries at midpoints between remaining photon energies
    for (std::size_t i = 2; i < rindex->GetVectorLength(); ++i) {
      vg = c_light / (0.5 * (n0 + n1) + (n1 - n0) / G4Log(E1 / E0));

      // allow only for 'normal dispersion' -> dn/d(logE) > 0
      if ((vg < 0) || (vg > c_light / (0.5 * (n0 + n1)))) {
        vg = c_light / (0.5 * (n0 + n1));
      }
      groupvel->InsertValues(0.5 * (E0 + E1), vg);

      // get next energy/value pair, or exit loop
      E0 = E1;
      n0 = n1;
      E1 = rindex->Energy(i);
      n1 = (*rindex)[i];

      if (E1 <= 0.) {
        G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat213", FatalException,
          "Optical Photon Energy <= 0");
      }
    }

    // add entry at last photon energy
    vg = c_light / (n1 + (n1 - n0) / G4Log(E1 / E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    if ((vg < 0) || (vg > c_light / n1)) {
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
