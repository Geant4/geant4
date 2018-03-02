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
// $Id: G4MaterialPropertiesTable.cc 108517 2018-02-16 08:18:55Z gcosmo $
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
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Log.hh"

#include <algorithm>

/////////////////
// Constructors
/////////////////

G4MaterialPropertiesTable::G4MaterialPropertiesTable()
{
  // elements of these 2 vectors must be in same order as
  // the corresponding enums in G4MaterialPropertiesIndex.hh
  G4MaterialPropertyName.push_back(G4String("RINDEX"));
  G4MaterialPropertyName.push_back(G4String("REFLECTIVITY"));
  G4MaterialPropertyName.push_back(G4String("REALRINDEX"));
  G4MaterialPropertyName.push_back(G4String("IMAGINARYRINDEX"));
  G4MaterialPropertyName.push_back(G4String("EFFICIENCY"));
  G4MaterialPropertyName.push_back(G4String("TRANSMITTANCE"));
  G4MaterialPropertyName.push_back(G4String("SPECULARLOBECONSTANT"));
  G4MaterialPropertyName.push_back(G4String("SPECULARSPIKECONSTANT"));
  G4MaterialPropertyName.push_back(G4String("BACKSCATTERCONSTANT"));
  G4MaterialPropertyName.push_back(G4String("GROUPVEL"));
  G4MaterialPropertyName.push_back(G4String("MIEHG"));
  G4MaterialPropertyName.push_back(G4String("RAYLEIGH"));
  G4MaterialPropertyName.push_back(G4String("WLSCOMPONENT"));
  G4MaterialPropertyName.push_back(G4String("WLSABSLENGTH"));
  G4MaterialPropertyName.push_back(G4String("ABSLENGTH"));
  G4MaterialPropertyName.push_back(G4String("FASTCOMPONENT"));
  G4MaterialPropertyName.push_back(G4String("SLOWCOMPONENT"));
  G4MaterialPropertyName.push_back(G4String("PROTONSCINTILLATIONYIELD"));
  G4MaterialPropertyName.push_back(G4String("DEUTERONSCINTILLATIONYIELD"));
  G4MaterialPropertyName.push_back(G4String("TRITONSCINTILLATIONYIELD"));
  G4MaterialPropertyName.push_back(G4String("ALPHASCINTILLATIONYIELD"));
  G4MaterialPropertyName.push_back(G4String("IONSCINTILLATIONYIELD"));
  G4MaterialPropertyName.push_back(G4String("ELECTRONSCINTILLATIONYIELD"));

  G4MaterialConstPropertyName.push_back(G4String("SURFACEROUGHNESS"));
  G4MaterialConstPropertyName.push_back(G4String("ISOTHERMAL_COMPRESSIBILITY"));
  G4MaterialConstPropertyName.push_back(G4String("RS_SCALE_FACTOR"));
  G4MaterialConstPropertyName.push_back(G4String("WLSMEANNUMBERPHOTONS"));
  G4MaterialConstPropertyName.push_back(G4String("WLSTIMECONSTANT"));
  G4MaterialConstPropertyName.push_back(G4String("MIEHG_FORWARD"));
  G4MaterialConstPropertyName.push_back(G4String("MIEHG_BACKWARD"));
  G4MaterialConstPropertyName.push_back(G4String("MIEHG_FORWARD_RATIO"));
  G4MaterialConstPropertyName.push_back(G4String("SCINTILLATIONYIELD"));
  G4MaterialConstPropertyName.push_back(G4String("RESOLUTIONSCALE"));
  G4MaterialConstPropertyName.push_back(G4String("FASTTIMECONSTANT"));
  G4MaterialConstPropertyName.push_back(G4String("FASTSCINTILLATIONRISETIME"));
  G4MaterialConstPropertyName.push_back(G4String("SLOWTIMECONSTANT"));
  G4MaterialConstPropertyName.push_back(G4String("SLOWSCINTILLATIONRISETIME"));
  G4MaterialConstPropertyName.push_back(G4String("YIELDRATIO"));
  G4MaterialConstPropertyName.push_back(G4String("FERMIPOT"));
  G4MaterialConstPropertyName.push_back(G4String("DIFFUSION"));
  G4MaterialConstPropertyName.push_back(G4String("SPINFLIP"));
  G4MaterialConstPropertyName.push_back(G4String("LOSS"));
  G4MaterialConstPropertyName.push_back(G4String("LOSSCS"));
  G4MaterialConstPropertyName.push_back(G4String("ABSCS"));
  G4MaterialConstPropertyName.push_back(G4String("SCATCS"));
  G4MaterialConstPropertyName.push_back(G4String("MR_NBTHETA"));
  G4MaterialConstPropertyName.push_back(G4String("MR_NBE"));
  G4MaterialConstPropertyName.push_back(G4String("MR_RRMS"));
  G4MaterialConstPropertyName.push_back(G4String("MR_CORRLEN"));
  G4MaterialConstPropertyName.push_back(G4String("MR_THETAMIN"));
  G4MaterialConstPropertyName.push_back(G4String("MR_THETAMAX"));
  G4MaterialConstPropertyName.push_back(G4String("MR_EMIN"));
  G4MaterialConstPropertyName.push_back(G4String("MR_EMAX"));
  G4MaterialConstPropertyName.push_back(G4String("MR_ANGNOTHETA"));
  G4MaterialConstPropertyName.push_back(G4String("MR_ANGNOPHI"));
  G4MaterialConstPropertyName.push_back(G4String("MR_ANGCUT"));
}

////////////////
// Destructor
////////////////

G4MaterialPropertiesTable::~G4MaterialPropertiesTable()
{
  MPTiterator i;
  for (i = MPT.begin(); i != MPT.end(); ++i)
  {
    delete (*i).second;
  }
  MPT.clear();
  MPTC.clear();

  MPiterator it;
  for (it = MP.begin(); it != MP.end(); ++it)
  {
    delete (*it).second;
  }
  MP.clear();
  MCP.clear();

  G4MaterialPropertyName.clear();
  G4MaterialConstPropertyName.clear();
}

////////////
// Methods
////////////

G4int G4MaterialPropertiesTable::GetConstPropertyIndex(const G4String& key,
                                                       G4bool warning) const
{
  // Returns the constant material property index corresponding to a key

  size_t index = std::distance(G4MaterialConstPropertyName.begin(),
                     std::find(G4MaterialConstPropertyName.begin(), 
                                     G4MaterialConstPropertyName.end(), key));
  if(index < G4MaterialConstPropertyName.size()) return index;
  if (warning) {
    G4ExceptionDescription ed;
    ed << "Constant Material Property Index for key " << key << " not found.";
    G4Exception("G4MaterialPropertiesTable::GetConstPropertyIndex()","mat206",
                JustWarning, ed);
  }
  return -1;
} 

G4int G4MaterialPropertiesTable::GetPropertyIndex(const G4String& key,
                                                  G4bool warning) const
{
  // Returns the material property index corresponding to a key
  size_t index = std::distance(G4MaterialPropertyName.begin(),
                     std::find(G4MaterialPropertyName.begin(), 
                                         G4MaterialPropertyName.end(), key));
  if(index < G4MaterialPropertyName.size()) return index;
  if (warning) {
    G4ExceptionDescription ed;
    ed << "Material Property Index for key " << key << " not found.";
    G4Exception("G4MaterialPropertiesTable::GetPropertyIndex()","mat207",
                 JustWarning, ed);
  }
  return -1;
} 

G4double G4MaterialPropertiesTable::GetConstProperty(const G4int index) const
{
  // Returns the constant material property corresponding to an index
  // fatal exception if property not found

  MCPiterator j;
  j = MCP.find(index);
  if ( j != MCP.end() ) return j->second;
  G4ExceptionDescription ed;
  ed << "Constant Material Property Index " << index << " not found.";
  G4Exception("G4MaterialPropertiesTable::GetConstProperty()","mat202",
              FatalException, ed);
  return 0.;
}

G4double G4MaterialPropertiesTable::GetConstProperty(const char *key) const
{
  // Returns the constant material property corresponding to a key
  // fatal exception if property not found

  const G4int index = GetConstPropertyIndex(G4String(key));
  return GetConstProperty(index);
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const char *key) const
{
  // Returns true if a const property 'key' exists
  const G4int index = GetConstPropertyIndex(G4String(key));

  MCPiterator j;
  j = MCP.find(index);
  if ( j != MCP.end() ) return true;
  return false;
}

G4MaterialPropertyVector*
G4MaterialPropertiesTable::GetProperty(const char *key, G4bool warning)
{
  // Returns a Material Property Vector corresponding to a key
  const G4int index = GetPropertyIndex(G4String(key), warning);
  return GetProperty(index);
}

G4MaterialPropertyVector*
G4MaterialPropertiesTable::GetProperty(const G4int index, G4bool warning)
{
  // Returns a Material Property Vector corresponding to an index
  MPiterator i;
  i = MP.find(index);
  if ( i != MP.end() ) return i->second;
  if (warning) {
    G4ExceptionDescription ed;
    ed << "Material Property for index " << index << " not found.";
    G4Exception("G4MaterialPropertiesTable::GetPropertyIndex()","mat208",
                 JustWarning, ed);
  }
  return nullptr;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::AddProperty(
                                            const char *key,
                                            G4double   *PhotonEnergies,
                                            G4double   *PropertyValues,
                                            G4int      NumEntries)
{
  // Provides a way of adding a property to the Material Properties
  // Table given a pair of numbers and a key
  G4String k(key);
  // if the key doesn't exist, add it
  if (std::find(G4MaterialPropertyName.begin(),
                G4MaterialPropertyName.end(), k) ==
       G4MaterialPropertyName.end()) {
    G4MaterialPropertyName.push_back(k);
  }
  G4int index = GetPropertyIndex(k);

  G4MaterialPropertyVector *mpv = new G4MaterialPropertyVector(PhotonEnergies, 
                                                   PropertyValues, NumEntries);
  MP[index] = mpv;

  // if key is RINDEX, we calculate GROUPVEL - 
  // contribution from Tao Lin (IHEP, the JUNO experiment) 
  if (k=="RINDEX") {
      CalculateGROUPVEL();
  }

  return mpv;
}

void G4MaterialPropertiesTable::
AddProperty(const char *key, G4MaterialPropertyVector *mpv)
{
  //  Provides a way of adding a property to the Material Properties
  //  Table given an G4MaterialPropertyVector Reference and a key
  G4String k(key);
  // if the key doesn't exist, add it
  if (std::find(G4MaterialPropertyName.begin(),
                G4MaterialPropertyName.end(), k) == 
       G4MaterialPropertyName.end()) {
    G4MaterialPropertyName.push_back(k);
  }
  G4int index = GetPropertyIndex(k);
  MP[ index ] = mpv;

  // if key is RINDEX, we calculate GROUPVEL -
  // contribution from Tao Lin (IHEP, the JUNO experiment) 
  if (k=="RINDEX") {
      CalculateGROUPVEL();
  }
} 

void G4MaterialPropertiesTable::AddEntry(const char *key,
                                         G4double    aPhotonEnergy,
                                         G4double    aPropertyValue)
{
  // Allows to add an entry pair directly to the Material Property Vector
  // given a key
  G4String k(key);
  if (std::find(G4MaterialPropertyName.begin(),
                G4MaterialPropertyName.end(), k) == 
       G4MaterialPropertyName.end()) {
    G4MaterialPropertyName.push_back(k);
  }
  G4int index = GetPropertyIndex(k);

  G4MaterialPropertyVector *targetVector=MP[index];
  if (targetVector != nullptr)
  {
    targetVector->InsertValues(aPhotonEnergy, aPropertyValue);
  }
  else
  {
    G4Exception("G4MaterialPropertiesTable::AddEntry()", "mat203",
                FatalException, "Material Property Vector not found.");
  }
}

void G4MaterialPropertiesTable::DumpTable()
{
  // material properties
  MPiterator i;
  for (i = MP.begin(); i != MP.end(); ++i)
  {
    G4cout << (*i).first << ": "<< G4MaterialPropertyName[(*i).first] <<G4endl;
    if ( (*i).second != 0 )
    {
      (*i).second->DumpValues();
    }
    else
    {
      G4Exception("G4MaterialPropertiesTable::DumpTable()", "mat204",
                  JustWarning, "NULL Material Property Vector Pointer.");
    }
  }
  // material constant properties
  MCPiterator j;
  for (j = MCP.begin(); j != MCP.end(); ++j)
  {
    G4cout << j->first <<": "<< G4MaterialConstPropertyName[j->first] <<G4endl;
    if ( j->second != 0 )
    {
      G4cout << j->second << G4endl;
    }
    else
    {
      G4Exception("G4MaterialPropertiesTable::DumpTable()", "mat202",
                  JustWarning, "No Material Constant Property.");
    }
  }
}

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace {
 G4Mutex materialPropertyTableMutex = G4MUTEX_INITIALIZER;
}
#endif

G4MaterialPropertyVector* G4MaterialPropertiesTable::CalculateGROUPVEL()
{
#ifdef G4MULTITHREADED
  G4AutoLock mptm(&materialPropertyTableMutex);
#endif
  // reconsider (i.e, remove) above mutex as this method is likely called only
  // when RINDEX is added during the detector construction phase (i.e., adding 
  // meterials properties into geometry) on the master thread (Oct. 2017, SYJ)

  // check if "GROUPVEL" already exists
  MPiterator itr;
  itr = MP.find(kGROUPVEL);
  if(itr != MP.end()) return itr->second;

  // fetch RINDEX data, give up if unavailable
  //
  G4MaterialPropertyVector *rindex = this->GetProperty(kRINDEX);
  if (rindex==0)  { return 0; }

  // RINDEX exists but has no entries, give up
  //
  if ( rindex->GetVectorLength() == 0 ) { return 0; }

  // add GROUPVEL vector
  //
  G4MaterialPropertyVector* groupvel = new G4MaterialPropertyVector();

  // fill GROUPVEL vector using RINDEX values
  // rindex built-in "iterator" was advanced to first entry above
  //
  G4double E0 = rindex->Energy(0);
  G4double n0 = (*rindex)[0];

  if (E0 <= 0.)
  {
    G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat205",
                FatalException, "Optical Photon Energy <= 0");
  }
                                                                                
  if ( rindex->GetVectorLength() >= 2 )
  {
    // good, we have at least two entries in RINDEX
    // get next energy/value pair

    G4double E1 = rindex->Energy(1);
    G4double n1 = (*rindex)[1];

    if (E1 <= 0.)
    {
      G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat205",
                  FatalException, "Optical Photon Energy <= 0");
    }

    G4double vg;

    // add entry at first photon energy
    //
    vg = c_light/(n0+(n1-n0)/G4Log(E1/E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    //
    if((vg<0) || (vg>c_light/n0))  { vg = c_light/n0; }

    groupvel->InsertValues( E0, vg );

    // add entries at midpoints between remaining photon energies
    //

    for (size_t i = 2; i < rindex->GetVectorLength(); i++)
    {
      vg = c_light/( 0.5*(n0+n1)+(n1-n0)/G4Log(E1/E0));

      // allow only for 'normal dispersion' -> dn/d(logE) > 0
      //
      if((vg<0) || (vg>c_light/(0.5*(n0+n1))))  { vg = c_light/(0.5*(n0+n1)); }
      groupvel->InsertValues( 0.5*(E0+E1), vg );

      // get next energy/value pair, or exit loop
      //
      E0 = E1;
      n0 = n1;
      E1 = rindex->Energy(i);
      n1 = (*rindex)[i];

      if (E1 <= 0.)
      {
        G4Exception("G4MaterialPropertiesTable::CalculateGROUPVEL()", "mat205",
                    FatalException, "Optical Photon Energy <= 0");
      }
    }

    // add entry at last photon energy
    //
    vg = c_light/(n1+(n1-n0)/G4Log(E1/E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    //
    if((vg<0) || (vg>c_light/n1))  { vg = c_light/n1; }
    groupvel->InsertValues( E1, vg );
  }
  else // only one entry in RINDEX -- weird!
  {
    groupvel->InsertValues( E0, c_light/n0 );
  }
                                                                                
  this->AddProperty( "GROUPVEL", groupvel );
                                                                                
  return groupvel;
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::SetGROUPVEL()
{
  G4String message("SetGROUPVEL will be obsolete from the next release ");
  message += "Use G4MaterialPropertiesTable::CalculateGROUPVEL() instead";

  G4Exception("G4MaterialPropertiesTable::SetGROUPVEL()", "Obsolete",
               JustWarning, message);
  return CalculateGROUPVEL();
}

std::map< G4String, G4MaterialPropertyVector*, std::less<G4String> >*
G4MaterialPropertiesTable::GetPropertiesMap() 
{ 
  // warning message
  G4String message("GetPropertiesMap will be obsolete from the next release ");
  message += "Use G4MaterialPropertiesTable::GetPropertyMap() instead";
  G4Exception("G4MaterialPropertiesTable::GetPropertiesMap()", "Obsolete",
               JustWarning, message);

  for (MPiterator miter = MP.begin(); miter != MP.end(); miter++)
  {
    if(miter->second) {
      MPT [ G4MaterialPropertyName[miter->first] ] = miter->second;
    }
    else {
      G4Exception("G4MaterialPropertiesTable::GetPropertiesMap()","NullPointer",
                   JustWarning, "Null Pointer for Material Property");
      continue;
    }
  }
  return &MPT; 
}

std::map< G4String, G4double, std::less<G4String> >* G4MaterialPropertiesTable::GetPropertiesCMap()
{ 
  // warning message
  G4String message("GetPropertiesCMap will be obsolete from the next release ");
  message += "Use G4MaterialPropertiesTable::GetConstPropertyMap() instead";
  G4Exception("G4MaterialPropertiesTable::GetPropertiesCMap()", "Obsolete",
               JustWarning, message);

  for (MCPiterator miter = MCP.begin(); miter != MCP.end(); miter++) {
    MPTC[ G4MaterialConstPropertyName[miter->first] ] = miter->second;
  }
  return &MPTC; 
}

std::vector<G4String> G4MaterialPropertiesTable::GetMaterialPropertyNames() const
{
  return G4MaterialPropertyName;;
}

std::vector<G4String> G4MaterialPropertiesTable::GetMaterialConstPropertyNames() const
{
  return G4MaterialConstPropertyName;
}
