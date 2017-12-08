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
// $Id: G4MaterialPropertiesTable.cc 106997 2017-10-31 10:22:36Z gcosmo $
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

}

////////////
// Methods
////////////

G4MCPindex G4MaterialPropertiesTable::GetConstPropertyIndex(const G4String& key) const
{
  // Returns the constant material property index corresponding to a key

  size_t index = std::distance(G4MaterialConstPropertyName.begin(),
                     std::find(G4MaterialConstPropertyName.begin(), 
                                         G4MaterialConstPropertyName.end(), key));
  if((G4MCPindex)index != kNumberOfConstPropertyIndex) return (G4MCPindex)index;
  G4cout << "key: " << key << G4endl;
  G4Exception("G4MaterialPropertiesTable::GetConstPropertyIndex()","mat206",
              FatalException, "Constant Material Property Index not found.");
  return kNullConstPropertyIndex;
} 

G4MPindex G4MaterialPropertiesTable::GetPropertyIndex(const G4String& key) const
{
  // Returns the material property index corresponding to a key
  size_t index = std::distance(G4MaterialPropertyName.begin(),
                     std::find(G4MaterialPropertyName.begin(), 
                                         G4MaterialPropertyName.end(), key));
  if((G4MPindex)index != kNumberOfPropertyIndex) return (G4MPindex)index;
  G4cout << "key: " << key << G4endl;
  G4Exception("G4MaterialPropertiesTable::GetPropertyIndex()","mat207",
               FatalException, "Material Property Index not found.");
  return kNullPropertyIndex;
} 

G4double G4MaterialPropertiesTable::GetConstProperty(const G4MCPindex index) const
{
  // Returns the constant material property corresponding to an index

  MCPiterator j;
  j = MCP.find(index);
  if ( j != MCP.end() ) return j->second;
  G4cout << "index: " << index << G4endl;
  G4Exception("G4MaterialPropertiesTable::GetConstProperty()","mat202",
              FatalException, "Constant Material Property Index not found.");
  return 0.;
}

G4double G4MaterialPropertiesTable::GetConstProperty(const char *key) const
{
  // Returns the constant material property corresponding to a key
  const G4MCPindex index = GetConstPropertyIndex(G4String(key));
  return GetConstProperty(index);
}

G4bool G4MaterialPropertiesTable::ConstPropertyExists(const char *key) const
{
  // Returns true if a const property 'key' exists
  const G4MCPindex index = GetConstPropertyIndex(G4String(key));

  MCPiterator j;
  j = MCP.find(index);
  if ( j != MCP.end() ) return true;
  return false;
}

G4MaterialPropertyVector*
G4MaterialPropertiesTable::GetProperty(const char *key)
{
  // Returns a Material Property Vector corresponding to a key
  const G4MPindex index = GetPropertyIndex(G4String(key));
  return GetProperty(index);
}

G4MaterialPropertyVector*
G4MaterialPropertiesTable::GetProperty(const G4MPindex index)
{
  // Returns a Material Property Vector corresponding to an index
  MPiterator i;
  i = MP.find(index);
  if ( i != MP.end() ) return i->second;
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
  G4MPindex index = GetPropertyIndex(G4String(key));

  G4MaterialPropertyVector *mpv = new G4MaterialPropertyVector(PhotonEnergies, 
                                                   PropertyValues, NumEntries);
  MP[index] = mpv;

  // if key is RINDEX, we calculate GROUPVEL - 
  // contribution from Tao Lin (IHEP, the JUNO experiment) 
  if (G4String(key)=="RINDEX") {
      CalculateGROUPVEL();
  }

  return mpv;
}

void G4MaterialPropertiesTable::
AddProperty(const char *key, G4MaterialPropertyVector *mpv)
{
  //  Provides a way of adding a property to the Material Properties
  //  Table given an G4MaterialPropertyVector Reference and a key
  G4MPindex index = GetPropertyIndex(G4String(key));
  MP[ index ] = mpv;

  // if key is RINDEX, we calculate GROUPVEL -
  // contribution from Tao Lin (IHEP, the JUNO experiment) 
  if (G4String(key)=="RINDEX") {
      CalculateGROUPVEL();
  }
} 

void G4MaterialPropertiesTable::AddEntry(const char *key,
                                         G4double    aPhotonEnergy,
                                         G4double    aPropertyValue)
{
  // Allows to add an entry pair directly to the Material Property Vector
  // given a key
  const G4MPindex index = GetPropertyIndex(G4String(key));

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
  MPiterator i;
  for (i = MP.begin(); i != MP.end(); ++i)
  {
    G4cout << (*i).first << G4endl;
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
  MCPiterator j;
  for (j = MCP.begin(); j != MCP.end(); ++j)
  {
    G4cout << j->first << G4endl;
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
