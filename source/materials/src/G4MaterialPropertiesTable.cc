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
// $Id: G4MaterialPropertiesTable.cc 80573 2014-04-29 07:40:08Z gcosmo $
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
}

////////////
// Methods
////////////

void G4MaterialPropertiesTable::DumpTable()
{
  MPTiterator i;
  for (i = MPT.begin(); i != MPT.end(); ++i)
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
  MPTCiterator j;
  for (j = MPTC.begin(); j != MPTC.end(); ++j)
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

G4MaterialPropertyVector* G4MaterialPropertiesTable::SetGROUPVEL()
{
#ifdef G4MULTITHREADED
  G4AutoLock mptm(&materialPropertyTableMutex);
#endif

  // check if "GROUPVEL" already exists
  MPTiterator itr;
  itr = MPT.find("GROUPVEL");
  if(itr != MPT.end()) return itr->second;

  // fetch RINDEX data, give up if unavailable
  //
  G4MaterialPropertyVector *rindex = this->GetProperty("RINDEX");
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
    G4Exception("G4MaterialPropertiesTable::SetGROUPVEL()", "mat205",
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
      G4Exception("G4MaterialPropertiesTable::SetGROUPVEL()", "mat205",
                  FatalException, "Optical Photon Energy <= 0");
    }

    G4double vg;

    // add entry at first photon energy
    //
    vg = c_light/(n0+(n1-n0)/std::log(E1/E0));

    // allow only for 'normal dispersion' -> dn/d(logE) > 0
    //
    if((vg<0) || (vg>c_light/n0))  { vg = c_light/n0; }

    groupvel->InsertValues( E0, vg );

    // add entries at midpoints between remaining photon energies
    //

    for (size_t i = 2; i < rindex->GetVectorLength(); i++)
    {
      vg = c_light/( 0.5*(n0+n1)+(n1-n0)/std::log(E1/E0));

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
        G4Exception("G4MaterialPropertiesTable::SetGROUPVEL()", "mat205",
                    FatalException, "Optical Photon Energy <= 0");
      }
    }

    // add entry at last photon energy
    //
    vg = c_light/(n1+(n1-n0)/std::log(E1/E0));

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
