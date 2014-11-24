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
// $Id: G4CascadeChannelTables.cc 69632 2013-05-09 01:17:48Z dwright $
//
// Factory to return pointer to Bertini cross-section table based on
// collision initial state (hadron type codes).
//
// Author:  Michael Kelsey (SLAC)
//
// 20110729  M. Kelsey -- Use static instance() function to work around
//		"disappearance" bug on Linux (GCC 4.1.2).  Add diagnostics.
// 20110916  M. Kelsey -- Add "load on demand" to GetTable(), with full set
//		of channel .hh files for use with LoadTable().
// 20110923  M. Kelsey -- Add optional stream& argument to printTable(),
//		pass through.
// 20111007  M. Kelsey -- Add new gamma-n and gamma-p tables.
// 20130129  M. Kelsey -- Drop load-on-demand interfaces, fill in ctor
// 20130429  M. Kelsey -- Change instance to thread-local pointer.
// 20141121  Use G4AutoDelete to avoid end-of-thread memory leaks

#include "G4CascadeChannelTables.hh"
#include "G4AutoDelete.hh"
#include "G4CascadeChannel.hh"
#include "G4CascadeGamNChannel.hh"
#include "G4CascadeGamPChannel.hh"
#include "G4CascadeKminusNChannel.hh"
#include "G4CascadeKminusPChannel.hh"
#include "G4CascadeKplusNChannel.hh"
#include "G4CascadeKplusPChannel.hh"
#include "G4CascadeKzeroBarNChannel.hh"
#include "G4CascadeKzeroBarPChannel.hh"
#include "G4CascadeKzeroNChannel.hh"
#include "G4CascadeKzeroPChannel.hh"
#include "G4CascadeLambdaNChannel.hh"
#include "G4CascadeLambdaPChannel.hh"
#include "G4CascadeNNChannel.hh"
#include "G4CascadeNPChannel.hh"
#include "G4CascadePPChannel.hh"
#include "G4CascadePiMinusNChannel.hh"
#include "G4CascadePiMinusPChannel.hh"
#include "G4CascadePiPlusNChannel.hh"
#include "G4CascadePiPlusPChannel.hh"
#include "G4CascadePiZeroNChannel.hh"
#include "G4CascadePiZeroPChannel.hh"
#include "G4CascadeSigmaMinusNChannel.hh"
#include "G4CascadeSigmaMinusPChannel.hh"
#include "G4CascadeSigmaPlusNChannel.hh"
#include "G4CascadeSigmaPlusPChannel.hh"
#include "G4CascadeSigmaZeroNChannel.hh"
#include "G4CascadeSigmaZeroPChannel.hh"
#include "G4CascadeXiMinusNChannel.hh"
#include "G4CascadeXiMinusPChannel.hh"
#include "G4CascadeXiZeroNChannel.hh"
#include "G4CascadeXiZeroPChannel.hh"
#include "G4CascadeOmegaMinusNChannel.hh"
#include "G4CascadeOmegaMinusPChannel.hh"
#include "G4CascadeMuMinusPChannel.hh"
#include "G4InuclParticleNames.hh"
#include <iostream>
#include <map>
using namespace G4InuclParticleNames;


// Singleton is created at first invocation

G4ThreadLocal G4CascadeChannelTables* G4CascadeChannelTables::theInstance = 0;

const G4CascadeChannelTables& G4CascadeChannelTables::instance() {
  if (!theInstance) {
    theInstance = new G4CascadeChannelTables;
    G4AutoDelete::Register(theInstance);
  }

  return *theInstance;
}


// Constructor and destructor fully populate tables

G4CascadeChannelTables::G4CascadeChannelTables() {
  tables.clear();
  tables[gam*neu] = new G4CascadeGamNChannel;
  tables[gam*pro] = new G4CascadeGamPChannel;
  tables[k0*neu]  = new G4CascadeKzeroNChannel;
  tables[k0*pro]  = new G4CascadeKzeroPChannel;
  tables[k0b*neu] = new G4CascadeKzeroBarNChannel;
  tables[k0b*pro] = new G4CascadeKzeroBarPChannel;
  tables[kmi*neu] = new G4CascadeKminusNChannel;
  tables[kmi*pro] = new G4CascadeKminusPChannel;
  tables[kpl*neu] = new G4CascadeKplusNChannel;
  tables[kpl*pro] = new G4CascadeKplusPChannel;
  tables[lam*neu] = new G4CascadeLambdaNChannel;
  tables[lam*pro] = new G4CascadeLambdaPChannel;
  tables[neu*neu] = new G4CascadeNNChannel;
  tables[neu*pro] = new G4CascadeNPChannel;
  tables[pi0*neu] = new G4CascadePiZeroNChannel;
  tables[pi0*pro] = new G4CascadePiZeroPChannel;
  tables[pim*neu] = new G4CascadePiMinusNChannel;
  tables[pim*pro] = new G4CascadePiMinusPChannel;
  tables[pip*neu] = new G4CascadePiPlusNChannel;
  tables[pip*pro] = new G4CascadePiPlusPChannel;
  tables[pro*pro] = new G4CascadePPChannel;
  tables[s0*neu]  = new G4CascadeSigmaZeroNChannel;
  tables[s0*pro]  = new G4CascadeSigmaZeroPChannel;
  tables[sm*neu]  = new G4CascadeSigmaMinusNChannel;
  tables[sm*pro]  = new G4CascadeSigmaMinusPChannel;
  tables[sp*neu]  = new G4CascadeSigmaPlusNChannel;
  tables[sp*pro]  = new G4CascadeSigmaPlusPChannel;
  tables[xi0*neu] = new G4CascadeXiZeroNChannel;
  tables[xi0*pro] = new G4CascadeXiZeroPChannel;
  tables[xim*neu] = new G4CascadeXiMinusNChannel;
  tables[xim*pro] = new G4CascadeXiMinusPChannel;
  tables[om*neu]  = new G4CascadeOmegaMinusNChannel;
  tables[om*pro]  = new G4CascadeOmegaMinusPChannel;
  tables[mum*pro] = new G4CascadeMuMinusPChannel;
}

G4CascadeChannelTables::~G4CascadeChannelTables() {
  TableMap::iterator entry;
  for (entry = tables.begin(); entry != tables.end(); ++entry) {
    delete entry->second; entry->second = 0;
  }

  tables.clear();
}


// Argument is interaction code, product of G4InuclEP types

const G4CascadeChannel* G4CascadeChannelTables::GetTable(G4int initialState) {
  return instance().FindTable(initialState);
}

// Arguments are individual G4InuclElementaryParticle types

const G4CascadeChannel* 
G4CascadeChannelTables::GetTable(G4int had1, G4int had2) {
  return GetTable(had1*had2);
}

// Return cross-section table requested by user

const G4CascadeChannel* 
G4CascadeChannelTables::FindTable(G4int initialState) const {
#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << "G4CascadeChannelTables::FindTable " << initialState << G4endl;
#endif
  TableMap::const_iterator entry = tables.find(initialState);
  return (entry != tables.end()) ? entry->second : 0;
}


// Convenience functions for diagnostic output

void G4CascadeChannelTables::Print(std::ostream& os) {
  const TableMap& theTables = instance().tables;	// For convenience
  TableMap::const_iterator entry;
  for (entry = theTables.begin(); entry != theTables.end(); ++entry) {
    if (entry->second) entry->second->printTable(os);
  }
}

void G4CascadeChannelTables::PrintTable(G4int initialState, std::ostream& os) {
  const G4CascadeChannel* tbl = GetTable(initialState);
  if (tbl) tbl->printTable(os);
}
