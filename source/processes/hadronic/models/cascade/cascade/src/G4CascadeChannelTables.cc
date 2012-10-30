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

#include "G4CascadeChannelTables.hh"
#include "G4CascadeChannel.hh"
#include <iostream>
#include <map>


// Singleton is created at first invocation

G4CascadeChannelTables& G4CascadeChannelTables::instance() {
  static G4CascadeChannelTables theInstance;
  return theInstance;
}


// Argument is interaction code, product of G4InuclEP types

const G4CascadeChannel* G4CascadeChannelTables::GetTable(G4int initialState) {
  const G4CascadeChannel* theTable = instance().FindTable(initialState);
  if (!theTable) theTable = instance().LoadTable(initialState);
  return theTable;
}

// Arguments are individual G4InuclElementaryParticle types

const G4CascadeChannel* 
G4CascadeChannelTables::GetTable(G4int had1, G4int had2) {
  return GetTable(had1*had2);
}

// Register cross-section table for later lookup

void  
G4CascadeChannelTables::AddTable(G4int initialState, G4CascadeChannel* table) {
  instance().SaveTable(initialState, table);
}

// Return cross-section table requested by user

const G4CascadeChannel* 
G4CascadeChannelTables::FindTable(G4int initialState) {
#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << "G4CascadeChannelTables::FindTable " << initialState << G4endl;
#endif
  return (tables.find(initialState)!=tables.end()) ? tables[initialState] : 0;
}


// Register specified table in list, replacing previous version

void 
G4CascadeChannelTables::SaveTable(G4int initialState, G4CascadeChannel* table) {
#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << "G4CascadeChannelTables::SaveTable " << initialState << G4endl;
#endif
  if (!table) return;		// Avoid unnecessary work

  if (FindTable(initialState)) delete tables[initialState];
  tables[initialState] = table;
}


// Convenience function for diagnostic output

void G4CascadeChannelTables::PrintTable(G4int initialState, std::ostream& os) {
  const G4CascadeChannel* tbl = GetTable(initialState);
  if (tbl) tbl->printTable(os);
}


// Special function to create and store table for specified interaction

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
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

const G4CascadeChannel* G4CascadeChannelTables::LoadTable(G4int initialState) {
#ifdef G4CASCADE_DEBUG_SAMPLER
  G4cout << "G4CascadeChannelTables::LoadTable " << initialState << G4endl;
#endif

  G4CascadeChannel* tbl = 0;
  switch (initialState) {
  case gam*neu: tbl = new G4CascadeGamNChannel; break;
  case gam*pro: tbl = new G4CascadeGamPChannel; break;
  case k0*neu:  tbl = new G4CascadeKzeroNChannel; break;
  case k0*pro:  tbl = new G4CascadeKzeroPChannel; break;
  case k0b*neu: tbl = new G4CascadeKzeroBarNChannel; break;
  case k0b*pro: tbl = new G4CascadeKzeroBarPChannel; break;
  case kmi*neu: tbl = new G4CascadeKminusNChannel; break;
  case kmi*pro: tbl = new G4CascadeKminusPChannel; break;
  case kpl*neu: tbl = new G4CascadeKplusNChannel; break;
  case kpl*pro: tbl = new G4CascadeKplusPChannel; break;
  case lam*neu: tbl = new G4CascadeLambdaNChannel; break;
  case lam*pro: tbl = new G4CascadeLambdaPChannel; break;
  case neu*neu: tbl = new G4CascadeNNChannel; break;
  case neu*pro: tbl = new G4CascadeNPChannel; break;
  case pi0*neu: tbl = new G4CascadePiZeroNChannel; break;
  case pi0*pro: tbl = new G4CascadePiZeroPChannel; break;
  case pim*neu: tbl = new G4CascadePiMinusNChannel; break;
  case pim*pro: tbl = new G4CascadePiMinusPChannel; break;
  case pip*neu: tbl = new G4CascadePiPlusNChannel; break;
  case pip*pro: tbl = new G4CascadePiPlusPChannel; break;
  case pro*pro: tbl = new G4CascadePPChannel; break;
  case s0*neu:  tbl = new G4CascadeSigmaZeroNChannel; break;
  case s0*pro:  tbl = new G4CascadeSigmaZeroPChannel; break;
  case sm*neu:  tbl = new G4CascadeSigmaMinusNChannel; break;
  case sm*pro:  tbl = new G4CascadeSigmaMinusPChannel; break;
  case sp*neu:  tbl = new G4CascadeSigmaPlusNChannel; break;
  case sp*pro:  tbl = new G4CascadeSigmaPlusPChannel; break;
  case xi0*neu: tbl = new G4CascadeXiZeroNChannel; break;
  case xi0*pro: tbl = new G4CascadeXiZeroPChannel; break;
  case xim*neu: tbl = new G4CascadeXiMinusNChannel; break;
  case xim*pro: tbl = new G4CascadeXiMinusPChannel; break;
  case om*neu:  tbl = new G4CascadeOmegaMinusNChannel; break;
  case om*pro:  tbl = new G4CascadeOmegaMinusPChannel; break;
  default: tbl = 0;
  }

  SaveTable(initialState, tbl);
  return tbl;
}
