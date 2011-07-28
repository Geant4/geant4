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

#include "G4CascadeChannelTables.hh"

#include "G4CascadeChannel.hh"
#include "G4CascadePiMinusNChannel.hh"
#include "G4CascadePiMinusPChannel.hh"
#include "G4CascadePiPlusNChannel.hh"
#include "G4CascadePiPlusPChannel.hh"
#include "G4CascadePiZeroNChannel.hh"
#include "G4CascadePiZeroPChannel.hh"
#include "G4CascadeKplusPChannel.hh"
#include "G4CascadeKplusNChannel.hh"
#include "G4CascadeKzeroPChannel.hh"
#include "G4CascadeKzeroNChannel.hh"
#include "G4CascadeKminusPChannel.hh"
#include "G4CascadeKminusNChannel.hh"
#include "G4CascadeKzeroBarPChannel.hh"
#include "G4CascadeKzeroBarNChannel.hh"
#include "G4CascadeNNChannel.hh"
#include "G4CascadeNPChannel.hh"
#include "G4CascadePPChannel.hh"
#include "G4CascadeLambdaPChannel.hh"
#include "G4CascadeLambdaNChannel.hh"
#include "G4CascadeSigmaPlusPChannel.hh"
#include "G4CascadeSigmaPlusNChannel.hh"
#include "G4CascadeSigmaZeroPChannel.hh"
#include "G4CascadeSigmaZeroNChannel.hh"
#include "G4CascadeSigmaMinusPChannel.hh"
#include "G4CascadeSigmaMinusNChannel.hh"
#include "G4CascadeXiZeroPChannel.hh"
#include "G4CascadeXiZeroNChannel.hh"
#include "G4CascadeXiMinusPChannel.hh"
#include "G4CascadeXiMinusNChannel.hh"

#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

#include <map>


// Return cross-section table requested by user, or create it

const G4CascadeChannel* 
G4CascadeChannelTables::GetTable(G4int initialState) {
  static std::map<G4int, G4CascadeChannel*> tables;

  if (tables.find(initialState) == tables.end()) {	// Create only if needed
    switch (initialState) {
    case pro*pro: tables[initialState] = new G4CascadePPChannel; break;
    case pro*neu: tables[initialState] = new G4CascadeNPChannel; break;
    case pip*pro: tables[initialState] = new G4CascadePiPlusPChannel; break;
    case neu*neu: tables[initialState] = new G4CascadeNNChannel; break;
    case pim*pro: tables[initialState] = new G4CascadePiMinusPChannel; break;
    case pip*neu: tables[initialState] = new G4CascadePiPlusNChannel; break;
    case pi0*pro: tables[initialState] = new G4CascadePiZeroPChannel; break;
    case pim*neu: tables[initialState] = new G4CascadePiMinusNChannel; break;
    case kpl*pro: tables[initialState] = new G4CascadeKplusPChannel; break;
    case kmi*pro: tables[initialState] = new G4CascadeKminusPChannel; break;
    case pi0*neu: tables[initialState] = new G4CascadePiZeroNChannel; break;
    case k0*pro:  tables[initialState] = new G4CascadeKzeroPChannel; break;
    case k0b*pro: tables[initialState] = new G4CascadeKzeroBarPChannel; break;
    case lam*pro: tables[initialState] = new G4CascadeLambdaPChannel; break;
    case kpl*neu: tables[initialState] = new G4CascadeKplusNChannel; break;
    case sp*pro:  tables[initialState] = new G4CascadeSigmaPlusPChannel; break;
    case s0*pro:  tables[initialState] = new G4CascadeSigmaZeroPChannel; break;
    case kmi*neu: tables[initialState] = new G4CascadeKminusNChannel; break;
    case sm*pro:  tables[initialState] = new G4CascadeSigmaMinusPChannel; break;
    case xi0*pro: tables[initialState] = new G4CascadeXiZeroPChannel; break;
    case k0*neu:  tables[initialState] = new G4CascadeKzeroNChannel; break;
    case xim*pro: tables[initialState] = new G4CascadeXiMinusPChannel; break;
    case k0b*neu: tables[initialState] = new G4CascadeKzeroBarNChannel; break;
    case lam*neu: tables[initialState] = new G4CascadeLambdaNChannel; break;
    case sp*neu:  tables[initialState] = new G4CascadeSigmaPlusNChannel; break;
    case s0*neu:  tables[initialState] = new G4CascadeSigmaZeroNChannel; break;
    case sm*neu:  tables[initialState] = new G4CascadeSigmaMinusNChannel; break;
    case xi0*neu: tables[initialState] = new G4CascadeXiZeroNChannel; break;
    case xim*neu: tables[initialState] = new G4CascadeXiMinusNChannel; break;
    default:      tables[initialState] = 0;
    }
  }

  return tables[initialState];
}


// Convenience function for diagnostic output

void G4CascadeChannelTables::PrintTable(G4int initialState) {
  const G4CascadeChannel* tbl = GetTable(initialState);
  if (tbl) tbl->printTable();
}
