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
// Convert particles FLUKA <-> G4 worlds.
//
// NB 1: Conversion tables are created at initialization time, 
// and only accessed at run time.
//
// NB 2: At run time, if ever a particle returned from the FLUKA interface 
// is not found in the G4 particle table, an exception is thrown.
//
// NB 3: There are FLUKA particles categories,
// which are in the FLUKA particles tables, and which do not make sense to convert here,
// since they should not be retuned from the FLUKA interface as final state secondaries anyway!
//
// Author: V.Vlachoudis, G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef FLUKA_PARTICLE_TABLE_HH
#define FLUKA_PARTICLE_TABLE_HH


#include <map>
#include <unordered_map>

// G4
#include "globals.hh"


class G4ParticleDefinition;


namespace fluka_particle_table {

  void initialize();

  // FLUKA
  const G4String& fluka2name(const G4int id);

  // FLUKA -> G4
  const G4ParticleDefinition* fluka2geant(const G4int ij);
  const G4String& fluka2geantName(const G4String& name);
  const G4ParticleDefinition* fluka2geant(const G4String& name);

  // G4 -> FLUKA
  G4int geant2fluka(const G4ParticleDefinition* def);
  const G4String& geant2flukaName(const G4ParticleDefinition* def);
}


#endif
#endif // G4_USE_FLUKA
