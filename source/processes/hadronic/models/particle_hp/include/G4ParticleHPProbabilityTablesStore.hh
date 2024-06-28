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
// -------------------------------------------------------------------
//
//      Geant4 header file 
//
//      File name: G4ParticleHPProbabilityTablesStore.hh
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//	         Loic Thulliez (CEA France)
// 
//      Creation date: 4 June 2024
//
//      Description: Class to store all probability tables for different isotopes
//                   and in future also for different temperatures.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//
#ifndef G4ParticleHPProbabilityTablesStore_h
#define G4ParticleHPProbabilityTablesStore_h 1

#include "globals.hh"
#include <map>
#include <vector>
#include <thread>

class G4Material;
class G4Element;
class G4Isotope;
class G4DynamicParticle;
class G4ParticleHPIsoProbabilityTable;


class G4ParticleHPProbabilityTablesStore {
public:
  static G4ParticleHPProbabilityTablesStore * GetInstance();

  void Init();
  void InitURRlimits();

  std::vector< std::map< G4int, G4ParticleHPIsoProbabilityTable* > >* GetProbabilityTables() { return ProbabilityTables; };
  std::vector< std::pair< G4double, G4double > >* GetURRlimits(){ return URRlimits; };
  G4double GetIsoCrossSectionPT( const G4DynamicParticle*, G4int, const G4Isotope*, const G4Element*, const G4Material* );

  std::vector< std::map< std::thread::id, G4double > > random_number_cache;

private:
  static G4ParticleHPProbabilityTablesStore* instance;

  G4ParticleHPProbabilityTablesStore();
  G4ParticleHPProbabilityTablesStore( const G4ParticleHPProbabilityTablesStore& ){};
  ~G4ParticleHPProbabilityTablesStore();

  std::vector< std::vector< G4int > >* Temperatures;
  std::vector< std::map< G4int, G4ParticleHPIsoProbabilityTable* > >* ProbabilityTables;
  std::vector< std::pair< G4double, G4double > >* URRlimits;
  std::vector< std::map< std::thread::id, G4double > > energy_cache;
  G4String dirName;
  G4String filename;
  G4int numIso;
  G4bool usedNjoy;
  G4bool usedCalendf;
};

#endif
