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
//      File name: G4ParticleHPIsoProbabilityTable.hh
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//	         Loic Thulliez (CEA France)
// 
//      Creation date: 4 June 2024
//
//      Description: Class for the probability table of the given isotope
//                   and for the given temperature.
//                   It reads the files with probability tables and
//                   finds the correct cross-section.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//
#ifndef G4ParticleHPIsoProbabilityTable_h
#define G4ParticleHPIsoProbabilityTable_h 1

#include "globals.hh"

#include <vector>
#include <thread>
#include <map>

class G4ParticleHPVector;
class G4DynamicParticle;
class G4Element;

class G4ParticleHPIsoProbabilityTable
{
 public:

  G4ParticleHPIsoProbabilityTable() = default;
  virtual ~G4ParticleHPIsoProbabilityTable();
  virtual void Init( G4int, G4int, G4int, G4double, const G4String& );
  virtual G4double GetCorrelatedIsoCrossSectionPT( const G4DynamicParticle*, G4int, const G4Element*, 
                                                   G4double&, G4double&, std::thread::id& );
  virtual G4double GetIsoCrossSectionPT( const G4DynamicParticle*, G4int, const G4Element*, G4double&, 
                                         std::map< std::thread::id, G4double >&, std::thread::id& );

 protected:

  G4double GetDopplerBroadenedElasticXS( const G4DynamicParticle*, G4int, G4int );
  G4double GetDopplerBroadenedCaptureXS( const G4DynamicParticle*, G4int, G4int );
  G4double GetDopplerBroadenedFissionXS( const G4DynamicParticle*, G4int, G4int );
  G4double GetDopplerBroadenedInelasticXS( const G4DynamicParticle*, G4int, G4int );

  G4int Z = 0;
  G4int A = 0;
  G4int m = -1;
  G4double T = -1.;

  G4double Emin = DBL_MAX;
  G4double Emax = 0.;
  G4int nEnergies = 0;

  std::map< std::thread::id, G4double > energy_cache;
  std::map< std::thread::id, G4double > xsela_cache;
  std::map< std::thread::id, G4double > xscap_cache;
  std::map< std::thread::id, G4double > xsfiss_cache;

  G4ParticleHPVector* theEnergies = nullptr;
  std::vector< std::vector< G4double >* >* theProbabilities = nullptr;
  std::vector< std::vector< G4double >* >* theElasticData = nullptr;
  std::vector< std::vector< G4double >* >* theCaptureData = nullptr;
  std::vector< std::vector< G4double >* >* theFissionData = nullptr;
  std::vector< std::vector< G4double >* >* theInelasticData = nullptr;

  G4String filename;
};

#endif
