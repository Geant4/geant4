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
//      File name: G4ParticleHPIsoProbabilityTable_NJOY.hh
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//	         Loic Thulliez (CEA France)
// 
//      Creation date: 4 June 2024
//
//      Description: Class for the probability table of the given isotope
//                   and for the given temperature generated with NJOY.
//                   It reads the files with probability tables and
//                   finds the correct cross-section.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//
#ifndef G4ParticleHPIsoProbabilityTable_NJOY_h
#define G4ParticleHPIsoProbabilityTable_NJOY_h 1

#include "globals.hh"
#include "G4ParticleHPInterpolator.hh"
#include "G4ParticleHPIsoProbabilityTable.hh"
#include <vector>
#include <thread>
#include <map>

class G4DynamicParticle;
class G4Element;


class G4ParticleHPIsoProbabilityTable_NJOY : public G4ParticleHPIsoProbabilityTable {
public:
  G4ParticleHPIsoProbabilityTable_NJOY();
  ~G4ParticleHPIsoProbabilityTable_NJOY();
  void Init( G4int, G4int, G4int, G4double, G4String ) override;
  G4double GetCorrelatedIsoCrossSectionPT( const G4DynamicParticle*, G4int, const G4Element*, G4double&, G4double&, 
                                           std::thread::id& ) override;
  G4double GetIsoCrossSectionPT( const G4DynamicParticle*, G4int, const G4Element*, G4double&, 
                                 std::map< std::thread::id, G4double >&, std::thread::id& ) override;
private:
  G4int tableOrder;
  G4int lssf_flag;
  G4ParticleHPInterpolator theInt;
};

#endif
