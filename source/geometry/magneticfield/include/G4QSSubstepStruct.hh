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
// Structs used by G4QSStepper.

// Author: Mattias Portnoy (Univ. Buenos Aires), 2024
// --------------------------------------------------------------------
#ifndef G4QSS_SUBSTEPSTRUCT_HH
#define G4QSS_SUBSTEPSTRUCT_HH

#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4qss_misc.hh"

#include <map>
#include <cmath>

constexpr G4int MAX_QSS_ORDER=3;
constexpr G4int NUMBER_OF_VARIABLES_QSS = 6;

typedef G4double QSStateVector[6];

/**
 * @brief Structs used by G4QSStepper.
 */

struct Substep
{
  QSStateVector state_x[MAX_QSS_ORDER+1];
  QSStateVector state_q[MAX_QSS_ORDER];
  QSStateVector state_tx;
  QSStateVector state_tq;
  QSStateVector sync_t;
  G4double t{0.0};
  // simple id method so that substeps can have different orders in same step
  G4int extrapolation_method =0;  // To be overwritten in run (expect: 1,2 or 3)
  G4double b_field[3] = {0.0, 0.0, 0.0 };

  inline Substep();
};

// --------------------------------------------------------------------------------

struct Substeps
{
  G4int _arrlength = 30;
  Substep* _substeps = nullptr;
  G4int current_substep_index = -1;

  Substeps();
  ~Substeps();

  /** Mimics the functionality of GNU method reallocarray. */
  void* safe_reallocarray(void* ptr, size_t numMembers, size_t size);
   
  inline void resize();

  inline Substep* create_substep();

  inline void save_substep(Substep* substep);

  inline void reset();

};

#include "G4QSSubstepStruct.icc"

#endif
