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
// Structs used by G4QSStepper
//
// Author: Mattias Portnoy (Univ. Buenos Aires) - 2024
// --------------------------------------------------------------------
#ifndef G4QSS_SUBSTEPSTRUCT_HH
#define G4QSS_SUBSTEPSTRUCT_HH 1

#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4qss_misc.hh"

#include <map>
#include <cmath>

constexpr G4int MAX_QSS_ORDER=3;
typedef G4double QSStateVector[6];

struct Substep
{
  QSStateVector state_x[MAX_QSS_ORDER+1];
  QSStateVector state_q[MAX_QSS_ORDER];
  QSStateVector state_tx;
  QSStateVector state_tq;
  QSStateVector sync_t;
  G4double t;
  // simple id method so that substeps can have different orders in same step
  G4int extrapolation_method;
  G4double b_field[3];
};

struct Substeps
{
  G4int _arrlength = 30;
  Substep* _substeps = static_cast<Substep *>(malloc((_arrlength) * sizeof(Substep)));
  G4int current_substep_index = -1;

  // Mimics the functionality of GNU method reallocarray 
  void* safe_reallocarray(void* ptr, size_t numMembers, size_t size)
  {
    if (size != 0 && numMembers > std::numeric_limits<size_t>::max() / size)
    {
      return nullptr;
    }
    return realloc(ptr, numMembers * size);
  }
   
  inline void resize()
  {
    _arrlength = fmax(_arrlength*2, 1500);
    _substeps = static_cast<Substep *>(safe_reallocarray(_substeps, _arrlength, sizeof(Substep)));
    if( _substeps == nullptr )
    { 
       G4ExceptionDescription ermsg;
       ermsg << "QSS2: Size of state exceed available memory : number of elemets = " << _arrlength
             << " size of each element= " << sizeof(Substep) << G4endl;
       G4Exception( "G4QSSubstepStruct::resize", "GeomField0008", FatalException, ermsg ); 
    }
  }

  inline Substep* create_susbtep()
  {
    current_substep_index++;

    if (unlikely( current_substep_index >= _arrlength ))
    {
      resize();
    }
    return &(_substeps[current_substep_index]);
  }
  inline void save_substep(Substep* substep)
  {
    memcpy(create_susbtep(), substep, sizeof(Substep));
  }

  inline void reset()
  {
    current_substep_index = -1;
  }

};

#endif
