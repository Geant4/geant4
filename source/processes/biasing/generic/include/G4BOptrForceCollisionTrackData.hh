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
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//    Extends G4Track properties with information needed for the
//    force collision biasing operator.
//    The G4BOptrForceCollision class is made friend of this one in
//    order to keep unexposed to public most of the data members, as
//    they are used to control the logic.
//
// ------------------ G4BOptrForceCollisionTrackData ------------------
//
// Author: M.Verderi (LLR), October 2015
//
// --------------------------------------------------------------------

#ifndef G4BOptrForceCollisionTrackData_hh
#define G4BOptrForceCollisionTrackData_hh

class G4BOptrForceCollision;
#include "G4VAuxiliaryTrackInformation.hh"

enum class ForceCollisionState { free, toBeCloned, toBeForced, toBeFreeFlight };

class G4BOptrForceCollisionTrackData : public G4VAuxiliaryTrackInformation {

friend class G4BOptrForceCollision;
  
public:
  G4BOptrForceCollisionTrackData( const G4BOptrForceCollision* );
  ~G4BOptrForceCollisionTrackData();
  
  // -- from base class:
  void Print() const;

  // -- Get methods:
  G4bool                             IsFreeFromBiasing() const
  { return ( fForceCollisionState == ForceCollisionState::free);}
  // -- no set methods are provided : sets are made under exclusive control of G4BOptrForceCollision objects through friendness.
  
private:
  const G4BOptrForceCollision* fForceCollisionOperator;
  ForceCollisionState             fForceCollisionState;

  void Reset()
  {
    fForceCollisionOperator = nullptr;
    fForceCollisionState    = ForceCollisionState::free;
  }
  
};

#endif
