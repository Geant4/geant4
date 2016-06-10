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
// $Id: G4TouchableHistory.hh 86527 2014-11-13 15:06:24Z gcosmo $
//
// 
// class G4TouchableHistory
//
// Class description:
//
// Object representing a touchable detector element, and its history in the
// geometrical hierarchy, including its net resultant local->global transform.

// History:
// - Created. Paul Kent, August 1996
// ----------------------------------------------------------------------
#ifndef G4TOUCHABLEHISTORY_HH
#define G4TOUCHABLEHISTORY_HH

#include "G4VTouchable.hh"

#include "G4NavigationHistory.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "geomwdefs.hh"

class G4TouchableHistory : public G4VTouchable
{

 public:  // with description

  G4TouchableHistory(); 
    // The default constructor produces a touchable-history of 
    // 'zero-depth', ie an "unphysical" and not very unusable one.
    // It is for initialisation only.  

  G4TouchableHistory( const G4NavigationHistory& history );
    // Copy constructor

  ~G4TouchableHistory();
    // Destructor

  inline G4VPhysicalVolume* GetVolume( G4int depth=0 ) const;
  inline G4VSolid* GetSolid( G4int depth=0 ) const;
  const G4ThreeVector& GetTranslation( G4int depth=0 ) const;
  const G4RotationMatrix* GetRotation( G4int depth=0 ) const;

  inline G4int GetReplicaNumber( G4int depth=0 ) const;
  inline G4int GetHistoryDepth()  const;
  G4int MoveUpHistory( G4int num_levels = 1 );
    // Access methods for touchables with history

  void  UpdateYourself( G4VPhysicalVolume*   pPhysVol,
                        const G4NavigationHistory* history=0 ); 
    // Update methods for touchables with history

 public:  // without description

  inline const G4NavigationHistory* GetHistory() const;
    // Should this method be "deprecated" ?
    // it is used now in G4Navigator::LocateGlobalPointAndSetup

  inline void *operator new(size_t);
  inline void operator delete(void *aTH);
    // Override "new" and "delete" to use "G4Allocator".

 private:

  inline G4int CalculateHistoryIndex( G4int stackDepth ) const;

  G4RotationMatrix frot;
  G4ThreeVector ftlate;
  G4NavigationHistory fhistory;
};

#include "G4TouchableHistory.icc"
#endif
