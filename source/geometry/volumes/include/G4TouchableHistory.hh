// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TouchableHistory.hh,v 1.4 2000-11-01 16:51:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4TouchableHistory
//
// Class description:
//
// Object representing a touchable detector element, and its history in the
// geomtrical hierarchy, including its net resultant local->global transform.

// History:
// - Created. Paul Kent, August 1996

#ifndef G4TOUCHABLEHISTORY_HH
#define G4TOUCHABLEHISTORY_HH

#include "G4VTouchable.hh"

#include "G4NavigationHistory.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4TouchableHistory : public G4VTouchable
{

 public:  // with description

  G4TouchableHistory(const G4NavigationHistory& history);
  G4TouchableHistory(); 
    // The default constructor produces a touchable-history of 
    // 'zero-depth', ie an "unphysical" and not very unusable one.
    // It is for initialisation only.  

  virtual ~G4TouchableHistory();

  inline G4VPhysicalVolume* GetVolume(G4int depth=0) const;
  inline G4VSolid* GetSolid(G4int depth=0) const;
  const G4ThreeVector& GetTranslation(G4int depth=0) const;
  const G4RotationMatrix* GetRotation(G4int depth=0) const;

  inline G4int GetReplicaNumber(G4int depth=0) const;
  inline G4int GetHistoryDepth()  const;
  G4int MoveUpHistory( G4int num_levels = 1 );
    // Access methods for touchables with history

  virtual void  UpdateYourself(G4VPhysicalVolume*   pPhysVol,
			       const G4NavigationHistory* history=0); 
    // Update methods for touchables with history

 public:  // without description

  inline const G4NavigationHistory* GetHistory() const;
    // Should this method be "deprecated" ?
    // it is used now in G4Navigator::LocateGlobalPointAndSetup

  // void  ResetLevel();
       // Set the level to the top level.

 private:

  inline G4int CalculateHistoryIndex(G4int stackDepth) const;

  G4RotationMatrix frot;
  G4ThreeVector ftlate;
  G4NavigationHistory fhistory;
};

#include "G4TouchableHistory.icc"
#endif
