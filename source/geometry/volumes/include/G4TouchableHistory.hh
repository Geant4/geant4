// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TouchableHistory.hh,v 1.1 1999-01-07 16:08:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4TouchableHistory     Paul Kent August 1996
//
// Object representing a touchable detector element, and
// its history in the geomtrical hierarchy, including
// its net resultant local->global transform

#ifndef G4TOUCHABLEHISTORY_HH
#define G4TOUCHABLEHISTORY_HH

#include "G4VTouchable.hh"

#include "G4NavigationHistory.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4TouchableHistory : public G4VTouchable
{
public:
  G4TouchableHistory(const G4NavigationHistory &history);
  // The default constructor produces a touchable-history of 
  //  'zero-depth', ie an "unphysical" and not very unusable one.
  //  It is for initialisation only .  
  G4TouchableHistory(); 
  ~G4TouchableHistory();
  G4VPhysicalVolume* GetVolume(G4int depth=0) const;
  G4VSolid* GetSolid(G4int depth=0) const;
  const G4ThreeVector& GetTranslation(G4int depth=0) const;
  const G4RotationMatrix*  GetRotation(G4int depth=0) const;

  // New access methods for touchables with history
  G4int GetReplicaNumber(G4int depth=0) const;
  G4int GetHistoryDepth()  const;
  G4int MoveUpHistory( G4int num_levels = 1 );
  // void  ResetLevel();  // Set the level to the top level.

  // Update methods for touchables with history
  virtual void  UpdateYourself(     G4VPhysicalVolume*   pPhysVol,
			      const G4NavigationHistory* history=NULL); 

  // Should this method be "depricated" ?
  //    it is used now in G4Navigator::LocateGlobalPointAndSetup
  //
  const G4NavigationHistory* GetHistory() const;

private:
  G4int CalculateHistoryIndex(G4int stackDepth) const;

private:
  G4RotationMatrix frot;
  G4ThreeVector ftlate;
  G4NavigationHistory fhistory;
};

#include "G4TouchableHistory.icc"
#endif
