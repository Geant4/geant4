// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITreeSceneHandler.hh,v 1.4 2001-06-04 15:02:11 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy to standard output as
//   ASCII stream.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#ifndef G4ASCIITREESCENEHANDLER_HH
#define G4ASCIITREESCENEHANDLER_HH

#include "G4VTreeSceneHandler.hh"

#include "g4std/set"

class G4VPhysicalVolume;
class G4ModelingParameters;

class G4ASCIITreeSceneHandler: public G4VTreeSceneHandler {
public:
  G4ASCIITreeSceneHandler(G4VGraphicsSystem& system,
			  const G4String& name);
  virtual ~G4ASCIITreeSceneHandler();
  void BeginModeling();
  void EndModeling();
protected:
  void Dump(const G4VSolid&);
  const G4ModelingParameters* fpOriginalMP;  // Keeps pointer to original.
  G4ModelingParameters* fpNonCullingMP;      // For temporary non-culling.
  G4std::set<G4LogicalVolume*> fLVSet;
  typedef G4std::set<G4LogicalVolume*>::iterator LVSetIterator;
  G4std::set<G4VPhysicalVolume*> fReplicaSet;
  typedef G4std::set<G4VPhysicalVolume*>::iterator ReplicaSetIterator;
};

#endif
