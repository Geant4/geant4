// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VFlavoredParallelWorld.hh,v 1.1 1999-04-15 15:43:11 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Abstract interface for GEANT4 Flavored Parallel World.
// P. Mora de Freitas & M. Verderi 14/April/1999.
//

#ifndef G4VFLAVOREDPARALLELWORLD_HH
#define G4VFLAVOREDPARALLELWORLD_HH

class G4VFlavoredParallelWorld {

public:

  virtual ~G4VFlavoredParallelWorld () {}

  // G4VFlavoredParallelWorld Interface for visualisation.

  virtual 
  G4VPhysicalVolume* GetThePhysicalVolumeWorld() const =0;
};
#endif
