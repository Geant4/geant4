// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FlavoredParallelWorldModel.hh,v 1.2 1999-01-11 00:48:42 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// P. Mora de Freitas et M.Verderi - 19 June 1998.
// Model for flavored parallel world volumes.

#ifndef G4FLAVOREDPARALLELWORLDMODEL_HH
#define G4FLAVOREDPARALLELWORLDMODEL_HH

#include "G4PhysicalVolumeModel.hh"

class G4FlavoredParallelWorld;

class G4FlavoredParallelWorldModel : public G4PhysicalVolumeModel {
  
public:
  
  G4FlavoredParallelWorldModel
  (G4FlavoredParallelWorld* FPW,
   G4int soughtDepth = G4PhysicalVolumeModel::UNLIMITED,
   const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* mp = 0);

  ~G4FlavoredParallelWorldModel ();
  
private:

  G4FlavoredParallelWorld* theFlavoredParallelWorld;

};

#endif
