// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FlavoredParallelWorldModel.hh,v 1.4 1999-11-25 14:18:44 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// P. Mora de Freitas et M.Verderi - 19 June 1998.
//
// Class Description:
//
// Model for flavored parallel world volumes.  Inherits from
// G4PhysicalVolumeModel; for more information see that class
// description.

#ifndef G4FLAVOREDPARALLELWORLDMODEL_HH
#define G4FLAVOREDPARALLELWORLDMODEL_HH

#include "G4PhysicalVolumeModel.hh"

class G4VFlavoredParallelWorld;

class G4FlavoredParallelWorldModel : public G4PhysicalVolumeModel {
  
public: // With description
  
  G4FlavoredParallelWorldModel
  (G4VFlavoredParallelWorld* FPW,
   G4int soughtDepth = G4PhysicalVolumeModel::UNLIMITED,
   const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* mp = 0);

  ~G4FlavoredParallelWorldModel ();
  
private:

  G4VFlavoredParallelWorld* theFlavoredParallelWorld;

};

#endif
