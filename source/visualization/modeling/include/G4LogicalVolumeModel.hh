// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeModel.hh,v 1.1 1999-10-04 15:37:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th July 1999.
// Model for logical volumes.

// It describes a logical volume and its daughters to any depth -
// usually only the first by default.

#ifndef G4LOGICALVOLUMEMODEL_HH
#define G4LOGICALVOLUMEMODEL_HH

#include "G4PhysicalVolumeModel.hh"

#include "globals.hh"
#include "G4Transform3D.hh"

class G4LogicalVolume;
class G4ModelingParameters;

class G4LogicalVolumeModel: public G4PhysicalVolumeModel {

public:

  G4LogicalVolumeModel
  (G4LogicalVolume*,
   G4int soughtDepth = 1,
   const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* = 0);

  virtual ~G4LogicalVolumeModel ();

  G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  void DescribeYourselfTo (G4VGraphicsScene&);

protected:

  /////////////////////////////////////////////////////////
  // Data members...

  G4LogicalVolume* fpLV;

};

#endif
