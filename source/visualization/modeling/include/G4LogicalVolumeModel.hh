// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeModel.hh,v 1.3 1999-12-15 14:54:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th July 1999.
//
// Class Description:
//
// Model for logical volumes.  It describes a logical volume and its
// daughters to any depth - usually only the first by default.
//
// Inherits from G4PhysicalVolumeModel; for more information see that
// class description.

#ifndef G4LOGICALVOLUMEMODEL_HH
#define G4LOGICALVOLUMEMODEL_HH

#include "G4PhysicalVolumeModel.hh"

#include "globals.hh"
#include "G4Transform3D.hh"

class G4LogicalVolume;
class G4ModelingParameters;

class G4LogicalVolumeModel: public G4PhysicalVolumeModel {

public: // With description

  G4LogicalVolumeModel
  (G4LogicalVolume*,
   G4int soughtDepth = 1,
   const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* = 0);

  virtual ~G4LogicalVolumeModel ();

  void DescribeYourselfTo (G4VGraphicsScene&);

  G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

protected:

  /////////////////////////////////////////////////////////
  // Data members...

  G4LogicalVolume* fpLV;

};

#endif
