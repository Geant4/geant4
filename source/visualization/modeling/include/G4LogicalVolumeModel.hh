// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeModel.hh,v 1.4 2000-04-12 13:02:36 johna Exp $
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

  // This called from G4PhysicalVolumeModel::DescribeAndDescend by the
  // virtual function mechanism.
  void DescribeSolid
  (const G4Transform3D& theAT,
   G4VSolid* pSol,
   const G4VisAttributes* pVisAttribs,
   G4VGraphicsScene& sceneHandler);

  /////////////////////////////////////////////////////////
  // Data members...

  G4LogicalVolume* fpLV;

};

#endif
