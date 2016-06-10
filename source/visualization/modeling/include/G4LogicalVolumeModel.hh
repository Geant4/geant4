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
// $Id: G4LogicalVolumeModel.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
   G4bool booleans = true,
   G4bool voxels = true,
   G4bool readout = true,
   const G4Transform3D& modelTransformation = G4Transform3D(),
   const G4ModelingParameters* = 0);

  virtual ~G4LogicalVolumeModel ();

  void DescribeYourselfTo (G4VGraphicsScene&);

  G4bool Validate (G4bool) {return true;}

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
  G4bool fBooleans;  // Flag for drawing boolean components.
  G4bool fVoxels;    // Flag for drawing voxels.
  G4bool fReadout;   // Flag for drawing readout geometry.

};

#endif
