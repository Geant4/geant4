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
/// @file VoxelParam.hh
/// @brief Define voxel parameterization

#ifndef VOXEL_PARAM_H
#define VOXEL_PARAM_H

#include "G4VPVParameterisation.hh"

class VoxelParam : public G4VPVParameterisation {
public:
  VoxelParam();
  ~VoxelParam();

  virtual void ComputeTransformation(const G4int id,
                                     G4VPhysicalVolume* vol) const;

  virtual void ComputeDimensions(G4Box& box, const G4int id,
                                 const G4VPhysicalVolume* vol) const;

  virtual void ComputeDimensions
  (G4Trd&,const G4int,const G4VPhysicalVolume*) const  {}

  virtual void ComputeDimensions
  (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Ellipsoid&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Para&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}

  virtual void ComputeDimensions
  (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}

};

#endif

