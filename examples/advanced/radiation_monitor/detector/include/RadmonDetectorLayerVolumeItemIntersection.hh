//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonDetectorLayerVolumeItemIntersection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemIntersection.hh,v 1.2 2006-06-28 13:47:41 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Intersection of two solids
//

#ifndef   RADMONDETECTORLAYERVOLUMEITEMINTERSECTION_HH
 #define  RADMONDETECTORLAYERVOLUMEITEMINTERSECTION_HH
 
 // Include files
 #include "RadmonVDetectorLayerVolumeItemOperation.hh"

 class RadmonDetectorLayerVolumeItemIntersection : public RadmonVDetectorLayerVolumeItemOperation
 {
  public:
   inline                                       RadmonDetectorLayerVolumeItemIntersection();
   inline                                       RadmonDetectorLayerVolumeItemIntersection(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemIntersection(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemIntersection(G4VSolid * solid, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemIntersection(const RadmonDetectorLayerVolumeItem * item);
   inline                                      ~RadmonDetectorLayerVolumeItemIntersection();

  protected:
   virtual G4VSolid *                           Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition);
 };

 // Inline implementations
 #include "RadmonDetectorLayerVolumeItemIntersection.icc"
#endif /* RADMONDETECTORLAYERVOLUMEITEMINTERSECTION_HH */
