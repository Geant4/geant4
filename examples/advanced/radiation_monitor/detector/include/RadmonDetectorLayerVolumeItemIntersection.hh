//
// File name:     RadmonDetectorLayerVolumeItemIntersection.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemIntersection.hh,v 1.1 2005-09-27 13:53:30 capra Exp $
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
