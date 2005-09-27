//
// File name:     RadmonDetectorLayerVolumeItemUnion.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemUnion.hh,v 1.1 2005-09-27 13:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Union of two solids
//

#ifndef   RADMONDETECTORLAYERVOLUMEITEMUNION_HH
 #define  RADMONDETECTORLAYERVOLUMEITEMUNION_HH
 
 // Include files
 #include "RadmonVDetectorLayerVolumeItemOperation.hh"

 class RadmonDetectorLayerVolumeItemUnion : public RadmonVDetectorLayerVolumeItemOperation
 {
  public:
   inline                                       RadmonDetectorLayerVolumeItemUnion();
   inline                                       RadmonDetectorLayerVolumeItemUnion(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemUnion(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemUnion(G4VSolid * solid, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemUnion(const RadmonDetectorLayerVolumeItem * item);
   inline                                      ~RadmonDetectorLayerVolumeItemUnion();

  protected:
   virtual G4VSolid *                           Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition);
 };

 // Inline implementations
 #include "RadmonDetectorLayerVolumeItemUnion.icc"
#endif /* RADMONDETECTORLAYERVOLUMEITEMUNION_HH */
