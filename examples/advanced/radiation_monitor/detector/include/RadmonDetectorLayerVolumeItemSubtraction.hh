//
// File name:     RadmonDetectorLayerVolumeItemSubtraction.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemSubtraction.hh,v 1.1 2005-09-27 13:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Subtraction of two solids
//

#ifndef   RADMONDETECTORLAYERVOLUMEITEMSUBTRACTION_HH
 #define  RADMONDETECTORLAYERVOLUMEITEMSUBTRACTION_HH
 
 // Include files
 #include "RadmonVDetectorLayerVolumeItemOperation.hh"

 class RadmonDetectorLayerVolumeItemSubtraction : public RadmonVDetectorLayerVolumeItemOperation
 {
  public:
   enum Mode
   {
    leftMinusRight,
    rightMinusLeft
   };
  
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, G4VSolid * solid, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, const RadmonDetectorLayerVolumeItem * item);
   inline                                      ~RadmonDetectorLayerVolumeItemSubtraction();

  protected:
   virtual G4VSolid *                           Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition);
   
  private:
   Mode                                         opMode;
 };

 // Inline implementations
 #include "RadmonDetectorLayerVolumeItemSubtraction.icc"
#endif /* RADMONDETECTORLAYERVOLUMEITEMSUBTRACTION_HH */
