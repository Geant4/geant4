//
// File name:     RadmonTDetectorLayerConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorLayerConstructor.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Generates a layer using a component with the following methods
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORLAYERCONSTRUCTOR_HH
 #define  RADMONTDETECTORLAYERCONSTRUCTOR_HH

 // Include files
 #include "RadmonVDetectorLabelledEntityConstructor.hh"

 #include "globals.hh" 

 // Forward declaration
 class RadmonDetectorLayerVolumesList;

 template <class LayerVolumesComponent> class RadmonTDetectorLayerConstructor : public RadmonVDetectorLabelledEntityConstructor
 {
  public:
                                                RadmonTDetectorLayerConstructor(const G4String & label);
                                               ~RadmonTDetectorLayerConstructor();

   virtual G4LogicalVolume *                    ConstructLogicalVolume(void);
   virtual void                                 DestructLogicalVolume(void);

   G4double                                     GetWidth(void) const;
   G4double                                     GetHeight(void) const;
   G4double                                     GetThickness(void) const;

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorLayerConstructor();
                                                RadmonTDetectorLayerConstructor(const RadmonTDetectorLayerConstructor & copy);
   RadmonTDetectorLayerConstructor &            operator=(const RadmonTDetectorLayerConstructor & copy);

  // Private attributes
   RadmonDetectorLayerVolumesList *             volumesList;
 };
#endif /* RADMONTDETECTORLAYERCONSTRUCTOR_HH */
