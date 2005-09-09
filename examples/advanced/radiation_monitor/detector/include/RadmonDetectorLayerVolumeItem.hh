//
// File name:     RadmonDetectorLayerVolumesItem.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItem.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   solids, physical and logical volumes
//

#ifndef   RADMONDETECTORLAYERVOLUMEITEM_HH
 #define  RADMONDETECTORLAYERVOLUMEITEM_HH
 
 // Include files
 #include "G4RotationMatrix.hh"
 #include "G4ThreeVector.hh"
 #include "G4String.hh"
 
 // Forward declarations
 class G4VSolid;
 class G4Material;
 class G4LogicalVolume;
 class G4VPhysicalVolume;

 class RadmonDetectorLayerVolumeItem
 {
  public:
                                                RadmonDetectorLayerVolumeItem();
                                               ~RadmonDetectorLayerVolumeItem();

   void                                         SetSolid(G4VSolid * solid);
   void                                         SetMaterial(G4Material * material);
   void                                         SetName(const G4String & name);
   void                                         SetPosition(const G4ThreeVector & position);
   void                                         SetRotation(const G4RotationMatrix & rotation);
   void                                         SetMotherVolumeItem(RadmonDetectorLayerVolumeItem * item);

   G4VSolid *                                   GetSolid(void) const;
   G4Material *                                 GetMaterial(void) const;
   const G4String &                             GetName(void) const;
   const G4ThreeVector &                        GetPosition(void) const;
   const G4RotationMatrix &                     GetRotation(void) const;
   RadmonDetectorLayerVolumeItem *              GetMotherVolumeItem(void) const;

   G4LogicalVolume *                            GetLogicalVolume(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorLayerVolumeItem(const RadmonDetectorLayerVolumeItem & copy);
   RadmonDetectorLayerVolumeItem &              operator=(const RadmonDetectorLayerVolumeItem & copy);

  // Private attributes
   G4VSolid *                                   volumeSolid;
   G4Material *                                 volumeMaterial;
   G4String                                     volumeName;
   G4RotationMatrix                             volumeRotation;
   G4ThreeVector                                volumePosition;
   G4LogicalVolume *                            volumeLogical;
   G4VPhysicalVolume *                          volumePhysical;
   RadmonDetectorLayerVolumeItem *              volumeMother;
 };
#endif /* RADMONDETECTORLAYERVOLUMEITEM_HH */
