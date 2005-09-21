//
// File name:     RadmonDetectorLayerVolumesList.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.hh,v 1.2 2005-09-21 14:52:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Collection of solids, physical and logical volumes
//

#ifndef   RADMONDETECTORLAYERVOLUMESLIST_HH
 #define  RADMONDETECTORLAYERVOLUMESLIST_HH

 // Include files
 #include "globals.hh"
 #include <vector>
 
 // Forward declaration
 class RadmonDetectorLayerVolumeItem;

 class RadmonDetectorLayerVolumesList
 {
  public:
   inline                                       RadmonDetectorLayerVolumesList();
   inline                                      ~RadmonDetectorLayerVolumesList();

   inline G4int                                 GetNItems(void) const;
   inline G4bool                                Empty(void) const;

   inline const RadmonDetectorLayerVolumeItem * GetItem(G4int index) const;
   inline RadmonDetectorLayerVolumeItem *       GetItem(G4int index);

   RadmonDetectorLayerVolumeItem *              AppendItem(void);

   void                                         RemoveItem(G4int index);
   void                                         RemoveItemsByRange(G4int first, G4int last);
   void                                         RemoveAllItems(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorLayerVolumesList(const RadmonDetectorLayerVolumesList & copy);
   RadmonDetectorLayerVolumesList &             operator=(const RadmonDetectorLayerVolumesList & copy);

  // Private data types
   typedef std::vector<RadmonDetectorLayerVolumeItem *> VolumesVector;

  // Private methods
   void                                         RemoveRange(VolumesVector::iterator begin, VolumesVector::iterator end);
   
  // Private attributes
   VolumesVector                                volumesVector;
 };
 
 // Inline implementations
 #include "RadmonDetectorLayerVolumesList.icc"
#endif /* RADMONDETECTORLAYERVOLUMESLIST_HH */
