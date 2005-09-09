//
// File name:     RadmonDetectorLayerVolumesList.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Collection of solids, physical and logical volumes
//

#ifndef   RADMONDETECTORLAYERVOLUMESLIST_HH
 #define  RADMONDETECTORLAYERVOLUMESLIST_HH

 // Include files
 #include "globals.hh"
 
 // Forward declaration
 class RadmonDetectorLayerVolumeItem;

 class RadmonDetectorLayerVolumesList
 {
  public:
                                                RadmonDetectorLayerVolumesList();
                                               ~RadmonDetectorLayerVolumesList();

   G4Int                                        GetNItems(void) const;
   G4bool                                       Empty(void) const;

   const RadmonDetectorLayerVolumeItem *        GetItem(G4int index) const;
   RadmonDetectorLayerVolumeItem *              GetItem(G4int index);

   RadmonDetectorLayerVolumeItem *              CreateItem(void);

   void                                         RemoveItem(G4int index);
   void                                         RemoveItemsByRange(G4int first, G4int last);
   void                                         RemoveAllItems(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorLayerVolumesList(const RadmonDetectorLayerVolumesList & copy);
   RadmonDetectorLayerVolumesList &             operator=(const RadmonDetectorLayerVolumesList & copy);

  // Private attributes
   std::list<RadmonDetectorLayerVolumeItem *>   volumesList;
 };
#endif /* RADMONDETECTORLAYERVOLUMESLIST_HH */
