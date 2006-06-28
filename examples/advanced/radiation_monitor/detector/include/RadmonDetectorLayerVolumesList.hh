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
// File name:     RadmonDetectorLayerVolumesList.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.hh,v 1.3 2006-06-28 13:47:53 gunter Exp $
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
