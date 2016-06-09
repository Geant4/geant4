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
// File name:     RadmonDetectorLayerVolumesList.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.hh,v 1.4 2006/06/29 16:10:33 gunter Exp $
// Tag:           $Name: geant4-08-02 $
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
