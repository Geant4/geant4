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
// File name:     RadmonDetectorLayerVolumesItem.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItem.hh,v 1.5 2006/06/29 16:10:17 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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
 class G4VisAttributes;
 class G4VSolid;
 class G4Material;
 class G4LogicalVolume;
 class G4VPhysicalVolume;
 class G4VSensitiveDetector;

 class RadmonDetectorLayerVolumeItem
 {
  public:
   inline                                       RadmonDetectorLayerVolumeItem();
                                               ~RadmonDetectorLayerVolumeItem();

   inline void                                  SetSolid(G4VSolid * solid);
   inline void                                  SetAttributes(G4VisAttributes * attrs);
   inline void                                  SetMaterial(G4Material * material);
   inline void                                  SetName(const G4String & name);
   inline void                                  SetPosition(const G4ThreeVector & position);
   inline void                                  SetRotation(const G4RotationMatrix & rotation);
   inline void                                  SetSensitiveDetector(G4VSensitiveDetector * detector);
   inline void                                  SetMotherVolumeItem(RadmonDetectorLayerVolumeItem * item);

   inline G4VSolid *                            GetSolid(void) const;
   inline G4VisAttributes *                     GetAttributes(void) const;
   inline G4Material *                          GetMaterial(void) const;
   inline const G4String &                      GetName(void) const;
   inline const G4ThreeVector &                 GetPosition(void) const;
   inline const G4RotationMatrix &              GetRotation(void) const;
   inline G4VSensitiveDetector *                GetSensitiveDetector(void) const;
   inline RadmonDetectorLayerVolumeItem *       GetMotherVolumeItem(void) const;

   G4LogicalVolume *                            GetLogicalVolume(void);

  private:
  // Private methods
   void                                         Assertion(void);
  
  // Hidden constructors and operators
                                                RadmonDetectorLayerVolumeItem(const RadmonDetectorLayerVolumeItem & copy);
   RadmonDetectorLayerVolumeItem &              operator=(const RadmonDetectorLayerVolumeItem & copy);

  // Private attributes
   G4VSolid *                                   volumeSolid;
   G4VisAttributes *                            volumeAttributes;
   G4Material *                                 volumeMaterial;
   G4String                                     volumeName;
   G4RotationMatrix                             volumeRotation;
   G4ThreeVector                                volumePosition;
   G4LogicalVolume *                            volumeLogical;
   G4VPhysicalVolume *                          volumePhysical;
   G4VSensitiveDetector *                       volumeSensitiveDetector;
   RadmonDetectorLayerVolumeItem *              volumeMother;
 };
 
 // Inline implementations
 #include "RadmonDetectorLayerVolumeItem.icc"
#endif /* RADMONDETECTORLAYERVOLUMEITEM_HH */
