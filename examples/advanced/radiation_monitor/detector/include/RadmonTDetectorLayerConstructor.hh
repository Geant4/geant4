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
// File name:     RadmonTDetectorLayerConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorLayerConstructor.hh,v 1.4 2006/06/29 16:12:17 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//
// Description:   Generates a layer using a component with the following methods
//                 - Constructor with RadmonVDetectorLabelledEntityConstructor or its base classes as argument
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

 template <class LayerVolumesComponent>
 class RadmonTDetectorLayerConstructor : public RadmonVDetectorLabelledEntityConstructor
 {
  public:
                                                RadmonTDetectorLayerConstructor(const G4String & label);
                                               ~RadmonTDetectorLayerConstructor();

   virtual G4LogicalVolume *                    ConstructLogicalVolume(void);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void) const;

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorLayerConstructor();
                                                RadmonTDetectorLayerConstructor(const RadmonTDetectorLayerConstructor & copy);
   RadmonTDetectorLayerConstructor &            operator=(const RadmonTDetectorLayerConstructor & copy);

  // Private attributes
   RadmonDetectorLayerVolumesList *             volumesList;
   LayerVolumesComponent                        component;
 };
 
 // Inline implementations
 #include "RadmonTDetectorLayerConstructor.icc"
#endif /* RADMONTDETECTORLAYERCONSTRUCTOR_HH */
