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
// File name:     RadmonTDetectorLayerConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorLayerConstructor.hh,v 1.3 2006-06-28 13:49:05 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
