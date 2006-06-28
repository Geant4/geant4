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
// File name:     RadmonTDetectorVolumesWithHoleDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithHoleDecorator.hh,v 1.4 2006-06-28 13:49:31 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with a hole a component with the following method
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemSubtraction.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithHoleDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithHoleDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithHoleDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithHoleDecorator();
                                                RadmonTDetectorVolumesWithHoleDecorator(const  RadmonTDetectorVolumesWithHoleDecorator & copy);
   RadmonTDetectorVolumesWithHoleDecorator &    operator=(const  RadmonTDetectorVolumesWithHoleDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   RadmonDetectorLayerVolumeItemSubtraction     operation;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithHoleDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH */
