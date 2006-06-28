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
// File name:     RadmonTDetectorVolumesWithPinsDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithPinsDecorator.hh,v 1.2 2006-06-28 13:49:43 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with the pin out
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHPINSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHPINSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithPinsDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithPinsDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithPinsDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithPinsDecorator();
                                                RadmonTDetectorVolumesWithPinsDecorator(const  RadmonTDetectorVolumesWithPinsDecorator & copy);
   RadmonTDetectorVolumesWithPinsDecorator &    operator=(const  RadmonTDetectorVolumesWithPinsDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   G4VisAttributes *                            visAttributes;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithPinsDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHPINSDECORATOR_HH */
