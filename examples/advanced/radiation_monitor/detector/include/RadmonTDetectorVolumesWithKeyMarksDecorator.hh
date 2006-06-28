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
// File name:     RadmonTDetectorVolumesWithKeyMarksDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithKeyMarksDecorator.hh,v 1.2 2006-06-28 13:49:35 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with key marks
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHKEYMARKSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHKEYMARKSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithKeyMarksDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithKeyMarksDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithKeyMarksDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithKeyMarksDecorator();
                                                RadmonTDetectorVolumesWithKeyMarksDecorator(const  RadmonTDetectorVolumesWithKeyMarksDecorator & copy);
   RadmonTDetectorVolumesWithKeyMarksDecorator & operator=(const  RadmonTDetectorVolumesWithKeyMarksDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   G4VisAttributes *                            visAttributes;
   G4VSolid *                                   box;
   G4VSolid *                                   tub;
   G4RotationMatrix                             identity;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithKeyMarksDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHKEYMARKSDECORATOR_HH */
