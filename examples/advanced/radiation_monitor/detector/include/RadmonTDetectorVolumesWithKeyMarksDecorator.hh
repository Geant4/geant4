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
// File name:     RadmonTDetectorVolumesWithKeyMarksDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithKeyMarksDecorator.hh,v 1.3 2006/06/29 16:12:36 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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
