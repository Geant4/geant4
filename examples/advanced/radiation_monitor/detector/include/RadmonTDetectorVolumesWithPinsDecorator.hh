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
// File name:     RadmonTDetectorVolumesWithPinsDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithPinsDecorator.hh,v 1.3 2006/06/29 16:12:44 gunter Exp $
// Tag:           $Name: geant4-09-02 $
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
