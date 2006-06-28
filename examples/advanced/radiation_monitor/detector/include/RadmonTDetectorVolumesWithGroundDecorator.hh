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
// File name:     RadmonTDetectorVolumesWithGroundDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithGroundDecorator.hh,v 1.3 2006-06-28 13:49:27 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with a ground plane
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHGROUNDDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHGROUNDDECORATOR_HH

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4Box;
 class G4VisAttributes;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithGroundDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithGroundDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithGroundDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithGroundDecorator();
                                                RadmonTDetectorVolumesWithGroundDecorator(const  RadmonTDetectorVolumesWithGroundDecorator & copy);
   RadmonTDetectorVolumesWithGroundDecorator &  operator=(const  RadmonTDetectorVolumesWithGroundDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   G4Box *                                      box;
   G4VisAttributes *                            visAttributes;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithGroundDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHGROUNDDECORATOR_HH */
