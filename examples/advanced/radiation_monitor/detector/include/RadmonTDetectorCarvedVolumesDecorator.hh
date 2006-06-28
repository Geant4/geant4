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
// File name:     RadmonTDetectorCarvedVolumesDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorCarvedVolumesDecorator.hh,v 1.4 2006-06-28 13:49:01 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Carves the borders of a component with the following method
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH
 #define  RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH

 // Include files
 #include <stack>
 #include "RadmonDetectorLayerVolumeItemSubtraction.hh"
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"
 #include "G4RotationMatrix.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;
 class G4VAttributes;

 template <class LayerVolumesComponent>
 class RadmonTDetectorCarvedVolumesDecorator
 {
  public:
                                                RadmonTDetectorCarvedVolumesDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorCarvedVolumesDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorCarvedVolumesDecorator();
                                                RadmonTDetectorCarvedVolumesDecorator(const RadmonTDetectorCarvedVolumesDecorator & copy);
   RadmonTDetectorCarvedVolumesDecorator &      operator=(const RadmonTDetectorCarvedVolumesDecorator & copy);

  // Private data types
   typedef std::stack<G4VSolid *>               OwnedSolids;
  
  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   OwnedSolids                                  ownedSolids;
   G4VisAttributes *                            visAttributes;
   RadmonDetectorLayerVolumeItemSubtraction     subtraction;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
   G4RotationMatrix                             identity;
  };
 
 // Inline implementations
 #include "RadmonTDetectorCarvedVolumesDecorator.icc"
#endif /* RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH */
