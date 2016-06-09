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
// File name:     RadmonTDetectorCarvedVolumesDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorCarvedVolumesDecorator.hh,v 1.5 2006/06/29 16:12:13 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
