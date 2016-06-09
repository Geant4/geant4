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
// File name:     RadmonTDetectorVolumesWithPadsDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithPadsDecorator.hh,v 1.3 2006/06/29 16:12:40 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Decorates with a ground plane
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHPADSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHPADSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"
 #include "RadmonDetectorPadsData.hh"
 #include <stack>

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithPadsDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithPadsDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithPadsDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithPadsDecorator();
                                                RadmonTDetectorVolumesWithPadsDecorator(const  RadmonTDetectorVolumesWithPadsDecorator & copy);
   RadmonTDetectorVolumesWithPadsDecorator &    operator=(const  RadmonTDetectorVolumesWithPadsDecorator & copy);

  // Private datatypes
   typedef std::stack<G4VisAttributes *>        AttributesStack;
   typedef std::list<RadmonDetectorPadsData>    PadsDataList;

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
   AttributesStack                              attributesStack;
   PadsDataList                                 padsDataList;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithPadsDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHPADSDECORATOR_HH */
