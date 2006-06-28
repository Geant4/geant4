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
// File name:     RadmonTDetectorVolumesWithPadsDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithPadsDecorator.hh,v 1.2 2006-06-28 13:49:39 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
