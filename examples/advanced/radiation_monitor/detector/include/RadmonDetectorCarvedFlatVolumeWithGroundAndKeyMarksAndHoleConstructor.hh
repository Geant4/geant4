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
// File name:     RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.hh,v 1.2 2006-06-28 13:46:06 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a carved box with a ground pad and marks and hole
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSANDHOLECONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSANDHOLECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonTDetectorVolumesWithHoleDecorator.hh"
 #include "RadmonTDetectorVolumesWithKeyMarksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithHoleDecorator<RadmonTDetectorVolumesWithKeyMarksDecorator<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> > > > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor(const RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSANDHOLECONSTRUCTOR_HH */
