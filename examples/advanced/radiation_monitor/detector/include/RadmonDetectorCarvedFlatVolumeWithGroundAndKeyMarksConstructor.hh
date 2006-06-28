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
// File name:     RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.hh,v 1.2 2006-06-28 13:46:10 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a carved box with a ground pad and marks
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonTDetectorVolumesWithKeyMarksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithKeyMarksDecorator<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> > > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor(const RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH */
