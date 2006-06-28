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
// File name:     RadmonDetectorFlatVolumeWithPinsConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithPinsConstructor.hh,v 1.2 2006-06-28 13:47:14 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with pins
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHPINSCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHPINSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithPinsDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithPinsConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithPinsDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithPinsConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithPinsConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithPinsConstructor(const RadmonDetectorFlatVolumeWithPinsConstructor & copy);
   RadmonDetectorFlatVolumeWithPinsConstructor & operator=(const RadmonDetectorFlatVolumeWithPinsConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithPinsConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHPINSCONSTRUCTOR_HH */
