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
// File name:     RadmonDetectorCarvedFlatVolumeWithPinsConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithPinsConstructor.hh,v 1.2 2006-06-28 13:46:28 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders and pins
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHPINSCONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHPINSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithPinsDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithPinsConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithPinsDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithPinsConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithPinsConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithPinsConstructor(const RadmonDetectorCarvedFlatVolumeWithPinsConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithPinsConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithPinsConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithPinsConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHPINSCONSTRUCTOR_HH */
