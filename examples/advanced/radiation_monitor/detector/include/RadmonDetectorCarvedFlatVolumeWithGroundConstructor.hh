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
// File name:     RadmonDetectorCarvedFlatVolumeWithGroundConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithGroundConstructor.hh,v 1.2 2006-06-28 13:46:19 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithGroundConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithGroundConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithGroundConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithGroundConstructor(const RadmonDetectorCarvedFlatVolumeWithGroundConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithGroundConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithGroundConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithGroundConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH */
