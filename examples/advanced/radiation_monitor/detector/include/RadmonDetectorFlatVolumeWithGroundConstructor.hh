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
// File name:     RadmonDetectorFlatVolumeWithGroundConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithGroundConstructor.hh,v 1.2 2006-06-28 13:47:02 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithGroundConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithGroundConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithGroundConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithGroundConstructor(const RadmonDetectorFlatVolumeWithGroundConstructor & copy);
   RadmonDetectorFlatVolumeWithGroundConstructor & operator=(const RadmonDetectorFlatVolumeWithGroundConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithGroundConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH */
