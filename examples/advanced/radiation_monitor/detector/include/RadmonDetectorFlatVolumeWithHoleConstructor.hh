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
// File name:     RadmonDetectorFlatVolumeWithHoleConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithHoleConstructor.hh,v 1.2 2006-06-28 13:47:06 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with a hole
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHHOLECONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHHOLECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithHoleDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithHoleConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithHoleDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithHoleConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithHoleConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithHoleConstructor(const RadmonDetectorFlatVolumeWithHoleConstructor & copy);
   RadmonDetectorFlatVolumeWithHoleConstructor & operator=(const RadmonDetectorFlatVolumeWithHoleConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithHoleConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHHOLECONSTRUCTOR_HH */
