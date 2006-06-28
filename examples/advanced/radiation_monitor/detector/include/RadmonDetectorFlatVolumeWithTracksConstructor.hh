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
// File name:     RadmonDetectorFlatVolumeWithTracksConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithTracksConstructor.hh,v 1.2 2006-06-28 13:47:22 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with tracks
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithTracksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithTracksConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithTracksDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithTracksConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithTracksConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithTracksConstructor(const RadmonDetectorFlatVolumeWithTracksConstructor & copy);
   RadmonDetectorFlatVolumeWithTracksConstructor & operator=(const RadmonDetectorFlatVolumeWithTracksConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithTracksConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH */
