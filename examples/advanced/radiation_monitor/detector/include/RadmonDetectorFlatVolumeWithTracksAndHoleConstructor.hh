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
// File name:     RadmonDetectorFlatVolumeWithTracksAndHoleConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithTracksAndHoleConstructor.hh,v 1.3 2006/06/29 16:09:59 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Builds a box with tracks and a hole
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHTRACKSANDHOLECONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHTRACKSANDHOLECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithTracksDecorator.hh"
 #include "RadmonTDetectorVolumesWithHoleDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithTracksAndHoleConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithTracksDecorator<RadmonTDetectorVolumesWithHoleDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithTracksAndHoleConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithTracksAndHoleConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithTracksAndHoleConstructor(const RadmonDetectorFlatVolumeWithTracksAndHoleConstructor & copy);
   RadmonDetectorFlatVolumeWithTracksAndHoleConstructor & operator=(const RadmonDetectorFlatVolumeWithTracksAndHoleConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithTracksAndHoleConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHTRACKSANDHOLECONSTRUCTOR_HH */
