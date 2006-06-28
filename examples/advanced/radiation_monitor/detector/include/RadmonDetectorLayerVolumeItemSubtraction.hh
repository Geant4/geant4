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
// File name:     RadmonDetectorLayerVolumeItemSubtraction.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemSubtraction.hh,v 1.2 2006-06-28 13:47:45 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Subtraction of two solids
//

#ifndef   RADMONDETECTORLAYERVOLUMEITEMSUBTRACTION_HH
 #define  RADMONDETECTORLAYERVOLUMEITEMSUBTRACTION_HH
 
 // Include files
 #include "RadmonVDetectorLayerVolumeItemOperation.hh"

 class RadmonDetectorLayerVolumeItemSubtraction : public RadmonVDetectorLayerVolumeItemOperation
 {
  public:
   enum Mode
   {
    leftMinusRight,
    rightMinusLeft
   };
  
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, G4VSolid * solid, G4bool ownSolid=true);
   inline                                       RadmonDetectorLayerVolumeItemSubtraction(Mode mode, const RadmonDetectorLayerVolumeItem * item);
   inline                                      ~RadmonDetectorLayerVolumeItemSubtraction();

  protected:
   virtual G4VSolid *                           Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition);
   
  private:
   Mode                                         opMode;
 };

 // Inline implementations
 #include "RadmonDetectorLayerVolumeItemSubtraction.icc"
#endif /* RADMONDETECTORLAYERVOLUMEITEMSUBTRACTION_HH */
