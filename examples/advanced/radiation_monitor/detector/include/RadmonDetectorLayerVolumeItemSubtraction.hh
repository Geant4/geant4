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
// File name:     RadmonDetectorLayerVolumeItemSubtraction.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItemSubtraction.hh,v 1.3 2006/06/29 16:10:25 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
