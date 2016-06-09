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
// File name:     RadmonVDetectorLayerVolumeItemOperation.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayerVolumeItemOperation.hh,v 1.3 2006/06/29 16:13:20 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Operates on solids
//

#ifndef   RADMONVDETECTORLAYERVOLUMEITEMOPERATION_HH
 #define  RADMONVDETECTORLAYERVOLUMEITEMOPERATION_HH
 
 // Include files
 #include "G4RotationMatrix.hh"
 #include "G4ThreeVector.hh"
 #include <stack>

 // Forward declarations
 class RadmonDetectorLayerVolumeItem;
 class G4VSolid;

 class RadmonVDetectorLayerVolumeItemOperation
 {
  public:
   void                                         Initialize(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true);
   void                                         Initialize(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true);
   void                                         Initialize(G4VSolid * solid, G4bool ownSolid=true);
   void                                         Initialize(const RadmonDetectorLayerVolumeItem * item);

   G4VSolid *                                   ApplyTo(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true, G4bool ownResult=true);
   inline G4VSolid *                            ApplyTo(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true, G4bool ownResult=true);
   inline G4VSolid *                            ApplyTo(G4VSolid * solid, G4bool ownSolid=true, G4bool ownResult=true);
   void                                         ApplyTo(RadmonDetectorLayerVolumeItem * item, G4bool ownResult=true);

  protected:
   inline                                       RadmonVDetectorLayerVolumeItemOperation();
                                                RadmonVDetectorLayerVolumeItemOperation(G4VSolid * solid, const G4RotationMatrix & rotation, const G4ThreeVector & position, G4bool ownSolid=true);
                                                RadmonVDetectorLayerVolumeItemOperation(G4VSolid * solid, const G4ThreeVector & position, G4bool ownSolid=true);
                                                RadmonVDetectorLayerVolumeItemOperation(G4VSolid * solid, G4bool ownSolid=true);
                                                RadmonVDetectorLayerVolumeItemOperation(const RadmonDetectorLayerVolumeItem * item);

                                               ~RadmonVDetectorLayerVolumeItemOperation();
   virtual G4VSolid *                           Operate(G4VSolid * left, G4VSolid * right, G4RotationMatrix * relativeRotation, const G4ThreeVector & relativePosition) = 0;
   
   inline G4RotationMatrix *                    AllocateMatrix();
   inline void                                  OwnSolid(G4VSolid * solid);

  private:
   inline void                                  Validate();

   void                                         Merge(const G4RotationMatrix & rightRotation, const G4ThreeVector & rightPosition, G4RotationMatrix & relativeRotation, G4ThreeVector & relativePosition) const;
   void                                         Absolute(const RadmonDetectorLayerVolumeItem * item, G4RotationMatrix & rotation, G4ThreeVector & position) const;
   
  // Hidden constructors and operators
                                                RadmonVDetectorLayerVolumeItemOperation(const RadmonVDetectorLayerVolumeItemOperation & copy);
   RadmonVDetectorLayerVolumeItemOperation &    operator=(const RadmonVDetectorLayerVolumeItemOperation & copy);

  // Private data types
   typedef std::stack<G4VSolid *>               OwnedSolids;
   typedef std::stack<G4RotationMatrix>         OwnedMatrices;
  
  // Private attributes
   OwnedSolids                                  ownedSolids;
   OwnedMatrices                                ownedMatrices;
   G4RotationMatrix                             leftRotation;
   G4RotationMatrix                             invLeftRotation;
   G4ThreeVector                                leftPosition;
   G4VSolid *                                   leftSolid;
 };
 
 // Inline implementations
 #include "RadmonVDetectorLayerVolumeItemOperation.icc"
#endif /* RADMONVDETECTORLAYERVOLUMEITEMOPERATION_HH */
