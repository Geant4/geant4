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
// File name:     RadmonVDetectorLayerVolumeItemOperation.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayerVolumeItemOperation.hh,v 1.2 2006-06-28 13:50:08 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
