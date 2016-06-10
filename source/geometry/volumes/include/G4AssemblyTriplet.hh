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
// $Id: G4AssemblyTriplet.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// Class G4AssemblyTriplet
//
// Class description:
//
// A class to help place logical or assembly volumes inside a generic 
// containers (like STL vector ) together with information about its rotation, 
// placement and eventually reflection.
// How to interpret the rotation and translation depends on the class which
// uses a container of these triplets. The first class using G4AssemblyTriplet
// is G4AssemblyVolume class.
// The pointer to the logical or assembly volume is copied so this class 
// does not take its ownership and does not delete the object behind.

// Author:      Radovan Chytracek
// Version:     1.0
// Date:        November 2000
//
// History:
// March 2006, I.Hrivnacova - Extended to support assembly of assemblies
//             of volumes and reflections
// ----------------------------------------------------------------------
#ifndef G4_ASSEMBLYTRIPLET_H
#define G4_ASSEMBLYTRIPLET_H 

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4LogicalVolume;
class G4AssemblyVolume;

class G4AssemblyTriplet
{
 public:  // with description

    G4AssemblyTriplet();
      // Default constructor

    G4AssemblyTriplet( G4LogicalVolume* pVolume,
                       G4ThreeVector& translation,
                       G4RotationMatrix* pRotation,
                       G4bool isReflection = false);
      // An explicit constructor for a logical volume
    
    G4AssemblyTriplet( G4AssemblyVolume* pAssembly,
                       G4ThreeVector& translation,
                       G4RotationMatrix* pRotation,
                       G4bool isReflection = false);
      // An explicit constructor for an assembly volume
    
    G4AssemblyTriplet( const G4AssemblyTriplet& second );
      // Copy constructor

    ~G4AssemblyTriplet();    
      // Destructor

    G4AssemblyTriplet& operator=( const G4AssemblyTriplet& second );
      // Assignment operator

    inline G4LogicalVolume* GetVolume() const;
      // Retrieve the logical volume reference

    inline void SetVolume( G4LogicalVolume* pVolume );
      // Update the logical volume reference

    inline G4AssemblyVolume* GetAssembly() const;
      // Retrieve the assembly volume reference

    inline void SetAssembly( G4AssemblyVolume* pAssembly );
      // Update the assembly volume reference

    inline G4ThreeVector GetTranslation() const;
      // Retrieve the logical volume translation

    inline void SetTranslation( G4ThreeVector& pVolume );
      // Update the logical volume translation

    inline G4RotationMatrix* GetRotation() const;
      // Retrieve the logical volume rotation

    inline void SetRotation( G4RotationMatrix* pVolume );
      // Update the logical volume rotation

    inline G4bool IsReflection() const;
      // Return true if the logical or assembly volume has reflection

 private:

    G4LogicalVolume*  fVolume;
      // A logical volume

    G4ThreeVector     fTranslation;
      // A logical volume translation

    G4RotationMatrix* fRotation;
      // A logical volume rotation

 private:

    // Member data for handling assemblies of assemblies and reflections

    G4AssemblyVolume* fAssembly;
      // An assembly volume

    G4bool            fIsReflection;
      // True if the logical or assembly volume has reflection  
};

#include "G4AssemblyTriplet.icc"

#endif // G4_ASSEMBLYTRIPLET_H
