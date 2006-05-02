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
// $Id: G4AssemblyTriplet.hh,v 1.1 2006-05-02 11:35:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

    G4LogicalVolume* GetVolume() const;
      // Retrieve the logical volume reference

    void SetVolume( G4LogicalVolume* pVolume );
      // Update the logical volume reference

    G4AssemblyVolume* GetAssembly() const;
      // Retrieve the assembly volume reference

    void SetAssembly( G4AssemblyVolume* pAssembly );
      // Update the assembly volume reference

    G4ThreeVector GetTranslation() const;
      // Retrieve the logical volume translation

    void SetTranslation( G4ThreeVector& pVolume );
      // Update the logical volume translation

    G4RotationMatrix* GetRotation() const;
      // Retrieve the logical volume rotation

    void SetRotation( G4RotationMatrix* pVolume );
      // Update the logical volume rotation

    G4bool IsReflection() const;
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
