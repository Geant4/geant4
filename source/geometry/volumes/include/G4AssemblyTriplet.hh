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
// G4AssemblyTriplet
//
// Class description:
//
// A class to help place logical or assembly volumes inside a generic 
// container (like STL vector) together with information about its rotation, 
// placement and eventually reflection.
// How to interpret the rotation and translation depends on the class which
// uses a container of these triplets. The first class using G4AssemblyTriplet
// is G4AssemblyVolume.
// The pointer to the logical or assembly volume is copied so this class 
// does not take its ownership and does not delete the object behind.

// Author: Radovan Chytracek (CERN), November 2000
// ----------------------------------------------------------------------
#ifndef G4_ASSEMBLYTRIPLET_HH
#define G4_ASSEMBLYTRIPLET_HH

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4LogicalVolume;
class G4AssemblyVolume;

/**
 * @brief G4AssemblyTriplet is a helper class for placing logical or assembly
 * volumes inside a generic container together with information about its
 * rotation, placement and eventually reflection.
 * How to interpret the rotation and translation depends on the class which
 * uses a container of these triplets. The first class using G4AssemblyTriplet
 * is G4AssemblyVolume.
 * The pointer to the logical or assembly volume is copied so this class 
 * does not take its ownership and does not delete the object behind.
 */

class G4AssemblyTriplet
{
  public:

    /**
       Default constructor.
     */
    G4AssemblyTriplet();

    /**
     * An explicit constructor for a logical volume.
     *  @param[in] pVolume Pointer to the logical volume.
     *  @param[in] translation Translation vector.
     *  @param[in] pRotation Pointer to the rotation matrix.
     *  @param[in] isReflection Flag to specify if volume is reflected.
     */
    G4AssemblyTriplet( G4LogicalVolume* pVolume,
                       G4ThreeVector& translation,
                       G4RotationMatrix* pRotation,
                       G4bool isReflection = false );
    
    /**
     * An explicit constructor for an assembly volume.
     *  @param[in] pAssembly Pointer to the assembly volume.
     *  @param[in] translation Translation vector.
     *  @param[in] pRotation Pointer to the rotation matrix.
     *  @param[in] isReflection Flag to specify if volume is reflected.
     */
    G4AssemblyTriplet( G4AssemblyVolume* pAssembly,
                       G4ThreeVector& translation,
                       G4RotationMatrix* pRotation,
                       G4bool isReflection = false );
    
    /**
       Copy constructor.
     */
    G4AssemblyTriplet( const G4AssemblyTriplet& second );

    /**
       Destructor.
     */
    ~G4AssemblyTriplet();    

    /**
       Assignment operator.
     */
    G4AssemblyTriplet& operator=( const G4AssemblyTriplet& second );

    /**
       Retrieves the logical volume reference.
     */
    inline G4LogicalVolume* GetVolume() const;

    /**
       Updates the logical volume reference.
     */
    inline void SetVolume( G4LogicalVolume* pVolume );

    /**
       Retrieves the assembly volume reference.
     */
    inline G4AssemblyVolume* GetAssembly() const;

    /**
       Updates the assembly volume reference.
     */
    inline void SetAssembly( G4AssemblyVolume* pAssembly );

    /**
       Retrieves the logical volume translation.
     */
    inline G4ThreeVector GetTranslation() const;

    /**
       Updates the logical volume translation.
     */
    inline void SetTranslation( G4ThreeVector& pVolume );

    /**
       Retrieves the logical volume rotation.
     */
    inline G4RotationMatrix* GetRotation() const;

    /**
       Updates the logical volume rotation.
     */
    inline void SetRotation( G4RotationMatrix* pVolume );

    /**
       Returns true if the logical or assembly volume has reflection.
     */
    inline G4bool IsReflection() const;

  private:

    /** A logical volume. */
    G4LogicalVolume* fVolume = nullptr;

    /** A logical volume translation. */
    G4ThreeVector fTranslation;

    /** A logical volume rotation. */
    G4RotationMatrix* fRotation = nullptr;

    // Member data for handling assemblies of assemblies and reflections

    /** An assembly volume. */
    G4AssemblyVolume* fAssembly = nullptr;

    /** True if the logical or assembly volume has reflection. */
    G4bool fIsReflection = false;
};

#include "G4AssemblyTriplet.icc"

#endif
