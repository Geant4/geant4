#ifndef ASSEMBLYTRIPLET_H
#define ASSEMBLYTRIPLET_H

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

/**
 * Class:       G4AssemblyTriplet
 * Description: A class to help place logical volumes inside a generic
 *              containers (like STL vector ) together with information
 *              about its rotation and placement. How to interpret the
 *              rotation and translation depends on the class which uses
 *              a container of these triplets. The first class using
 *              G4AssemblyTriplet is G4AssemblyVolume class.
 *              The pointer to the logical volume is copied so this class
 *              does not take its ownership and does not delete the object
 *              behind.
 *
 * Author:      Radovan Chytracek
 * Version:     1.0
 * Date:        November 2000
 */
class G4AssemblyTriplet {

public:

    // Default constructor
    G4AssemblyTriplet();    

    // An explicit constructor
    G4AssemblyTriplet( G4LogicalVolume* pVolume
                      ,G4ThreeVector& translation
                      ,G4RotationMatrix* pRotation
                     );
    
    // Copy constructor
    G4AssemblyTriplet( const G4AssemblyTriplet& second );

    // Destructor
    ~G4AssemblyTriplet();    

    // Retrieve the logical volume reference
    G4LogicalVolume* GetVolume() const;

    // Update the logical volume reference
    void SetVolume(G4LogicalVolume* pVolume);

    // Retrieve the logical volume translation
    G4ThreeVector GetTranslation() const;

    // Update the logical volume translation
    void SetTranslation(G4ThreeVector& pVolume);

    // Retrieve the logical volume rotation
    G4RotationMatrix* GetRotation() const;

    // Update the logical volume rotation
    void SetRotation(G4RotationMatrix* pVolume);
    
    // Assignment
    G4AssemblyTriplet& operator=( const G4AssemblyTriplet& second );
    

private:

    // A logical volume
    G4LogicalVolume*  fVolume;

    // A logical volume translation
    G4ThreeVector     fTranslation;

    // A logical volume rotation
    G4RotationMatrix* fRotation;
};

#include "G4AssemblyTriplet.icc"

#endif //ASSEMBLYTRIPLET_H

