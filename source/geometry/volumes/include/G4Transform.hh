// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Transform.hh,v 1.3 2000-04-25 16:15:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Transform
//
// Class description:
//
// A coordinate transform, consisting of a translation and a rotation.
// When applying the forwards transformations, rotation is applied first, then
// the translation.
//
// Member data:
//
// G4RotationMatrix* frotation
//   - Pointer to the transmations rotation component
// G4ThreeVector ftranslation
//   - Translation component

// History:
// 16.07.95 P.Kent Added SetTranslation
// 13.07.95 P.Kent Initial version

#ifndef G4TRANSFORM_HH
#define G4TRANSFORM_HH

#include "globals.hh"
#include <assert.h>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4Transform
{
  public:  // with description

    G4Transform( G4RotationMatrix* pRot=0) :
	frotation(pRot) {;}

    G4Transform( G4RotationMatrix* pRot,
		 const G4ThreeVector& pVec) :
	frotation(pRot),
	ftranslation(pVec) {;}

    ~G4Transform() {;}

    G4ThreeVector& Apply(G4ThreeVector& pVec) const;
      // Apply to specified vector

    G4ThreeVector& ApplyReverse(G4ThreeVector& pVec) const;
      // Apply reverse transformation to specified vector

    G4ThreeVector Transform(const G4ThreeVector& pVec) const;
      // Compute the transform of the specified vector

    G4ThreeVector ReverseTransform(const G4ThreeVector& pVec) const;
      // Compute the transform of the specified vector

    G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const;

    G4RotationMatrix* GetRotation() const;
    const G4ThreeVector& GetTranslation() const;
    void SetTranslation(const G4ThreeVector& pVec);
    inline void SetRotation(G4RotationMatrix* pRot);

    G4Transform& ComputeCompoundTransform(const G4Transform& t1,
					  const G4Transform& t2,
					  G4RotationMatrix* pRot);

    G4Transform& ComputeCompoundReverseTransform(const G4Transform& t1,
					         const G4Transform& t2,
						 G4RotationMatrix* pRot);

    G4bool IsRotated() const;
      // Return true if transform involves rotations


    G4bool operator == (const G4Transform& t) const;
      // Define equality as
      // 1) equal rotation matrix *pointers*, and equal translations
      // or
      // 2) equal rotation matrices and equal translations

  private:

    G4RotationMatrix* frotation;
    G4ThreeVector ftranslation;
};

#include "G4Transform.icc"

#endif
