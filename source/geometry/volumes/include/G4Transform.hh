// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Transform.hh,v 1.2 1999-12-15 14:50:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Transform
//
// A coordinate transform, consisting of a translation and a rotation.
// When applying the forwards transformations, rotation is applied first, then
// the translation.
//
//
// Member functions:
//
// G4Transform(G4RotationMatrix* pRot=0)
//   Create a translform consisting of a rotation by the specified matrix.
//   If NULL (the default) and indentity transform is constructed
// G4Transform(G4RotationMatrix* pRot,const G4ThreeVector& pVec)
//   Create a transform consisting of a rotation by the specified matrix, then
//   a translation by specified vector. The transform does not copy the
//   specified matrix: ownership (hence deletion responsibility) remains with
//   the caller.
//
//   If the rotation is 0 (NULL ptr), no rotation is performed.
//
// G4ThreeVector& Apply(G4ThreeVector& pVec)
//   Apply to the transformation to the specified vector.
// G4ThreeVector& ApplyReverse(G4ThreeVector& pVec)
//   Apply to the reverse transformation to the specified vector.
//
// G4ThreeVector Transform(const G4ThreeVector& pVec)
//   Compute the transformation of the specified vector.
// G4ThreeVector ReverseTransform(const G4ThreeVector& pVec)
//   Compute the reverse transformation to the specified vector.
//
// G4RotationMatrix* GetRotation()
//   Return ptr to rotation for the transformation (possibly null if no
//   rotation)
// const G4ThreeVector& GetTranslation() const
//   Return the translation component of the transformation
// G4Transform& SetTranslation(const G4ThreeVector& pVec)
//   Set the translation resulting from the transform.
//   Returns a self reference.
//
//   [Intended so that simple modifications can be made without creating
//    new transforms]
//
// G4Transform& SetRotation(G4RotationMatrix* pRot)
//   Set the translation resulting from the transform.
//   Returns a self reference.
//
// G4Transform& ComputeCompoundTransform(const G4Transform& t1,
//                         const G4Transform& t2,G4RotationMatrix* pRot)
//   Compute the compound transformation resulting from the transform t1
//   followed by the transform t2. The Rotation matrix pRot is modified
//   and set to the rotation resulting from the transformation. If no
//   rotation results, the rotation ptr of the compounded transform is set
//   to NULL. Returns self reference.
//
//   If transforms are considered as 4*4 matrices computes t1 . t2
//
// G4Transform& ComputeCompoundReverseTransform(const G4Transform& t1,
//                         const G4Transform& t2,G4RotationMatrix* pRot)
//   Compute the compound transformation resulting from the inverse of 
//   transform t1 followed by the transform t2. The Rotation matrix pRot
//   is modified and set to the rotation resulting from the transformation.
//   If no rotation results, the rotation ptr of the compounded transform
//   is set to NULL. Returns self reference.
//
//   If transforms are considered as 4*4 matrices computes [t1^-1].t2
//
// G4bool IsRotated() const
//   Return true if transform involves rotation
//
// G4bool operator == (const G4Transform& t) const
//   Define equality as
//    1) equal translations and equal rotation matrix *pointers*
//    or
//    2) equal translation and rotation matrices (contents)
//
// Member data:
//
// G4RotationMatrix* frotation
//   Pointer to the transmations rotation component
// G4ThreeVector ftranslation
//   Translation component
//
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
public:
    inline G4Transform( G4RotationMatrix* pRot=0) :
	frotation(pRot)
    {;}

    inline G4Transform( G4RotationMatrix* pRot,
		       const G4ThreeVector& pVec) :
	frotation(pRot),
	ftranslation(pVec)
    {;}

//    inline ~G4Transform() {;}

// Apply to specified vector
    inline G4ThreeVector& Apply(G4ThreeVector& pVec) const
    {
	if (frotation)
	    {
		return pVec=frotation->operator*(pVec)+ftranslation;
	    }
	else
	    {
		return pVec+=ftranslation;
	    }
    }

// Apply reverse transformation to specified vector
    inline G4ThreeVector& ApplyReverse(G4ThreeVector& pVec) const
    {
	if (frotation)
	    {
		return pVec=frotation->inverse()*(pVec-ftranslation);
	    }
	else
	    {
		return pVec-=ftranslation;
	    }
    }

// Compute the transform of the specified vector
    inline G4ThreeVector Transform(const G4ThreeVector& pVec) const
    {
	if (frotation)
	    {
		return frotation->operator*(pVec)+ftranslation;
	    }
	else
	    {
		return pVec+ftranslation;
	    }
    }

// Compute the transform of the specified vector
    inline G4ThreeVector ReverseTransform(const G4ThreeVector& pVec) const
    {
	if (frotation)
	    {
		return frotation->inverse()*(pVec-ftranslation);
	    }
	else
	    {
		return pVec-ftranslation;
	    }
    }


    inline G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const
    {
	if (!frotation)
	    {
		return pVec;
	    }
	else
	    {
		return  (*frotation)*pVec; // UNTESTED
	    }
    }

    inline G4RotationMatrix* GetRotation() const
    {
	return frotation;
    }

    inline const G4ThreeVector& GetTranslation() const
    {
	return ftranslation;
    }

    inline void SetTranslation(const G4ThreeVector& pVec)
    {
	ftranslation=pVec;
    }

    inline void SetRotation(G4RotationMatrix* pRot)
    {
	frotation=pRot;
    }

    inline G4Transform& ComputeCompoundTransform(const G4Transform& t1,
					  const G4Transform& t2,
					  G4RotationMatrix* pRot)
    {
	if (!(t1.GetRotation()&&t2.GetRotation()))
	    {
		SetTranslation(t1.GetTranslation()+t2.GetTranslation());
		SetRotation(0);
	    }
	else
	    {
	      pRot=pRot;	// Avoid unuse warning
		assert(0==1);	// Unimpl 13.7.95
	    }
	return *this;
    }

    inline G4Transform& ComputeCompoundReverseTransform(const G4Transform& t1,
							const G4Transform& t2,
							G4RotationMatrix* pRot)
    {
	if (!(t1.GetRotation()&&t2.GetRotation()))
	    {
		SetTranslation(t2.GetTranslation()-t1.GetTranslation());
		SetRotation(0);
	    }
	else
	    {
	      pRot=pRot;	// Avoid unused warning
		assert(0==1);	// Unimpl 13.7.95
	    }
	return *this;
    }

// Return true if transform involves rotations
    inline G4bool IsRotated() const
    {
	return (!frotation||frotation->isIdentity()) ? false : true;
    }

// Define equality as
// 1) equal rotation matrix *pointers*, and equal translations
// or
// 2) equal rotation matrices and equal translations

    inline G4bool operator == (const G4Transform& t) const
    {
	return ((GetTranslation()==t.GetTranslation())&&
		( (GetRotation()==t.GetRotation())
		  || (*GetRotation()==*t.GetRotation()) )) ? true : false;
    }

private:
    G4RotationMatrix* frotation;
    G4ThreeVector ftranslation;
};


#endif



