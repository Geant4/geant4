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
// class G4AffineTransform
//
// Class description:
//
// A class for geometric affine transformations [see, eg. Foley & Van Dam]
// Supports efficient arbitrary rotation & transformation of vectors and the
// computation of compound & inverse transformations. A 'rotation flag' is
// maintained internally for greater computational efficiency for transforms
// that do not involve rotation.
//
// Interfaces to the CLHEP classes G4ThreeVector & G4RotationMatrix
//
// For member function descriptions, see comments by declarations. For
// additional clarification, also check the `const' declarations for
// functions & their parameters.
//
// Member data:
//
//      G4double rxx,rxy,rxz; 
//      G4double ryx,ryy,ryz;  A 3x3 rotation matrix - net rotation
//      G4double rzx,rzy,rzz;
//      G4double tx,ty,tz;     Net translation 

// Author: Paul R C Kent (CERN), 06.08.1996 - Initial version
//         E.Tcherniaev (CERN), 19.09.1996, 06.05.2018 - Revised
// --------------------------------------------------------------------
#ifndef G4AFFINETRANSFORM_HH
#define G4AFFINETRANSFORM_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

/**
 * @brief G4AffineTransform is a class for geometric affine transformations.
 * It supports efficient arbitrary rotation & transformation of vectors and
 * the computation of compound & inverse transformations. A 'rotation flag'
 * is maintained internally for greater computational efficiency for transforms
 * that do not involve rotation.
 */

class G4AffineTransform
{
  public:

    /**
     * Constructor for G4AffineTransform. Initialises components to zero.
     */
    inline G4AffineTransform();

    /**
     * Constructor for Translation only: under t'form, translate point
     * at origin by 'tlate'.
     */
    inline G4AffineTransform(const G4ThreeVector& tlate);

    /**
     * Constructor for Rotation only: under t'form, rotate by 'rot'.
     */
    inline G4AffineTransform(const G4RotationMatrix& rot);

    /**
     * Constructor for Translation and Rotation: under t'form, rotate
     * by 'rot' then translate by 'tlate'.
     */
    inline G4AffineTransform(const G4RotationMatrix& rot,
                             const G4ThreeVector& tlate);

    /**
     * Alternative Constructor optionally rotating by 'rot' by pointer then
     * translate by 'tlate' - 'rot' may be null.
     */
    inline G4AffineTransform(const G4RotationMatrix* rot,
                             const G4ThreeVector& tlate);

    /**
     * Copy & move constructor.
     */
    inline G4AffineTransform(const G4AffineTransform& rhs) = default;
    inline G4AffineTransform(G4AffineTransform&& rhs) = default;

    /**
     * Assignment & move operators.
     */
    inline G4AffineTransform& operator=(const G4AffineTransform& rhs);
    inline G4AffineTransform& operator=(G4AffineTransform&& rhs) = default;

    /**
     * Default Destructor.
     */
    inline ~G4AffineTransform() = default;

    /**
     * Compound Transforms: tf2=tf2*tf1 equivalent to tf2*=tf1.
     *  @param[in] tf Transformation to combine.
     *  @returns The compound transformation of self*tf.
     */
    inline G4AffineTransform operator * (const G4AffineTransform& tf) const;

    /**
     * [Modifying] compound Transforms: Multiplies self by 'tf'.
     *  @param[in] tf Transformation to combine.
     *  @returns Returns self reference, i.e. A=AB for a*=b.
     */
    inline G4AffineTransform& operator *= (const G4AffineTransform& tf);

    /**
     * [Modifying] Product function, for avoiding (potential) temporaries:
     * c.Product(a,b) equivalent to c=a*b
     * c.InverseProduct(a*b,b ) equivalent to c=a
     * Sets self=tf1*tf2.
     *  @param[in] tf1 First transformation operand.
     *  @param[in] tf2 Second transformation operand.
     *  @returns Returns Self reference.
     */
    inline G4AffineTransform& Product(const G4AffineTransform& tf1,
                                      const G4AffineTransform& tf2);

    /**
     * [Modifying] Inverse Product function. Sets self=tf1*(tf2^-1).
     *  @param[in] tf1 First transformation operand.
     *  @param[in] tf2 Second transformation operand.
     *  @returns Returns Self reference.
     */
    inline G4AffineTransform& InverseProduct(const G4AffineTransform& tf1,
                                             const G4AffineTransform& tf2);

    /**
     * Transforms the specified point 'vec'.
     *  @returns vec*rot+tlate.
     */
    inline G4ThreeVector TransformPoint(const G4ThreeVector& vec) const;

    /**
     * Transforms the specified point 'vec' using inverse transformation.
     *  @returns The inverse transformation of the given point.
     */
    inline G4ThreeVector InverseTransformPoint(const G4ThreeVector& vec) const;

    /**
     * Transforms the specified 'axis'.
     *  @returns vec*rot.
     */
    inline G4ThreeVector TransformAxis(const G4ThreeVector& axis) const;

    /**
     * Transforms the specified 'axis' using inverse transformation.
     *  @returns The inverse transformation of the given axis.
     */
    inline G4ThreeVector InverseTransformAxis(const G4ThreeVector& axis) const;

    /**
     * Transforms the specified point 'vec' (in place): sets vec=vec*rot+tlate.
     *  @param[in,out] vec The point to transform.
     */
    inline void ApplyPointTransform(G4ThreeVector& vec) const;

    /**
     * Transforms the specified 'axis' (in place): sets axis=axis*rot.
     *  @param[in,out] axis The axis to transform.
     */
    inline void ApplyAxisTransform(G4ThreeVector& axis) const;

    /**
     * Returns the inverse of the current transform.
     */
    inline G4AffineTransform Inverse() const;

    /**
     * [Modifying] Sets self=inverse of self.
     *  @returns Self reference.
     */
    inline G4AffineTransform& Invert();

    /**
     * [Modifying] Adjust the net translation by the given vector.
     *  @returns Self reference.
     */
    inline G4AffineTransform& operator +=(const G4ThreeVector& tlate);
    inline G4AffineTransform& operator -=(const G4ThreeVector& tlate);

    /**
     * Equality and inequality operators.
     */
    inline G4bool operator == (const G4AffineTransform& tf) const;
    inline G4bool operator != (const G4AffineTransform& tf) const;

    /**
     * Access operator.
     */
    inline G4double operator [] (const G4int n) const;

    /**
     * Returns true if transform includes rotation.
     */
    inline G4bool IsRotated() const;

    /**
     * Returns true if transform includes translation.
     */
    inline G4bool IsTranslated() const;

    /**
     * Returns the net rotation matrix.
     */
    inline G4RotationMatrix NetRotation() const;

    /**
     * Returns the inverse net rotation matrix.
     */
    inline G4RotationMatrix InverseNetRotation() const;

    /**
     * Returns the net translation vector.
     */
    inline G4ThreeVector NetTranslation() const;

    /**
     * Returns the inverse net translation vector.
     */
    inline G4ThreeVector InverseNetTranslation() const;

    /**
     * Setters for rotation and translation.
     */
    inline void SetNetRotation(const G4RotationMatrix& rot);
    inline void SetNetTranslation(const G4ThreeVector& tlate);

    /**
     * Conversion operator (cast) to G4Transform3D.
     */
    inline operator G4Transform3D () const;

  private:

    /**
     * Private components Constructor.
     */
    inline G4AffineTransform(
                  const G4double prxx, const G4double prxy, const G4double prxz,
                  const G4double pryx, const G4double pryy, const G4double pryz,
                  const G4double przx, const G4double przy, const G4double przz,
                  const G4double ptx, const G4double pty, const G4double ptz);

  private:

    G4double rxx,rxy,rxz;
    G4double ryx,ryy,ryz;
    G4double rzx,rzy,rzz;
    G4double tx,ty,tz;
};

/**
 * Streaming operator.
 */
std::ostream& operator << (std::ostream& os, const G4AffineTransform& transf);

#include "G4AffineTransform.icc"

#endif
