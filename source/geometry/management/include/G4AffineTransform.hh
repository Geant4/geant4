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
// $Id: G4AffineTransform.hh 98309 2016-07-06 10:30:15Z gcosmo $
//
//
// class G4AffineTransform
//
// Class description:
//
// A class for geometric affine transformations [see, eg. Foley & Van Dam]
// Supports efficient arbitrary rotation & transformation of vectors and the
// computation of compound & inverse transformations. A `rotation flag' is
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

// History:
// Paul R C Kent 6 Aug 1996 - initial version
//
// 19.09.96 E.Chernyaev:
// - direct access to the protected members of the G4RotationMatrix class
//   replaced by access via public access functions            
// - conversion of the rotation matrix to angle & axis used to get
//   a possibility to remove "friend" from the G4RotationMatrix class
// --------------------------------------------------------------------
#ifndef G4AFFINETRANSFORM_HH
#define G4AFFINETRANSFORM_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

class G4AffineTransform
{

public:

  inline G4AffineTransform();

public: // with description

  inline G4AffineTransform(const G4ThreeVector& tlate);
    // Translation only: under t'form translate point at origin by tlate

  inline G4AffineTransform(const G4RotationMatrix& rot);
    // Rotation only: under t'form rotate by rot

  inline G4AffineTransform(const G4RotationMatrix& rot,
                           const G4ThreeVector& tlate);
    // Under t'form: rotate by rot then translate by tlate

  inline G4AffineTransform(const G4RotationMatrix* rot,
                           const G4ThreeVector& tlate);
    // Optionally rotate by *rot then translate by tlate - rot may be null
  
  inline G4AffineTransform operator * (const G4AffineTransform& tf) const;
    // Compound Transforms:
    //       tf2=tf2*tf1 equivalent to tf2*=tf1
    //       Returns compound transformation of self*tf

  inline G4AffineTransform& operator *= (const G4AffineTransform& tf);
    // (Modifying) Multiplies self by tf; Returns self reference
    //             ie. A=AB for a*=b


  inline G4AffineTransform& Product(const G4AffineTransform& tf1,
                                    const G4AffineTransform& tf2);
    // 'Products' for avoiding (potential) temporaries:
    //            c.Product(a,b) equivalent to c=a*b
    //            c.InverseProduct(a*b,b ) equivalent to c=a
    // (Modifying) Sets self=tf1*tf2; Returns self reference

  inline G4AffineTransform& InverseProduct(const G4AffineTransform& tf1,
                                           const G4AffineTransform& tf2);
    // (Modifying) Sets self=tf1*(tf2^-1); Returns self reference

  inline G4ThreeVector TransformPoint(const G4ThreeVector& vec) const;
    // Transform the specified point: returns vec*rot+tlate

  inline G4ThreeVector TransformAxis(const G4ThreeVector& axis) const;
    // Transform the specified axis: returns

  inline void ApplyPointTransform(G4ThreeVector& vec) const;
    // Transform the specified point (in place): sets vec=vec*rot+tlate

  inline void ApplyAxisTransform(G4ThreeVector& axis) const;
    // Transform the specified axis (in place): sets axis=axis*rot;

  inline G4AffineTransform Inverse() const;
    // Return inverse of current transform

  inline G4AffineTransform& Invert();
    // (Modifying) Sets self=inverse of self; Returns self reference

  inline G4AffineTransform& operator +=(const G4ThreeVector& tlate);
  inline G4AffineTransform& operator -=(const G4ThreeVector& tlate);
    // (Modifying) Adjust net translation by given vector;
    //             Returns self reference

  inline G4bool operator == (const G4AffineTransform& tf) const;
  inline G4bool operator != (const G4AffineTransform& tf) const;

  inline G4double operator [] (const G4int n) const;

  inline G4bool IsRotated() const;
    // True if transform includes rotation

  inline G4bool IsTranslated() const;
    // True if transform includes translation

  inline G4RotationMatrix NetRotation() const;

  inline G4ThreeVector NetTranslation() const;

  inline void SetNetRotation(const G4RotationMatrix& rot);

  inline void SetNetTranslation(const G4ThreeVector& tlate);

  inline operator G4Transform3D () const;
    // Conversion operator (cast) to G4Transform3D

private:

  inline G4AffineTransform(
                  const G4double prxx, const G4double prxy, const G4double prxz,
                  const G4double pryx, const G4double pryy, const G4double pryz,
                  const G4double przx, const G4double przy, const G4double przz,
                  const G4double ptx, const G4double pty, const G4double ptz);

  G4double rxx,rxy,rxz;
  G4double ryx,ryy,ryz;
  G4double rzx,rzy,rzz;
  G4double tx,ty,tz;
};

std::ostream& operator << (std::ostream& os, const G4AffineTransform& transf);

#include "G4AffineTransform.icc"

#endif
