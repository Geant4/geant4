// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AffineTransform.hh,v 1.2 1999-12-15 14:49:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4AffineTransform Header File
//
// A class for geometric affine transformations [see, eg. Foley & Van Dam]
// Supports efficient arbitrary rotation & transformation of vectors and the
// computation of compound & inverse transformations. A `rotation flag' is
// maintained internally for greater computational efficiency for transforms
// that do not involve rotation.
//
// Interfaces to the GEANT4 modified CLHEP classes G4ThreeVector &
// G4RotationMatrix
//
// For member function descriptions, see comments by declarations. For
// additional clarification, also check the `const' declarations for
// functions & their parameters.
//
// Private Member data:
//
//	G4double rxx,rxy,rxz; 
//	G4double ryx,ryy,ryz;  A 3x3 rotation matrix - net rotation
//	G4double rzx,rzy,rzz;
//	G4double tx,ty,tz;     Net translation 
//
// History:
// Paul R C Kent 6 Aug 1996 - initial version
//
// 19.09.96 E.Chernyaev:
// - direct access to the protected members of the G4RotationMatrix class
//   replaced by access via public access functions            
// - conversion of the rotation matrix to angle & axis used to get
//   a possibility to remove "friend" from the G4RotationMatrix class

#ifndef G4AFFINETRANSFORM_HH
#define G4AFFINETRANSFORM_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4AffineTransform
{
public:
	G4AffineTransform();

// Translation only: Under t'form translate point at origin by tlate

        G4AffineTransform(const G4ThreeVector &tlate);
 
// Rotation only: Under t'form rotate by rot

	G4AffineTransform(const G4RotationMatrix &rot);
	
// Under t'form: Rotate by rot then translate by tlate

	G4AffineTransform(const G4RotationMatrix &rot,	
			  const G4ThreeVector &tlate);
			  
// Optionally rotate by *rot then translate by tlate - rot may be null

        G4AffineTransform(const G4RotationMatrix *rot,
			  const G4ThreeVector &tlate);
			  
			  
// NOTE: Compound Transforms
// 
// tf2=tf2*tf1 equivalent to tf2*=tf1
//
// Returns compound transformation of self*tf

        G4AffineTransform operator * (const G4AffineTransform &tf) const;
	
// (Modifying) Multiplies self by tf; Returns self reference
//             ie. A=AB for a*=b

        G4AffineTransform& operator *= (const G4AffineTransform &tf);
	

// 'Products' for avoiding (potential) temporaries
//
// c.Product(a,b) equivalent to c=a*b
//
// c.InverseProduct(a*b,b ) equivalent to c=a
//
// (Modifying) Sets self=tf1*tf2; Returns self reference

	G4AffineTransform& Product(const G4AffineTransform &tf1,
				   const G4AffineTransform &tf2);
				   
// (Modifying) Sets self=tf1*(tf2^-1); Returns self reference

        G4AffineTransform& InverseProduct(const G4AffineTransform &tf1,
				          const G4AffineTransform &tf2);
	
// Transform the specified point: returns vec*rot+tlate

	G4ThreeVector TransformPoint(const G4ThreeVector &vec) const;
	
// Transform the specified axis: returns

	G4ThreeVector TransformAxis(const G4ThreeVector &axis) const;
	
// Transform the specified point (in place): sets vec=vec*rot+tlate

	void ApplyPointTransform(G4ThreeVector &vec) const;
	
// Transform the specified axis (in place): sets axis=axis*rot;

	void ApplyAxisTransform(G4ThreeVector &axis) const;

// Return inverse of current transform

	G4AffineTransform Inverse() const;

// (Modifying) Sets self=inverse of self; Returns self reference

	G4AffineTransform& Invert();

// (Modifying) Adjust net translation by given vector; Returns self reference

	G4AffineTransform& operator +=(const G4ThreeVector &tlate);
	G4AffineTransform& operator -=(const G4ThreeVector &tlate);
	
	G4bool operator == (const G4AffineTransform &tf) const;
	G4bool operator != (const G4AffineTransform &tf) const;
	
	G4double operator [] (const G4int n) const;

// True if transform includes rotation

	G4bool IsRotated() const;
	
// Ture if transform includes translation

	G4bool IsTranslated() const;

	G4RotationMatrix NetRotation() const;
	
	G4ThreeVector NetTranslation() const;
	
	void SetNetRotation(const G4RotationMatrix &rot);
	
	void SetNetTranslation(const G4ThreeVector &tlate);

private:
	G4AffineTransform( const G4double prxx,const G4double prxy,const G4double prxz,
	          	   const G4double pryx,const G4double pryy,const G4double pryz,
		           const G4double przx,const G4double przy,const G4double przz,
			   
	 	           const G4double ptx,const G4double pty,const G4double ptz     ) ;

	G4double rxx,rxy,rxz;
	G4double ryx,ryy,ryz;
	G4double rzx,rzy,rzz;
	
	G4double tx,ty,tz;
};

#include "G4AffineTransform.icc"

#endif




