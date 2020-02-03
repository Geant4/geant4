// -*- C++ -*-
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepBoostZ class for performing specialized
// Lorentz transformations which are pure boosts in the Z direction, on 
// objects of the HepLorentzVector class.
//
// HepLorentzRotation is a concrete implementation of Hep4RotationInterface.
//
// .SS See Also
// RotationInterfaces.h
// LorentzVector.h LorentzRotation.h 
// Boost.h 
//
// .SS Author
// Mark Fischler

#ifndef HEP_BOOSTZ_H
#define HEP_BOOSTZ_H

#ifdef GNUPRAGMA
#pragma interface
#endif

#include "CLHEP/Vector/RotationInterfaces.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace CLHEP  {

// Declarations of classes and global methods
class HepBoostZ;               
inline HepBoostZ inverseOf ( const HepBoostZ & b ); 
class HepBoost;
class HepRotation;

/**
 * @author
 * @ingroup vector
 */
class HepBoostZ {

public:

  // ----------  Constructors and Assignment:

  inline HepBoostZ();
  // Default constructor. Gives a boost of 0.  

  inline HepBoostZ(const HepBoostZ & b);
  inline HepBoostZ(HepBoostZ && b) = default;
  // Copy and move constructors.

  inline HepBoostZ & operator = (const HepBoostZ & m);
  inline HepBoostZ & operator = (HepBoostZ && m) = default;
  // Copy and move assignment operators.

         HepBoostZ & set (double beta);
  inline HepBoostZ       (double beta);
  // Constructor from beta 

  // ----------  Accessors:

  inline double  beta()  const;
  inline double  gamma() const;
  inline Hep3Vector boostVector() const;
  inline Hep3Vector getDirection() const;

  inline double xx() const;
  inline double xy() const;
  inline double xz() const;
  inline double xt() const;
  inline double yx() const;
  inline double yy() const;
  inline double yz() const;
  inline double yt() const;
  inline double zx() const;
  inline double zy() const;
  inline double zz() const;
  inline double zt() const;
  inline double tx() const;
  inline double ty() const;
  inline double tz() const;
  inline double tt() const;
  // Elements of the matrix.

  inline HepLorentzVector col1() const;
  inline HepLorentzVector col2() const;
  inline HepLorentzVector col3() const;
  inline HepLorentzVector col4() const;
  // orthosymplectic column vectors

  inline HepLorentzVector row1() const;
  inline HepLorentzVector row2() const;
  inline HepLorentzVector row3() const;
  inline HepLorentzVector row4() const;
  // orthosymplectic row vectors

  HepRep4x4 rep4x4() const;
  //   4x4 representation:

  HepRep4x4Symmetric rep4x4Symmetric() const;
  //   Symmetric 4x4 representation.

  // ----------  Decomposition:

  void decompose (HepRotation  & rotation, HepBoost   & boost) const;
  void decompose (HepAxisAngle & rotation, Hep3Vector & boost) const;
  // Find R and B such that L = R*B -- trivial, since R is identity

  void decompose (HepBoost   & boost, HepRotation  & rotation) const;
  void decompose (Hep3Vector & boost, HepAxisAngle & rotation) const;
  // Find R and B such that L = B*R -- trivial, since R is identity

  // ----------  Comparisons:

  inline int compare( const HepBoostZ & b  ) const;
  // Dictionary-order comparison, in order of beta. 
  // Used in operator<, >, <=, >=

  inline bool operator == (const HepBoostZ & b) const;
  inline bool operator != (const HepBoostZ & b) const;
  inline bool operator <= (const HepBoostZ & b) const;
  inline bool operator >= (const HepBoostZ & b) const;
  inline bool operator <  (const HepBoostZ & b) const;
  inline bool operator >  (const HepBoostZ & b) const;
  // Comparisons.

  inline bool isIdentity() const;
  // Returns true if a null boost.

  inline  double distance2( const HepBoostZ & b ) const;
  	  double distance2( const HepBoost & b ) const;
  // Defined as the distance2 between the vectors (gamma*betaVector)

  double distance2( const HepRotation & r  ) const;
  double distance2( const HepLorentzRotation & lt  ) const;
  // Decompose lt = B*R; add norm2 to distance2 to between boosts.

  inline double howNear(   const HepBoostZ & b ) const;
  inline double howNear(   const HepBoost  & b ) const;
  inline double howNear(   const HepRotation & r ) const;
  inline double howNear(   const HepLorentzRotation & lt ) const;

  inline bool isNear(   const HepBoostZ & b,
             double epsilon=Hep4RotationInterface::tolerance) const;
  inline bool isNear(   const HepBoost & b,
             double epsilon=Hep4RotationInterface::tolerance) const;
  bool isNear(   const HepRotation & r,
             double epsilon=Hep4RotationInterface::tolerance) const;
  bool isNear(   const HepLorentzRotation & lt,
             double epsilon=Hep4RotationInterface::tolerance) const;

  // ----------  Properties:

  inline double norm2() const;
  // distance2 (IDENTITY), which is beta^2 * gamma^2

  void rectify();
  // sets according to the stored beta

  // ---------- Application:

  inline HepLorentzVector operator()( const HepLorentzVector & w ) const;
  // Transform a Lorentz Vector.             

  inline HepLorentzVector operator* ( const HepLorentzVector & w ) const;
  // Multiplication with a Lorentz Vector.

  // ---------- Operations in the group of 4-Rotations

  HepBoostZ   operator * (const HepBoostZ & b) const;
  HepLorentzRotation operator * (const HepBoost & b) const;
  HepLorentzRotation operator * (const HepRotation & r) const;
  HepLorentzRotation operator * (const HepLorentzRotation & lt) const;
  // Product of two Lorentz Rotations (this) * lt - matrix multiplication
  // Notice that the product of two pure boosts in different directions
  // is no longer a pure boost.

  inline HepBoostZ inverse() const;
  // Return the inverse.

  inline friend HepBoostZ inverseOf ( const HepBoostZ & b );
  // global methods to invert.

  inline HepBoostZ & invert();
  // Inverts the Boost matrix.

  // ---------- I/O:

  std::ostream & print( std::ostream & os ) const;
  // Output form is BOOSTZ (beta=..., gamma=...);  

  // ---------- Tolerance

  static inline double getTolerance();
  static inline double setTolerance(double tol);

protected:

  inline HepLorentzVector vectorMultiplication
                                        ( const HepLorentzVector & w ) const;
  // Multiplication with a Lorentz Vector.

  HepLorentzRotation matrixMultiplication (const HepRep4x4 & m) const;
  HepLorentzRotation matrixMultiplication (const HepRep4x4Symmetric & m) const;
                   
  inline HepBoostZ (double beta, double gamma);

  double  beta_;
  double  gamma_;

};  // HepBoostZ

inline   
std::ostream & operator << 
	( std::ostream & os, const HepBoostZ& b ) {return b.print(os);}

}  // namespace CLHEP

#include "CLHEP/Vector/BoostZ.icc"

#endif /* HEP_BOOSTZ_H */
