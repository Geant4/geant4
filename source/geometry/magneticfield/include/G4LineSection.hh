// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LineSection.hh,v 1.3 2000-02-25 16:57:02 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  A utility class that calculates the distance of a point from a 
//   line section.
//
// History:
// ~1996 J. Apostolakis, first version
// 25.02.00 V. Grichine, cosmetics against /0. and sqrt(<0.)
//  

#include "globals.hh" 
#include "G4ThreeVector.hh"

typedef G4ThreeVector  POINT;
typedef POINT          Vector;


class G4LineSection {
  public:
     G4LineSection( const POINT& PntA, const POINT& PntB );

     G4double Dist( POINT OtherPnt ) const;
     G4double InvsqDistAB() const { return inverse_square_distAB; }  

     //
     static G4double Distline( const POINT& OtherPnt, 
			       const POINT& LinePntA, 
			       const POINT& LinePntB );
  private:
     POINT    EndpointA;
     Vector   VecAtoB;
     G4double inverse_square_distAB;
};
