// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LineSection.hh,v 1.1 1999-01-07 16:07:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  A utility class that calculates the distance of a point from a 
//   line section.
//
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
