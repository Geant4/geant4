//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LineSection.hh,v 1.5 2001-07-11 09:59:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4LineSection
//
// Class description:
//
// A utility class that calculates the distance of a point from a 
// line section.

// History:
// - Created. J. Apostolakis.
// - Cosmetics against /0. and sqrt(<0.). V. Grichine, 25.02.00.

#ifndef G4LineSection_hh
#define G4LineSection_hh

#include "globals.hh" 
#include "G4ThreeVector.hh"

class G4LineSection
{
  public:  // with description

     G4LineSection( const G4ThreeVector& PntA, const G4ThreeVector& PntB );

     G4double Dist( G4ThreeVector OtherPnt ) const;
     G4double InvsqDistAB() const;

     static G4double Distline( const G4ThreeVector& OtherPnt, 
			       const G4ThreeVector& LinePntA, 
			       const G4ThreeVector& LinePntB );
  private:

     G4ThreeVector    EndpointA;
     G4ThreeVector   VecAtoB;
     G4double inverse_square_distAB;
};

inline
G4double G4LineSection::InvsqDistAB() const
{
  return inverse_square_distAB;
}

#endif
