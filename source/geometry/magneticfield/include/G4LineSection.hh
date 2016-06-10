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
// $Id: G4LineSection.hh 85845 2014-11-05 15:43:58Z gcosmo $
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
// --------------------------------------------------------------------

#ifndef G4LineSection_hh
#define G4LineSection_hh

#include "G4Types.hh" 
#include "G4ThreeVector.hh"

class G4LineSection
{
  public:  // with description

     inline G4LineSection( const G4ThreeVector& PntA, const G4ThreeVector& PntB );

     G4double Dist( G4ThreeVector OtherPnt ) const;

     inline G4double GetABdistanceSq() const;

     inline static G4double Distline( const G4ThreeVector& OtherPnt, 
                                      const G4ThreeVector& LinePntA, 
                                      const G4ThreeVector& LinePntB );
  private:

     G4ThreeVector   EndpointA;
     G4ThreeVector   VecAtoB;
     G4double fABdistanceSq ;
};

// Inline methods implementations

inline
G4LineSection::G4LineSection( const G4ThreeVector& PntA, 
			      const G4ThreeVector& PntB )
  : EndpointA(PntA), VecAtoB(PntB-PntA)
{ 
  fABdistanceSq = VecAtoB.mag2();  
}

inline
G4double G4LineSection::GetABdistanceSq() const
{
  return fABdistanceSq;
}

inline
G4double G4LineSection::Distline( const G4ThreeVector& OtherPnt, 
				  const G4ThreeVector& LinePntA, 
				  const G4ThreeVector& LinePntB )
{
  G4LineSection LineAB( LinePntA, LinePntB );  // Line from A to B
  return LineAB.Dist( OtherPnt );
}

#endif
