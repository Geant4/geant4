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
// $Id: G4TwistedTrd.hh,v 1.2 2005/04/04 11:56:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrd
//
// Class description:
//
//  A G4TwistedTrd is a twisted trapezoid with the x and y dimensions 
//  varying along z 
// 
//
//   Member Data:
//
//     pDx1    Half-length along x at the surface positioned at -dz
//     pDx2    Half-length along x at the surface positioned at +dz
//     pDy1    Half-length along y at the surface positioned at -dz
//     pDy2    Half-length along y at the surface positioned at +dz
//     pDz     Half-length along z axis
//     pPhiTwist Twist angle

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRD__
#define __G4TWISTEDTRD__

#include "G4VTwistedFaceted.hh"

class G4TwistedTrd : public G4VTwistedFaceted
{
public:  // with description

  G4TwistedTrd( const G4String& pName,
                      G4double  pDx1,
                      G4double  pDx2,
                      G4double  pDy1, 
                      G4double  pDy2,
                      G4double  pDz,
                      G4double  pPhiTwist );

  virtual ~G4TwistedTrd();

  // accessors

  inline G4double GetX1HalfLength() const { return GetDx1() ; }
  inline G4double GetX2HalfLength() const { return GetDx3() ; }
  inline G4double GetY1HalfLength() const { return GetDy1() ; }
  inline G4double GetY2HalfLength() const { return GetDy2() ; }
  inline G4double GetZHalfLength()  const { return GetDz()  ; }
  inline G4double GetPhiTwist()     const { return GetTwistAngle() ; }

  G4GeometryType GetEntityType() const;
  G4Polyhedron*  CreatePolyhedron   () const;

  std::ostream&  StreamInfo(std::ostream& os) const;

} ;

#endif
