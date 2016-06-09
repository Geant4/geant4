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
// $Id: G4TwistedTrap.hh,v 1.5 2005/04/04 11:56:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrap
//
// Class description:
//
//  A G4TwistedTrap is a general twisted trapezoid: The faces perpendicular to the
//  z planes are trapezia, and their centres are not necessarily on
//  a line parallel to the z axis.
//
//      pDz     Half-length along the z-axis
//      pTheta  Polar angle of the line joining the centres of the faces
//              at -/+pDz
//      pPhi    Azimuthal angle of the line joing the centre of the face at
//              -pDz to the centre of the face at +pDz
//      pDy1    Half-length along y of the face at -pDz
//      pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
//      pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
//
//      pDy2    Half-length along y of the face at +pDz
//      pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
//      pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
//      pAlph   Angle with respect to the y axis from the centre of the side
//
//
//  A special regular case of a trapezoid with equal endcaps is available, 
//  with polar,azimuthal and tilt angles set to zero.
//
 
// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRAP__
#define __G4TWISTEDTRAP__

#include "G4VTwistedFaceted.hh"

class G4TwistedTrap : public G4VTwistedFaceted
{
public:  // with description

  G4TwistedTrap(const G4String &pName,
                      G4double  pPhiTwist,
                      G4double  pDx1,  // half x length at -pDz,-pDy
                      G4double  pDx2,  // half x length at -pDz,+pDy
                      G4double  pDy,
                      G4double  pDz);


  G4TwistedTrap(const G4String &pName,       // Name of instance
                      G4double  pPhiTwist,   // twist angle
                      G4double  pDz,     // half z length
                      G4double  pTheta,  // direction between end planes
                      G4double  pPhi,    // defined by polar and azim. angles
                      G4double  pDy1,    // half y length at -pDz
                      G4double  pDx1,    // half x length at -pDz,-pDy
                      G4double  pDx2,    // half x length at -pDz,+pDy
                      G4double  pDy2,    // half y length at +pDz
                      G4double  pDx3,    // half x length at +pDz,-pDy
                      G4double  pDx4,    // half x length at +pDz,+pDy
                      G4double  pAlph    // tilt angle
                );
 
  virtual ~G4TwistedTrap();

  // accessors

  inline G4double GetY1HalfLength() const { return GetDy1() ; }
  inline G4double GetX1HalfLength() const { return GetDx1() ; }
  inline G4double GetX2HalfLength() const { return GetDx2() ; }
  inline G4double GetY2HalfLength() const { return GetDy2() ; }
  inline G4double GetX3HalfLength() const { return GetDx3() ; }
  inline G4double GetX4HalfLength() const { return GetDx4() ; }
  inline G4double GetZHalfLength()  const { return GetDz()  ; }
  inline G4double GetPhiTwist()     const { return GetTwistAngle() ; }
  inline G4double GetTiltAngleAlpha()    const { return GetAlpha() ; }
  inline G4double GetPolarAngleTheta()   const { return GetTheta() ; }
  inline G4double GetAzimuthalAnglePhi() const { return GetPhi()   ; }

  G4GeometryType GetEntityType() const;
  G4Polyhedron*  CreatePolyhedron() const;

  std::ostream &StreamInfo(std::ostream& os) const;

} ;

#endif
