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
// $Id: G4TwistedBox.hh,v 1.5 2005/04/04 11:56:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedBox
//
// Class description:
//
//  A G4TwistedBox is a twisted cuboid of given half lengths pDx,pDy,pDz
//  and twist angle pPhiTwist. 
//  The Box is  centred on the origin with sides parallel to the x/y/z axes.
//
//   Member Data:
//
//     pDx    Half-length along x axis 
//     pDy    Half-length along y asis
//     pDz    Half-length along z axis
//     pPhiTwist Twist angle

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDBOX__
#define __G4TWISTEDBOX__

#include "G4VTwistedFaceted.hh"

class G4TwistedBox : public G4VTwistedFaceted
{
public:  // with description

  G4TwistedBox(const G4String& pName,
                     G4double  pPhiTwist,
                     G4double  pDx,
                     G4double  pDy, 
                     G4double  pDz );

  virtual ~G4TwistedBox();

  // accessors

  inline G4double GetXHalfLength() const { return GetDx1() ; }
  inline G4double GetYHalfLength() const { return GetDy1() ; }
  inline G4double GetZHalfLength() const { return GetDz()  ; }
  inline G4double GetPhiTwist()    const { return GetTwistAngle() ; }

  G4GeometryType GetEntityType()    const;
  G4Polyhedron*  CreatePolyhedron() const;

  std::ostream& StreamInfo(std::ostream& os) const;

};

#endif
