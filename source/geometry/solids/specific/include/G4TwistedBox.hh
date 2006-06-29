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
// $Id: G4TwistedBox.hh,v 1.10 2006-06-29 18:48:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//   27-Oct-2004 - O.Link (Oliver.Link@cern.ch)
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

    std::ostream& StreamInfo(std::ostream& os) const;

  public:  // without description

    G4TwistedBox(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.
};

#endif
