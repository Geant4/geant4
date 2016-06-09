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
// $Id: G4GeomTestPoint.hh,v 1.2 2003/11/03 17:15:20 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestPoint
//
// Class description:
//
// A point on a line segment intersecting a solid, used
// for geometry tests

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestPoint_hh
#define G4GeomTestPoint_hh

#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4VSolid;

class G4GeomTestPoint
{
  public:  // with description

    G4GeomTestPoint();
    G4GeomTestPoint( const G4GeomTestPoint &other );
    G4GeomTestPoint( const G4ThreeVector &thePoint,
                           G4double theS,
                           G4bool isEntering );
    virtual ~G4GeomTestPoint();
      // Constructors and virtual destructor

    G4bool operator==( const G4GeomTestPoint &other ) const;
    G4bool operator< ( const G4GeomTestPoint &other ) const;
    G4bool operator<=( const G4GeomTestPoint &other ) const;
      // Operators

    virtual const G4ThreeVector &GetPosition() const;
    virtual G4double GetDistance() const;
    virtual G4bool Entering() const;
      // Accessors
  
  protected:

    G4ThreeVector p;
    G4double s;      
    G4bool entering;
};

#endif
