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
// $Id: G4GeomTestVolPoint.hh,v 1.1 2001/10/17 12:59:56 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestVolPoint
//
// Class description:
//
// A point on a line segment intersecting a physical solid, used
// for geometry tests.
// A negative daughter index indicates the mother volume itself.

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)

#ifndef G4GeomTestVolPoint_hh
#define G4GeomTestVolPoint_hh

#include "G4GeomTestPoint.hh"
#include "G4RotationMatrix.hh"

class G4GeomTestVolPoint : public G4GeomTestPoint
{
  public:  // with description

    G4GeomTestVolPoint();
    G4GeomTestVolPoint( const G4ThreeVector &thePoint,
                              G4double theS,
                              G4bool isEntering,
                              G4int theDaughterIndex );
    G4GeomTestVolPoint( const G4GeomTestPoint &base,
                              G4int theDaughterIndex );
    G4GeomTestVolPoint( const G4GeomTestPoint &base,
                              G4int theDaughterIndex,
                        const G4ThreeVector &translation,
                        const G4RotationMatrix *rotation=0 );
    G4GeomTestVolPoint( const G4GeomTestVolPoint &other );
    virtual ~G4GeomTestVolPoint();
      // Constructors and virtual destructor

    G4int GetDaughterIndex() const;
      // Accessors

  protected:

    G4int daughterIndex;
};

#endif
