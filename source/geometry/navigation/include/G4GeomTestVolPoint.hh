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
// $Id$
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
// --------------------------------------------------------------------
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

    G4GeomTestVolPoint& operator=(const G4GeomTestVolPoint& other);
      // Assignment operator

    G4int GetDaughterIndex() const;
      // Accessors

  protected:

    G4int daughterIndex;
};

#endif
