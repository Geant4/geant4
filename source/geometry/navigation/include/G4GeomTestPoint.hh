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
// $Id: G4GeomTestPoint.hh,v 1.3 2006-06-29 18:35:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

    G4GeomTestPoint& operator=(const G4GeomTestPoint& other);
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
