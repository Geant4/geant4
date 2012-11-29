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
// G4GeomTestSegment
//
// Class description:
//
// Locate all points on a line that intersect a solid,
// and return them in no particular order

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestSegment_hh
#define G4GeomTestSegment_hh

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4GeomTestPoint.hh"

#include <vector>

class G4VSolid;
class G4GeomTestLogger;

class G4GeomTestSegment
{
  public:  // with description

    G4GeomTestSegment( const G4VSolid *theSolid,
                       const G4ThreeVector &theP,
                       const G4ThreeVector &theV,
                             G4GeomTestLogger *logger );
      // Constructor

    const G4VSolid *GetSolid() const;
    const G4ThreeVector &GetP() const;
    const G4ThreeVector &GetV() const;
    const G4GeomTestPoint &GetPoint( G4int i ) const;
    G4int GetNumberPoints() const;
      // Accessors

  private:

    void FindPoints( G4GeomTestLogger *logger );
    void FindSomePoints( G4GeomTestLogger *logger, G4bool forward );
    void PatchInconsistencies( G4GeomTestLogger *logger );
  
    const G4VSolid * solid;
    const G4ThreeVector p0,v;
  
    std::vector<G4GeomTestPoint> points;

    G4double kCarTolerance;
};

#endif
