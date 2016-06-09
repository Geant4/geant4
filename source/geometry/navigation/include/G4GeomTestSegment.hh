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
// $Id: G4GeomTestSegment.hh,v 1.2 2003/11/03 17:15:20 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
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
};

#endif
