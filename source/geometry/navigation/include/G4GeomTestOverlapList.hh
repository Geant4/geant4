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
// $Id: G4GeomTestOverlapList.hh,v 1.2 2003/11/03 17:15:20 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestOverlapList
//
// Class description:
//
// A list of line segments that are found inside two
// separate daughter volumes, indicating a geometry
// overlap error.
//
// This class relies on the compiler generated copy constructor
// and assignment operator.

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestOverlapList_hh
#define G4GeomTestOverlapList_hh

#include "G4GeomTestErrorList.hh"

class G4VPhysicalVolume;

class G4GeomTestOverlapList : public G4GeomTestErrorList
{
  public:  // with description
  
    G4GeomTestOverlapList();
    G4GeomTestOverlapList( const G4VPhysicalVolume *theMother,
                           G4int theDaughter1, G4int theDaughter2 );
    virtual ~G4GeomTestOverlapList();
       // Constructors and virtual destructor
  
    G4bool operator==( const G4GeomTestOverlapList &other ) const;
    G4bool operator< ( const G4GeomTestOverlapList &other ) const;
      // Comparison operators, based on daughter index numbers

    const G4VPhysicalVolume *GetDaughter1() const;
    const G4VPhysicalVolume *GetDaughter2() const;
    G4int GetDaughter1Index() const;
    G4int GetDaughter2Index() const;
      // Return pointers to volumes
  
    void GetDaught1Points( G4int, G4ThreeVector &, G4ThreeVector & ) const;
    void GetDaught2Points( G4int, G4ThreeVector &, G4ThreeVector & ) const;
      // Return start and end points in various coordinate systems

  private:

    G4int daughter1, daughter2;
};

#endif
