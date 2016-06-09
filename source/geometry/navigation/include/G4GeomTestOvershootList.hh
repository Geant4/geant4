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
// $Id: G4GeomTestOvershootList.hh,v 1.2 2003/11/03 17:15:20 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestOvershootList
//
// Class description:
//
// A list of line segments inside a specific volume daughter
// but outside its mother volume, indicating a geometry error

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestOvershootList_hh
#define G4GeomTestOvershootList_hh

#include "G4GeomTestErrorList.hh"

class G4GeomTestOvershootList : public G4GeomTestErrorList
{
  public:  // with description
  
    G4GeomTestOvershootList();
    G4GeomTestOvershootList( const G4VPhysicalVolume *theMother,
                                   G4int theDaughter );
    virtual ~G4GeomTestOvershootList();
      // Constructors and virtual destructor

    G4bool operator==( const G4GeomTestOvershootList &other ) const;
    G4bool operator< ( const G4GeomTestOvershootList &other ) const;
      // Comparison operators, based on daughter index

    const G4VPhysicalVolume *GetDaughter() const;
    G4int GetDaughterIndex() const;
      // Return pointers to volumes

    void GetDaughtPoints( G4int i, G4ThreeVector &s1, G4ThreeVector &s2 ) const;
      // Return start and end points in various
      // coordinate systems

  private:

    G4int daughter;
};

#endif
