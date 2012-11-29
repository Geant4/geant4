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
