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
// $Id: G4tgrPlaceSimple.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrPlaceSimple
//
// Class description:
//
// Class to describe a simple positioning of a G4tgrVolume inside
// another G4tgrVolume.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrPlaceSimple_h
#define G4tgrPlaceSimple_h

#include "globals.hh"

#include <vector>

#include "G4ThreeVector.hh"
#include "G4tgrPlace.hh"

class G4tgrPlaceSimple : public G4tgrPlace
{
  public:  // with description

    G4tgrPlaceSimple();
   ~G4tgrPlaceSimple();

    G4tgrPlaceSimple( const std::vector<G4String>& wl );

    // Accessors

    const G4String& GetRotMatName() const { return theRotMatName; }
    G4ThreeVector GetPlacement() const { return thePlace; }

    friend std::ostream& operator<<(std::ostream& os,
                                    const G4tgrPlaceSimple& obj);
  protected:

    G4ThreeVector thePlace;
      // The position with respect to parent

    G4String theRotMatName;
      // The rotation matrix (by name, as the rotations
      // matrices are not yet created)
};

#endif
