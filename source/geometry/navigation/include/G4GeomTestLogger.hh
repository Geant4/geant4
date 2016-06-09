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
// G4GeomTestLogger
//
// Class description:
//
// Abstract base class that defines the interface of a class that
// accepts error reports from the geometry test
//

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestLogger_hh
#define G4GeomTestLogger_hh

#include "G4Types.hh"
#include "G4ThreeVector.hh"

class G4VSolid;
class G4GeomTestOverlapList;
class G4GeomTestOvershootList;

class G4GeomTestLogger
{
  public:  // with description

    G4GeomTestLogger() {}
    virtual ~G4GeomTestLogger() {}
      // Constructor and virtual destructor

    virtual void SolidProblem( const G4VSolid *solid, 
                               const G4String &message,
                               const G4ThreeVector &point ) = 0;
      // Accept an error from a local calculation of
      // a solid, probably an inconsistency in geometry

    virtual void NoProblem( const G4String &message ) = 0;
      // Report a message of no error.

    virtual void OverlappingDaughters( const G4GeomTestOverlapList *) = 0;
      // Accept an error due to the overlap of two volumes

    virtual void OvershootingDaughter( const G4GeomTestOvershootList *) = 0;
      // Accept an error due to a daughter sticking outside
      // a mother volume
};

#endif
