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
// $Id: G4GeomTestLogger.hh,v 1.2 2003/11/03 17:15:20 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
