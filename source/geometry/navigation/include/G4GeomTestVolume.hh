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
// $Id: G4GeomTestVolume.hh 73253 2013-08-22 13:24:02Z gcosmo $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4GeomTestVolume
//
// Class description:
//
// Checks for inconsistencies in the geometric boundaries of a physical
// volume and the boundaries of all its immediate daughters.

// Author: G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4GeomTestVolume_hh
#define G4GeomTestVolume_hh

#include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4GeomTestLogger;

class G4GeomTestVolume
{
  public:  // with description

    G4GeomTestVolume( G4VPhysicalVolume *theTarget,
                      G4double theTolerance=0.0,    // mm
                      G4int numberOfPoints=10000,
                      G4bool theVerbosity=true);
    ~G4GeomTestVolume();
      // Constructor and destructor

    G4double GetTolerance() const;
    void SetTolerance(G4double tolerance);
      // Get/Set error tolerance (default set to 0*mm)
    G4int GetResolution() const;
    void SetResolution(G4int points);
      // Get/Set number of points to check (default set to 10000)
    G4bool GetVerbosity() const;
    void SetVerbosity(G4bool verbosity);
      // Get/Set verbosity mode (default set to true)
    G4int GetErrorsThreshold() const;
    void SetErrorsThreshold(G4int max);
      // Get/Set maximum number of errors to report (default set to 1)

    void TestRecursiveOverlap( G4int sLevel=0, G4int depth=-1 );
      // Activate overlaps check, propagating recursively to the daughters,
      // with possibility of specifying the initial level in the volume tree
      // and the depth (default is the whole tree).
      // Be careful: depending on the complexity of the geometry, this
      // could require long computational time

  private:

    G4VPhysicalVolume *target;        // Target volume
    G4double tolerance;               // Error tolerance
    G4int resolution;                 // Number of points to test
    G4int maxErr;                     // Maximum number of errors to report
    G4bool verbosity;                 // Verbosity level for overlaps check
};

#endif
