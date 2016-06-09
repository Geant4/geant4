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
// $Id: G4GeomTestVolume.hh,v 1.3 2006-06-29 18:35:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4GeomTestVolume_hh
#define G4GeomTestVolume_hh

#include "G4ThreeVector.hh"
#include "G4VisExtent.hh"
#include "G4GeomTestOverlapList.hh"
#include "G4GeomTestOvershootList.hh"

#include <map>

class G4VPhysicalVolume;
class G4GeomTestLogger;

class G4GeomTestVolume
{
  public:  // with description

    G4GeomTestVolume( const G4VPhysicalVolume *theTarget,
                            G4GeomTestLogger *theLogger,
                            G4double theTolerance=1E-4 );  // mm
    ~G4GeomTestVolume();
      // Constructor and destructor

    G4double GetTolerance() const;
    void SetTolerance(G4double tolerance);
      // Get/Set error tolerance (default set to 1E-4*mm)

    void TestCartGridXYZ( G4int nx=100, G4int ny=100, G4int nz=100 );
    void TestCartGridX( G4int ny=100, G4int nz=100 );
    void TestCartGridY( G4int nz=100, G4int nx=100 );
    void TestCartGridZ( G4int nx=100, G4int ny=100 );
      // Test using a grid of lines parallel to a cartesian axis
  
    void TestCartGrid( const G4ThreeVector &g1,
                       const G4ThreeVector &g2,
                       const G4ThreeVector &v,
                             G4int n1, 
                             G4int n2 );
      // Test using a grid of parallel lines
      //  g1 = First grid axis
      //  g2 = Second grid axis
      //  v  = Direction of lines
      //  n1 = Number of grid points along g1
      //  n2 = Number of grid points along g2
      // The spread of the grid points are automatically calculated
      // based on the extent of the solid

    void TestRecursiveCartGrid( G4int nx=100, G4int ny=100, G4int nz=100,
                                G4int sLevel=0, G4int depth=-1 );
      // Test using a grid, propagating recursively to the daughters, with
      // possibility of specifying the initial level in the volume tree and
      // the depth (default is the whole tree).
      // Be careful: depending on the complexity of the geometry, this
      // could require long computational time

    void TestCylinder( G4int nPhi=90, G4int nZ=50, G4int nRho=50,
                       G4double fracZ=0.8,  G4double fracRho=0.8,
                       G4bool usePhi=false    );
      // Test using a set of lines in a cylindrical
      // pattern of gradually increasing mesh size
      //       nPhi    = Number lines per phi
      //       nZ      = Number of z points
      //       nRho    = Number of rho points
      //       fracZ   = Fraction scale for points along z
      //       fracRho = Fraction scale for points along rho
      //       usePhi  = Include phi set of lines
      // Define a set of rho values such that:
      //       rho0 = Size of volume
      //       rho1 = frac*rho0
      //       rho2 = frac*rho1
      //       ... etc
      // And define a set of z values
      //       z0   = z size of volume
      //       z1   = fracZ*z0
      //       z2   = fracZ*z1
      //       .... etc
      // And define a set of nPhi phi values, evenly
      // distributed in phi
      //
      // Three sets of lines are tested:
      //   * Imagine the set of lines parallel to the z axis
      //     through each rho point, at a phi angle taken the
      //     set of phi angles
      //   * Imagine the set of lines running perpendicular
      //     to the z axis and through a point on the z axis
      //     at +/- each z point and at an angle taken from the
      //     set of phi values
      //   * If usePhi==true, now take each pair of lines from the 
      //     above two sets and imagine the line through the
      //     intersection and perpendicular to both

    void TestRecursiveCylinder( G4int nPhi=90, G4int nZ=50, G4int nRho=50,
                                G4double fracZ=0.8,  G4double fracRho=0.8,
                                G4bool usePhi=false,
                                G4int sLevel=0, G4int depth=-1 );
      // Test using a set of lines in a cylindrical pattern of gradually
      // increasing mesh size, propagating recursively to the daughters, with
      // possibility of specifying the initial level in the volume tree and
      // the depth (default is the whole tree).
      // Be careful: depending on the complexity of the geometry, this
      // could require long computational time

    void TestOneLine( const G4ThreeVector &p, const G4ThreeVector &v );
      // Test using a single line, specified by a point and direction,
      // in the coordinate system of the target volume

    void TestRecursiveLine( const G4ThreeVector &p, const G4ThreeVector &v,
                            G4int sLevel=0, G4int depth=-1 );
      // Test using a single line, specified by a point and direction,
      // propagating recursively to the daughters, with possibility of
      // specifying the initial level in the volume tree and
      // the depth (default is the whole tree).

    void ReportErrors();
      // Tabulate and report all errors so far encountered

    void ClearErrors();
      // Clear list of errors
  
  protected:

    const G4VPhysicalVolume *target;  // Target volume
    G4GeomTestLogger *logger;         // Error logger
    G4double tolerance;               // Error tolerance
    G4VisExtent extent;               // Geometric extent of volume
  
    std::map<G4long,G4GeomTestOverlapList> overlaps;
      // A list of overlap errors, keyed by the
      // daughter1*numDaughter+daughter2, where daughter1
      // is the smaller of the two daughter indices

    std::map<G4long,G4GeomTestOvershootList> overshoots;
      // A list of overshoot errors, keyed by the
      // daughter number
};

#endif
