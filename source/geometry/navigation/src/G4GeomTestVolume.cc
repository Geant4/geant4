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
// class G4GeomTestVolume implementation
//
// Author: G.Cosmo, CERN
// --------------------------------------------------------------------

#include <queue>
#include <set>

#include "G4GeomTestVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

//
// Constructor
//
G4GeomTestVolume::G4GeomTestVolume( G4VPhysicalVolume *theTarget,
                                    G4double theTolerance,
                                    G4int numberOfPoints,
                                    G4bool theVerbosity )
  : target(theTarget), tolerance(theTolerance),
    resolution(numberOfPoints), verbosity(theVerbosity)
{;}

//
// Destructor
//
G4GeomTestVolume::~G4GeomTestVolume() {;}

//
// Get error tolerance
//
G4double G4GeomTestVolume::GetTolerance() const
{
  return tolerance;
}

//
// Set error tolerance
//
void G4GeomTestVolume::SetTolerance(G4double tol)
{
  tolerance = tol;
}

//
// Get number of points to check (resolution)
//
G4int G4GeomTestVolume::GetResolution() const
{
  return resolution;
}

//
// Set number of points to check (resolution)
//
void G4GeomTestVolume::SetResolution(G4int np)
{
  resolution = np;
}

//
// Get verbosity
//
G4bool G4GeomTestVolume::GetVerbosity() const
{
  return verbosity;
}

//
// Set verbosity
//
void G4GeomTestVolume::SetVerbosity(G4bool verb)
{
  verbosity = verb;
}

//
// Get errors reporting threshold
//
G4int G4GeomTestVolume::GetErrorsThreshold() const
{
  return maxErr;
}

//
// Set maximum number of errors to report
//
void G4GeomTestVolume::SetErrorsThreshold(G4int max)
{
  maxErr = max;
}

//
// Test overlap in tree
//
void G4GeomTestVolume::TestOverlapInTree() const
{
  std::queue<G4VPhysicalVolume*> volumes;
  std::set<G4LogicalVolume*> checked;

  volumes.push(target);
  while (!volumes.empty())
  {
    G4VPhysicalVolume* current = volumes.front();
    volumes.pop();

    // check overlaps for daughters
    G4LogicalVolume* logical = current->GetLogicalVolume();
    std::size_t ndaughters = logical->GetNoDaughters();
    for (std::size_t i=0; i<ndaughters; ++i)
    {
      G4VPhysicalVolume* daughter = logical->GetDaughter(i);
      daughter->CheckOverlaps(resolution, tolerance, verbosity, maxErr);
    }

    // append the queue of volumes
    G4LogicalVolume* previousLogical = nullptr;
    for (std::size_t i=0; i<ndaughters; ++i)
    {
      G4VPhysicalVolume* daughter = logical->GetDaughter(i);
      G4LogicalVolume* daughterLogical = daughter->GetLogicalVolume();
      if (daughterLogical->GetNoDaughters() == 0) continue;
      G4bool found = (daughterLogical == previousLogical);
      if (!found) found = (checked.find(daughterLogical) != checked.cend());
      if (!found)
      {
        checked.emplace(daughterLogical);
        previousLogical = daughterLogical;
        volumes.push(daughter);
      }
      else
      {
        if (verbosity)
          G4cout << "Checking overlaps in tree of volume " << daughter->GetName()
                 << " (" << daughterLogical->GetSolid()->GetEntityType() << ")"
                 << " is omitted, to avoid duplication" << G4endl;
      }
    }
  }
}

//
// TestRecursiveOverlap
//
void G4GeomTestVolume::TestRecursiveOverlap( G4int slevel, G4int depth )
{
  // If reached requested level of depth (i.e. set to 0), exit.
  // If not depth specified (i.e. set to -1), visit the whole tree.
  // If requested initial level of depth is not zero, visit from beginning
  //
  if (depth == 0) return;
  if (depth != -1) depth--;
  if (slevel != 0) slevel--;

  //
  // As long as we reached the requested
  // initial level of depth, test ourselves
  //
  if ( slevel==0 )
  {
    target->CheckOverlaps(resolution, tolerance, verbosity, maxErr);
  }

  //
  // Loop over unique daughters
  //
  std::set<const G4LogicalVolume *> tested;

  const G4LogicalVolume *logical = target->GetLogicalVolume();
  G4int nDaughter = (G4int)logical->GetNoDaughters();
  for( auto iDaughter=0; iDaughter<nDaughter; ++iDaughter )
  {
    G4VPhysicalVolume *daughter = logical->GetDaughter(iDaughter);
    
    // Tested already?
    //
    // const G4LogicalVolume *daughterLogical =
    //      daughter->GetLogicalVolume();
    // std::pair<std::set<const G4LogicalVolume *>::iterator, G4bool>
    //       there = tested.insert(daughterLogical);
    // if (!there.second) continue;

    //
    // Recurse
    //
    G4GeomTestVolume vTest( daughter, tolerance, resolution, verbosity );
    vTest.SetErrorsThreshold(maxErr);
    vTest.TestRecursiveOverlap( slevel,depth );
  }
}
