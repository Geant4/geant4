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
// $Id: testG4NavigatorSafety.cc,v 1.1 2010-01-12 16:57:27 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//   Create a tubular "calorimeter". Generate random points along x, y & z
//   axes, printing location, steps & safeties. Compare results of standard
//   voxel safety calculations with "exact safety" computed values.
//
//   Arguments: Number of points to generate [Default: 1000]
//              Initial random seed [Default: CLHEP default], between [0,255]

#include <assert.h>
#include "G4ios.hh"
#include <stdlib.h>
#include <vector>

#include "globals.hh"

#include "G4Timer.hh"
#include "ApproxEqual.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

// Parameters for building a tubular calorimeter:
// an array of interlocking complete tubes inside a box
//
static const G4double kWorldhsize = 100*mm;
static const G4double kTubeHalfHeight = 10*mm;
static const G4double kTubeRadius = 5*mm;
static const G4double kTubeNoRow = 10;
static const G4double kTubeNoColumn = 11; // Odd for symmetrical array

static const G4double kBoxDx=kTubeNoRow*kTubeRadius;
static const G4double yDelta=2.0*kTubeRadius*std::sin(pi/3.0);
static const G4double kBoxDy=(kTubeNoColumn-1)*yDelta*0.5+kTubeRadius;
static const G4double kBoxDz=kTubeHalfHeight;

G4bool compare = false;
std::vector<G4ThreeVector> kPoints;
std::vector<std::pair<G4double, G4double> > kSafeties;

G4VPhysicalVolume* BuildGeometry()
{
  G4double bigXStart=-(kTubeNoRow-1)*kTubeRadius;
  G4double smallXStart=bigXStart+kTubeRadius;

  G4double bigYStart=-(kTubeNoColumn-1)*yDelta*0.5;
  G4double smallYStart=bigYStart+yDelta;

  G4int row,column;

  G4Box *worldBox = new G4Box ("World Box",kWorldhsize,kWorldhsize,kWorldhsize);
    // World box

  G4LogicalVolume *myWorldLog=
    new G4LogicalVolume(worldBox,0,"World",0,0,0);
    // Logical with no material,field, sensitive detector or user limits

  G4PVPlacement *myWorldPhys=
    new G4PVPlacement(0,G4ThreeVector(0,0,0),"World",myWorldLog,0,false,0);
    // World physical volume

  G4Box *calBox = new G4Box ("Cal Box",kBoxDx,kBoxDy,kBoxDz);
  G4Tubs *calTube = new G4Tubs("Cal Tube",0,kTubeRadius,kTubeHalfHeight,0,360);

  G4LogicalVolume *myDetectorLog=
    new G4LogicalVolume(calBox,0,"DetectorLog",0,0,0);
    // Logical with no material,field, sensitive detector or user limits

  G4PVPlacement *myDetectorPhys=
    new G4PVPlacement(0,G4ThreeVector(0,0,0),"DetectorPhys",
                      myDetectorLog,myWorldPhys,false,0);
    // Detector physical volume placed in the world

  G4LogicalVolume *calTubLog=
    new G4LogicalVolume(calTube,0,"Cal Crystal",0,0,0);

  G4String tname("Target");
  G4int copyNo=0;
  for (column=0;column<kTubeNoColumn;column+=2)
  {
    for (row=0;row<kTubeNoRow;row++)
    {    
//    G4PVPlacement *calPhys=
      new G4PVPlacement(0,G4ThreeVector(bigXStart+row*kTubeRadius*2.0,
                                        bigYStart+column*yDelta,0),
                        tname,calTubLog,myDetectorPhys,false,copyNo++);
    }
  }
  for (column=0;column<kTubeNoColumn-1;column+=2)
  {
    for (row=0;row<kTubeNoRow-1;row++)
    {
//    G4PVPlacement *calPhys=
      new G4PVPlacement(0,G4ThreeVector(smallXStart+row*kTubeRadius*2.0,
                                        smallYStart+column*yDelta),
                        tname,calTubLog,myDetectorPhys,false,copyNo++);
    }
  }
  return myDetectorPhys;
}

void generatePoints(G4int n)
{
  for (int i=0; i<n; i++)
  {
    G4ThreeVector p(CLHEP::RandFlat::shoot(-kWorldhsize,kWorldhsize),
                    CLHEP::RandFlat::shoot(-kWorldhsize,kWorldhsize),
                    CLHEP::RandFlat::shoot(-kWorldhsize,kWorldhsize));
    kPoints.push_back(p);
  }
}

void computeApproxSafeties(G4VPhysicalVolume *pTopNode)
{
  MyNavigator myNav;
  myNav.SetWorldVolume(pTopNode);
  std::vector<G4ThreeVector>::const_iterator pos;
  for (pos=kPoints.begin(); pos!=kPoints.end(); pos++)
  {
    G4double safety = myNav.ComputeSafety(*pos);
    std::pair<G4double,G4double> s(safety,0.);
    kSafeties.push_back(s);
  }
}

void computeExactSafeties(G4VPhysicalVolume *pTopNode)
{
/*
  G4int i=0;
  MyNavigator myNav;
  myNav.SetWorldVolume(pTopNode);
  std::vector<G4ThreeVector>::const_iterator pos;
  for (pos=kPoints.begin(); pos!=kPoints.end(); pos++)
  {
    G4double safety = myNav.ComputeExactSafety(*pos);
    kSafeties[i].second = safety;
    i++;
  }
  compare = true;
*/
}

void compareSafeties()
{
  G4int n = kPoints.size();
  if (!compare)
  {
    G4cout << "Printing out non-zero safety values computed ..." << G4endl;
    for (G4int i=0; i<n; i++)
    {
      G4double safety = kSafeties[i].first;
      if (safety)
      {
        G4cout << i << ". - Point: " << kPoints[i]
               << " - Safety: " << safety << G4endl;
      }
    }
  }
  else
  {
    G4int diffs=0;
    G4cout << "Printing out cases of different safety values ..." << G4endl;
    for (G4int i=0; i<n; i++)
    {
      if (kSafeties[i].first != kSafeties[i].second)
      {
        G4cout << i << ". - Point: " << kPoints[i]
               << " - Approx safety: " << kSafeties[i].first
               << " - Exact safety: " << kSafeties[i].second << G4endl;
        diffs++;
      }
      G4cout << "Total differences: " << diffs << G4endl;
    }
  }
}

int main(int argc, char* argv[])
{
  G4Timer timer;
  G4int iter=10000;

  if (argc==2)
  {
    G4int num = atoi(argv[1]);
    if (num>=0) { iter = num; }
    else { G4cout << ">>> Invalid number of iterations in input!" << G4endl
                  << "    Sticking to: " << iter << " ..." << G4endl; }
  }
  else if (argc==3)
  {
    G4int num = atoi(argv[1]);
    if (num>=0) { iter = num; }
    else { G4cout << ">>> Invalid number of iterations in input!" << G4endl
                  << "    Sticking to: " << iter << " ..." << G4endl; }
    G4long seed = atoi(argv[2]);
    if (seed>0) { CLHEP::HepRandom::setTheSeed(seed); }
    else { G4cout << ">>> Invalid negative random seed in input!" << G4endl
                  << "    Sticking to default: "
                  << CLHEP::HepRandom::getTheSeed() << " ..." << G4endl; }
  }

  G4VPhysicalVolume *myTopNode;
  myTopNode=BuildGeometry();  // Build the geometry
  G4GeometryManager::GetInstance()->CloseGeometry();  // Voxelise the geometry

  timer.Start();

  generatePoints(iter);
  computeApproxSafeties(myTopNode);
//  computeExactSafeties(myTopNode);
  compareSafeties();

  timer.Stop();

  G4cout << timer << G4endl;

  G4GeometryManager::GetInstance()->OpenGeometry();

  return 0;
}
