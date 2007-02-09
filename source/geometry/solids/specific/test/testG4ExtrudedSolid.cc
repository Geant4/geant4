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
// $Id: testG4ExtrudedSolid.cc,v 1.1 2007-02-09 12:05:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// testG4ExtrudedSolid
//
//  Test file for class G4ExtrudedSolid
//
//             Ensure asserts are compiled in
//
// --------------------------------------------------------------------

#include <assert.h>
#include <cmath>

#include "globals.hh"

#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4ExtrudedSolid.hh"

G4bool testExtrudedSolid0()
{
  // Xtru with concave polygon

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm, -30.*cm));
  
  G4ExtrudedSolid* XtruS0 
    = new G4ExtrudedSolid("XtruS0", polygon, 25.*cm, 
                 G4TwoVector(-20.*cm, 10.*cm), 1.5, G4TwoVector(), 0.5);

  G4ThreeVector testPoint(-37.5, -58.3, 242);
  G4cout << testPoint
         << "  Inside via G4ExtrudedSolid: " 
         << XtruS0->Inside(testPoint)
         << "  Inside via G4TessellatedSolid: "
         << XtruS0->G4TessellatedSolid::Inside(testPoint) << G4endl;

  // no offset, only scale
  //G4ExtrudedSolid* XtruS0 
  //  = new G4ExtrudedSolid("XtruS0", polygon, 25.*cm, 
  //               G4TwoVector(), 1.0, G4TwoVector(), 6.0);

  // no offset, no scale
  //G4ExtrudedSolid* XtruS0 
  //  = new G4ExtrudedSolid("XtruS0", polygon, 25.*cm, 
  //               G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  G4cout << *XtruS0 << G4endl;                     
  XtruS0->G4TessellatedSolid::StreamInfo(G4cout);

  return true;
}

// ---------------------------------------------------------------------------

G4bool testExtrudedSolid1()
{
  // Xtru with triangular polygon

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(  0.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  
  G4ExtrudedSolid* xtruS1 
    = new G4ExtrudedSolid("XtruS1", polygon, 30.*cm, 
                          G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  G4cout << *xtruS1 << G4endl;                     
  xtruS1->G4TessellatedSolid::StreamInfo(G4cout);;                     

  return true;
}

// ---------------------------------------------------------------------------

G4bool testExtrudedSolid2()
{
  // Box defined as Xtru

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  
  G4ExtrudedSolid* XtruS2 
    = new G4ExtrudedSolid("XtruS2", polygon, 30.*cm, 
                 G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  G4cout << *XtruS2 << G4endl;                     
  XtruS2->G4TessellatedSolid::StreamInfo(G4cout);;                     

  return true;
}

// ---------------------------------------------------------------------------

int main()
{
  assert(testExtrudedSolid0());
  assert(testExtrudedSolid1());
  assert(testExtrudedSolid2());
  return 0;
}
