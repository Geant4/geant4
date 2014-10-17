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
// Adapted by Hans Dierckx to test surface Area

// The MC method is applied to different solids that have a total surface
// area of 100 units. The method is repeated 20 times, after which the
// average and standard deviation is calculated. This generates the graphs
// in the report, if one would set the right values for N and gamma in
// G4VSolid.cc.

//             Ensure asserts are compiled in

#include <assert.h>
#include <cmath>
#include <ctime>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Cons.hh"
#include "G4BooleanSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"


#include "Randomize.hh"

/*
G4bool testG4Surf()
{
    G4double surf, MCsurf;

    G4cout << "\n\nCalculating Box:" << G4endl;
    G4Box box("Test Box",1,0.5,3); 
    // fDx fDy fDz
    surf = box.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid box2("Box 2", &box, &w);
    MCsurf= box2.GetSurfaceArea();
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation =" << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;

    return true;
}
*/

int main()
{

  // initialise random generator:
  CLHEP::RanluxEngine defaultEngine(1234567,4);
  CLHEP::HepRandom::setTheEngine(&defaultEngine);
  G4int seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);


  // testobj = Sphere
  G4Box w("Bigger Box",5,5,5);
  // G4double R = std::sqrt(25./pi);
  //G4Orb sol("Orb",R);
  
  // G4double R= std::sqrt(100./16./pi);
  // G4Tubs sol("Test Cyl",0.5*R,1.5*R,1.5*R,0,twopi);

  G4double R = std::sqrt(2.5);
  G4Trd sol("Test Trd", R,2*R, 2*R,R,std::sqrt(3.)*0.5*R);
  //Dx1,Dx2,Dy1,Dy2,Dz
  G4double exsurf = sol.GetSurfaceArea();
  G4cout << "Exact result:"<< exsurf << G4endl;
  G4IntersectionSolid sol2("Box 2", &sol, &w);
 
  for(G4int N = 100000; N<=10000001; N *=10){
    G4double scale;
    for(scale=2; scale>=0.00096; scale /=2){
      G4double sum = 0;
      G4double std = 0;
      G4double surf;
      G4double rep = 20;
      for(G4int p=0; p<rep;p++)
      {
         sol2.SetAreaStatistics(N);     // Inside loop to force area
         sol2.SetAreaAccuracy(scale*R); // recomputation for calculating
         surf = sol2.GetSurfaceArea();  // average result ...
         sum =sum+surf;
         std =std+(exsurf-surf)*(exsurf-surf);
      }
    G4double avsurf = sum/rep;
    std = std::sqrt(std/(rep-1));
    G4cout << "MC result: N= "<< N << " Average = "<< avsurf
           << " std= " << std << " scale = "<< scale <<  G4endl;
  }
}
  return 0;
}
