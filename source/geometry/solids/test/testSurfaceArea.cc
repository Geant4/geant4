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
// $Id: testSurfaceArea.cc,v 1.2 2007-02-12 11:29:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//  by Hans Dierckx to test surface Area

// Validation of the Monte Carlo method to calculate arbitrary surface areas. 
// For the different CSG solids available, the surface is calculated
// analytically and by the Monte Carlo method (after conversion into a
// boolean intersection solid). Both results are compared afterwards.

#include <assert.h>
#include <cmath>
#include <ctime>

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

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
#include "G4ExtrudedSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

#include "Randomize.hh"

G4bool testG4Surf()
{
    G4double surf, MCsurf;
    G4int tic, toc;
    
    G4Box w("Bigger Box",5,5,5);
    G4cout << "\n\nCalculating Box:" << G4endl;
    G4Box box("Test Box",1,0.5,3); 
    // fDx fDy fDz
    surf = box.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid box2("Box 2", &box, &w);
    tic=clock();
    MCsurf= box2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    
    
   

    G4cout << " Calculating Tube:" << G4endl;    
    //    G4Tubs cyl("Test Cyl",1,4,3.5,pi/5,3*pi/2);
    G4double R= std::sqrt(100./16./pi);
    G4Tubs cyl("Test Cyl",0.5*R,1.5*R,1.5*R,0,twopi);
    //pRmin pRmax pDz pSPhi pDPhi
    surf = cyl.GetSurfaceArea();
    G4cout << cyl.GetCubicVolume() << G4endl;
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid cyl2("Box 2", &cyl, &w);
    tic=clock();
    MCsurf=cyl2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
     
    G4cout << " Calculating Cone:" << G4endl;
    G4Cons cone("Test Cone",0.5,1,2,2.5,4,0,4*pi/3); 
    //pRmin1 pRmax1 pRmin2, pRmax2, pDz, pSphi pDphi
    surf = cone.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid cone2("Box 2", &cone, &w);
    tic=clock();
    MCsurf=cone2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    
    
    G4cout << " Calculating General G4Trap:" << G4endl;    
    G4Trap trap("Test Cyl",3, 20/180*pi,5/180*pi,  2,  1.5, 2,    pi/18,    0.8 ,  0.5,  0.7,   pi/18);
                          //Dz,Theta,   Phi,      Dy1, Dx1, Dx2, alpha1, Dy2, Dx3,Dx4, alpha2 
    surf = trap.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid trap2("Trap 2", &trap, &w);
    tic=clock();
    MCsurf=trap2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    

            
    G4cout << " Calculating Parallelepiped:" << G4endl;
    G4Para para("Test Para", 1.5, 2, 3, pi/18,    pi/9,  5/180*pi);
                          // Dx,Dy,Dz, alph, theta, phi
    surf = para.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid para2("Para 2", &para, &w);
    tic=clock();
    MCsurf=para2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    

    
    G4cout << " Calculating Trapezoid:" << G4endl;
    // G4Trd trd("Test Trd", 3,1, 4,1.5, 3);
    G4double a = std::sqrt(2.5);
    G4Trd trd("Test Trd", a,2*a, 2*a,a,std::sqrt(3.)*0.5*a);
    //Dx1,Dx2,Dy1,Dy2,Dz
    surf = trd.GetSurfaceArea();
    G4cout << trd.GetCubicVolume() << G4endl;
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid trd2("Trd 2", &trd, &w);
    tic=clock();
    MCsurf=trd2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    

    
    G4cout << " Calculating Sphere:" << G4endl;
    G4Sphere sp("Test Sphere",2,3, 3*pi/4,pi, pi/6, pi);
    //Rmin, Rmax, Sphi, Dphi, STheta, Dtheta
    surf = sp.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid sp2("Trd 2", &sp, &w);
    tic=clock();
    MCsurf=sp2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    
   
    
    G4cout << " Calculating Torus:" << G4endl;
    G4Torus tor("Test Torus", 0.2,2,2.5, 0,0.2);
    //Rmin Rmax Rtor Sphi Dphi
    surf = tor.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid tor2("Tor 2", &tor, &w);
    tic=clock();
    MCsurf=tor2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    


    G4cout << " Calculating Orb:" << G4endl;
    G4Orb orb("Test Orb",2);
    surf = orb.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    G4IntersectionSolid orb2("Box 2", &orb, &w);
    tic=clock();
    MCsurf=orb2.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;
    
  

    G4cout << " Calculating Extruded-Solid (against box):" << G4endl;
    // Box above defined as Xtru ...
    std::vector<G4TwoVector> polygon;
    polygon.push_back(G4TwoVector(-1.,-0.5));
    polygon.push_back(G4TwoVector(-1., 0.5));
    polygon.push_back(G4TwoVector( 1., 0.5));
    polygon.push_back(G4TwoVector( 1.,-0.5));

    G4ExtrudedSolid xtru("xtru-box", polygon, 3,
                         G4TwoVector(), 1.0, G4TwoVector(), 1.0);

    surf = xtru.GetSurfaceArea();
    G4cout <<"Surface = " << surf << G4endl;
    tic=clock();
    MCsurf=box.GetSurfaceArea();
    toc = clock()-tic;
    G4cout << "Elapsed time= "<< static_cast<G4double>(toc)/CLOCKS_PER_SEC << G4endl;
    G4cout <<"MC result = " << MCsurf << G4endl;
    G4cout << "deviation = " << (surf-MCsurf)/surf*100 << " %"<< G4endl;
    G4cout <<"*******************" << G4endl;

    /////////////////////////////////////////////////////
    return true;
}

int main()
{

// initialise random generator:
CLHEP::RanluxEngine defaultEngine(1234567,4);
CLHEP::HepRandom::setTheEngine(&defaultEngine);
G4int seed = time(NULL);
CLHEP::HepRandom::setTheSeed(seed);
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Surf());
    return 0;
}
