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
// $Id: testPointOnSurface.cc,v 1.2 2005-08-04 11:26:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// 
// Test for GetPointOnSurface() method of various solids.
// Returns 0 if generated random points are located on surface.
//
// Author: Dionisyos Anninos
//
// --------------------------------------------------------------
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"

#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4Hype.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh" 
#include "G4Tet.hh"

#include "Randomize.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool checkBox(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Box  t1("Solid Box #1", 
	    20*cm,
	    20*cm,
	    20*cm);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Box had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkEllipsoid(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Ellipsoid t1("Solid Ellipsoid #1", 
		 20*cm,
		 15*cm,
		 35*cm,
		 -10*cm,
		  10*cm);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Ellipsoid had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkHype(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Hype t1("Solid Hype #1", 
	    10*cm,
	    20*cm,
	    75*deg,
	    75*deg,
	    10*cm);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Hype had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkEllipticalTube(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4EllipticalTube t1("Solid Ell Tube #1", 
		      10*cm,
		      20*cm,
		      15*cm);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check EllTube had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkOrb(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Orb t1("Solid Orb #1",  
	   10*cm);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Orb had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkTorus(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;       
  G4int i=0,n=0;
  EInside surf;
 
  G4Torus t1("Torus  #1",   
	     5*cm,
	     8*cm,
	     20*cm,
	     33*deg,
	     270*deg);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
   
  G4cout <<"Check Torus had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkTubs(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Tubs t1("Tubs #1", 
	     5*cm,
	     8*cm,
	     11*cm,
	     42*deg,
	     120*deg);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Tubs had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkCons(G4int N)
{ 
 
  G4ThreeVector point;
  G4bool what = true;   
  G4int i=0,n=0; 
  EInside surf;
   
  G4Cons t1("Cons #1", 
	    6*cm, 10*cm,   
	    8*cm, 9*cm, 
	    10*cm,    
	    43*deg,  221*deg);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Cons had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkTrap(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Trap t1("Trap #1", 
	    20*mm,  10*deg, 5*deg, 
	    5*mm,   5*mm,   10*mm,    
	    20*deg, 4*mm,   7*mm, 11*mm, 20*deg);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Trap had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkPara(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Para t1("Para #1", 
	    5*mm,   10*mm,  20*mm, 
	    15*deg, 20*deg, 70*deg);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Para had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkTrd(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Trd t1("Trd #1", 
	    10*mm,  6*mm,  12*mm, 
	    8*mm,  15*mm);
  
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Trd had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}

G4bool checkSphere(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
    
  G4Sphere t1("Sphere", 
	      8*mm,  12*mm,     
	      43*deg,350*deg,  
	      21*deg,50*deg);  
   
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  G4cout <<"Check Sphere had "<<n<<" inconsistencies"<< G4endl; 
  return what;
} 

G4bool checkPolycone(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4double zPlanes[5] = {0., 1., 3.,5., 10.};
  G4double rInner[5]  = {6., 7.,  2., 2., 10.};
  G4double rOuter[5]  = {8., 8., 10.,10., 15.};     
  
  G4Polycone t1("aPcone",  
		269*deg, 342*deg,       
		5,zPlanes,rInner,rOuter); 
   
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ 
      n++; what = false;
      //G4cout <<" x "<<point.x()<<" y "<<point.y()<<" z "<<point.z()<<G4endl;
     }  
  }
  
  G4cout <<"Check PolyCone had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}  

G4bool checkTet(G4int N)
{

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Tet t1("aTet", 
	   G4ThreeVector(2.*cm,0.,10.*cm),
	   G4ThreeVector(0.,3.*cm,0.*cm),
	   G4ThreeVector(10.,10.*cm,5.*cm),
	   G4ThreeVector(0.*cm,5.*cm,5.*cm)); 
   
  for(i=0; i<N; i++){
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ 
      n++; what = false;
     }  
  }
  
  G4cout <<"Check Tet had "<<n<<" inconsistencies"<< G4endl; 
  return what;
}
  
  
int main() 
{ 
  G4bool what;     
   
   what = checkBox(1000000);        
   what = checkEllipsoid(1000000);    
   what = checkHype(1000000);       
   what = checkEllipticalTube(1000000);          
   what = checkOrb(1000000); 
   what = checkTorus(1000000);
   what = checkTubs(1000000);     
   what = checkCons(1000000);     
   what = checkTrap(1000000);      
   what = checkPara(1000000);  
   what = checkTrd(1000000); 
   what = checkSphere(1000000);  
   what = checkPolycone(1000000);    
   what = checkTet(1000000);
   return 0;
}
      
   
 
