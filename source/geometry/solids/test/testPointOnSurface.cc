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
// $Id: testPointOnSurface.cc,v 1.5 2007-02-12 11:29:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// 
// Test for GetPointOnSurface() method of various solids.
// Returns 0 if generated random points are located on surface.
//
// Author: Dionysios Anninos
//
// --------------------------------------------------------------
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "G4TwoVector.hh"
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
#include "G4Polyhedra.hh"
#include "G4EllipticalCone.hh"
#include "G4ExtrudedSolid.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool checkBox(G4int N)
{
  
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4BOX ******************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;

  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
  
  G4Timer time;
  
  G4Box  t1("Solid Box #1", 
	    20*cm,
	    20*cm,
	    20*cm);
  
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check Box had "<<n<<" inconsistencies..."<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkEllipsoid(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4ELLIPSOID ************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
  
  G4Timer time;
  
  G4Ellipsoid t1("Solid Ellipsoid #1", 
		 20*cm,
		 15*cm,
		 35*cm,
		 -10*cm,
		  10*cm);
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check Ellipsoid had "<<n<<" inconsistencies..."<< G4endl;
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkHype(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4HYPE *****************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Timer time;
 
  G4Hype t1("Solid Hype #1", 
	    10*cm,
	    20*cm,
	    75*deg,
	    75*deg,
	    10*cm);
  
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check Hype had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkEllipticalTube(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4ELLIPTICALTUBE *******"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Timer time;
  
  G4EllipticalTube t1("Solid Ell Tube #1", 
		      10*cm,
		      20*cm,
		      15*cm);
  
  time.Start();

  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check EllTube had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkOrb(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4ORB ******************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
 
  G4Timer time;

  G4Orb t1("Solid Orb #1",  
	   10*cm);

  time.Start();

  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();

  G4cout <<" Check Orb had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkTorus(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4TORUS ****************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;       
  G4int i=0,n=0;
  EInside surf;
  
  G4Timer time;
 
  G4Torus t1("Torus  #1",   
	     5*cm,
	     8*cm,
	     20*cm,
	     33*deg,
	     270*deg);
  
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
   
  time.Stop();
  
  G4cout <<" Check Torus had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkTubs(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4TUBS *****************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;

  G4Timer time;
 
  G4Tubs t1("Tubs #1", 
	     5*cm,
	     8*cm,
	     11*cm,
	     42*deg,
	     120*deg);
  
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check Tubs had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkCons(G4int N)
{ 
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4CONS *****************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;   
  G4int i=0,n=0; 
  EInside surf;
  
  G4Timer time;
   
  G4Cons t1("Cons #1", 
	    6*cm, 10*cm,   
	    8*cm, 9*cm, 
	    10*cm,    
	    43*deg,  221*deg);
  
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }

  time.Stop();
  
  G4cout <<" Check Cons had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkTrap(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4TRAP *****************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
  
  G4Timer time;
   
  G4Trap t1("Trap #1", 
	    20*mm,  10*deg, 5*deg, 
	    5*mm,   5*mm,   10*mm,    
	    20*deg, 4*mm,   7*mm, 11*mm, 20*deg);
  
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check Trap had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkPara(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4PARA *****************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Timer time;

  G4Para t1("Para #1", 
	    5*mm,   10*mm,  20*mm, 
	    15*deg, 20*deg, 70*deg);

  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();

  G4cout <<" Check Para had "<<n<<" inconsistencies"<< G4endl;
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkTrd(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4TRD ******************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
    
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Timer time;

  G4Trd t1("Trd #1", 
	    10*mm,  6*mm,  12*mm, 
	    8*mm,  15*mm);
  
  time.Start();

  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();

  G4cout <<" Check Trd had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkSphere(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4SPHERE ***************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;

  G4Timer time;
    
  G4Sphere t1("Sphere", 
	      8*mm,  12*mm,     
	      43*deg,350*deg,  
	      21*deg,50*deg);  
  
  time.Start();
   
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }

  time.Stop();
  
  G4cout <<" Check Sphere had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
} 

G4bool checkPolycone(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4POLYCONE *************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
       
  G4double zPlanes[5] = {0., 1., 3.,  5., 10.};
  G4double rInner[5]  = {6., 7., 2.,  2., 10.};
  G4double rOuter[5]  = {8., 8., 10.,10., 15.};     
  
  G4Polycone t1("aPcone",  
		269*deg, 342*deg,       
		5,zPlanes,rInner,rOuter); 
  
  G4Timer time;
  time.Start();
   
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }
  
  time.Stop();
  
  G4cout <<" Check PolyCone had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}  

G4bool checkTet(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4TET ******************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
   
  G4Tet t1("aTet", 
	   G4ThreeVector(2.*cm,0.,10.*cm),
	   G4ThreeVector(0.,3.*cm,0.*cm),
	   G4ThreeVector(10.,10.*cm,5.*cm),
	   G4ThreeVector(0.*cm,5.*cm,5.*cm)); 
  
  G4Timer time;
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }

  time.Stop();
  
  G4cout <<" Check Tet had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}

G4bool checkPolyhedra(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4POLYHEDRA ************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
 
  G4ThreeVector point;
  G4bool what = true;   
  G4int i=0,n=0;
  EInside surf;
       
  G4double zPlanes[5] = {-1., 10., 15., 25., 30.};
  G4double rInner[5]  = {0., 5., 0.,  7., 1.};
  G4double rOuter[5]  = {21., 6., 15., 15., 38.}; 
    
//   G4double z[10] = {30.,25.,15.,10.,-1.,-1.,10.,15.,25.,30.};
//   G4double r[10] = {1.,7.,0.,5.,0.,21.,6.,15.,15.,38.};
    
  G4Polyhedra t1("aPhedra",  
		53.*deg, 163.*deg,         
		8,5,zPlanes,rInner,rOuter); 
     
//   G4Polyhedra t1("aPhedra",  
// 		 53.*deg, 163.*deg,         
// 		 8, 10, r, z);
  
  G4Timer time;
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();  
    surf  = t1.Inside(point);
    if(surf != kSurface)
    { 
      n++; what = false;
      G4cout <<" x "<<point.x()<<" y "<<point.y()<<" z "<<point.z()<<G4endl;
    }  
  }
  
  time.Stop();

  G4cout <<" Check Polyhedra had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}   

G4bool checkEllipticalCone(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4ELLIPTICALCONE *******"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;   
  G4int i=0,n=0;
  EInside surf;
       
  G4EllipticalCone t1("aElliCone",  
		      20*cm,25*cm,20*cm, 10*cm); 
   
  G4Timer time;
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();  
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }

  time.Stop();
  
  G4cout <<" Check EllipticalCone had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}   

G4bool checkExtrudedSolid(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4EXTRUDEDSOLID ********"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;   
  G4int i=0,n=0;
  EInside surf;
       
  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm, -30.*cm));

  G4ExtrudedSolid t1("Concave_XTRU", polygon, 25.*cm,
                     G4TwoVector(-20.*cm, 10.*cm), 1.5, G4TwoVector(), 0.5);

  G4Timer time;
  time.Start();
  
  for(i=0; i<N; i++)
  {
    point = t1.GetPointOnSurface();  
    surf  = t1.Inside(point);
    if(surf != kSurface){ n++; what = false; }  
  }

  time.Stop();
  
  G4cout <<" Check ExtrudedSolid had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}   

int main() 
{ 
  G4bool what;    
  G4int N = 1000000;
  
  G4cout <<G4endl;
  G4cout <<"********************************************************************************"<<G4endl;
  G4cout <<"**************** TEST GET POINT ON SURFACE METHOD ******************************"<<G4endl;
  G4cout <<"******************** FOR "<<N<<" RANDOM POINTS *********************************"<<G4endl;
  G4cout <<"********************************************************************************"<<G4endl;
   
  G4cout <<G4endl<<G4endl;

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
  what = checkPolyhedra(1000000);      
  what = checkEllipticalCone(1000000);
  what = checkExtrudedSolid(1000000);
  
  G4cout <<G4endl;
    
  G4cout <<"********************************************************************************"<<G4endl;
  G4cout <<"********************** END OF POINT-ON-SURFACE TEST ****************************"<<G4endl;
  G4cout <<"********************************************************************************"<<G4endl;  
  G4cout <<G4endl;

  return 0;     
}
       
   
 
