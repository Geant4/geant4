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
// $Id: testPointOnSurfacePoly.cc,v 1.3 2008-04-25 08:50:00 gcosmo Exp $
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

#include "G4Polycone.hh" 
#include "G4Polyhedra.hh" 

#include "G4Timer.hh"
#include "Randomize.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool checkPolycone(G4int N)
{
  G4cout<<"**************************************"<<G4endl;
  G4cout<<"************* G4POLYCONE *************"<<G4endl;
  G4cout<<"**************************************"<<G4endl<<G4endl;
  
  G4ThreeVector point;
  G4bool what = true;  
  G4int i=0,n=0;
  EInside surf;
  //-------------------------------------------
  // Original Polycone Test
  //-------------------------------------------    
  // G4double zPlanes[5] = {0., 1., 3.,  5., 10.};
  // G4double rInner[5]  = {6., 7., 2.,  2., 10.};
  // G4double rOuter[5]  = {8., 8., 10.,10., 15.};    
  // G4Polycone t1("aPcone",  
  //		0., 2.*pi,       
  //		5,zPlanes,rInner,rOuter); 
  
  //-------------------------------------------
  // Polycone Test Z1=Z2
  //-------------------------------------------    
  // G4double zPlanes[4] = {0. ,27.5, 27.5, 59.};
  // G4double rInner[4]  = {555.,555., 572.,  572.};
  // G4double rOuter[4]  = {669., 669., 652.,652.};    
  // G4Polycone t1("aPcone",  
  //		0., 2.*pi,       
  //		4,zPlanes,rInner,rOuter);  
  //-------------------------------------------
  // More Complex  Polycone Test Z1=Z2
  //-------------------------------------------    
  // G4double zPlanes[6] = {0., 61., 61.,  72.3, 124.2,153.0};
  // G4double rInner[6]  = {291., 291., 2016.,  2016., 2044.,2044.};
  // //G4double rInner[6]  = {291., 2016., 2016.,  2016., 2044.,2044.};
  // G4double rOuter[6]  = {2070.,2070.,2070.,2070.,2070.,2070.};    
  // G4Polycone t1("aPcone",  
  // 		0., 2.*pi,       
  // 		6,zPlanes,rInner,rOuter); 
  //----------------------------------------------------------------
  // Complex Polycone Example from CMS DetectorDescription via GDML 
  //----------------------------------------------------------------    
    G4double zPlanes[6] = {-5541.,-3750.,-3750.,3750.,3750.,5541.};
    G4double rInner[6]  = {89.3,82.2452,1775.,1775.,82.2452,89.3};
    G4double rOuter[6]  = {2950.,2950.,2950.,2950.,2950.,2950.};    
    G4Polycone t1("aPcone",  
   		0., 2.*pi,       
   		6,zPlanes,rInner,rOuter); 
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
  
  G4cout <<" Check PolyCone had "<<n<<" inconsistencies"<< G4endl; 
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
  //-------------------------------------------
  // Original Polyhedra Test 
  //-------------------------------------------         
  // G4double zPlanes[5] = {-1., 10., 15., 25., 30.};
  // G4double rInner[5]  = {0., 5., 0.,  7., 1.};
  // G4double rOuter[5]  = {21., 6., 15., 15., 38.}; 
    
  ////   G4double z[10] = {30.,25.,15.,10.,-1.,-1.,10.,15.,25.,30.};
  ////   G4double r[10] = {1.,7.,0.,5.,0.,21.,6.,15.,15.,38.};
    
  // G4Polyhedra t1("aPhedra",  
  //		53.*deg, 163.*deg,         
  //    	8,5,zPlanes,rInner,rOuter); 
     
  ////   G4Polyhedra t1("aPhedra",  
  //// 		 53.*deg, 163.*deg,         
  //// 		 8, 10, r, z);
  //-------------------------------------------
  // More Complex Polyhedra with  Z1=Z2
  //-------------------------------------------  
  // G4double zPlanes[5] = {-1., 15., 15., 25., 30.};
  // G4double rInner[5]  = {0., 5., 0.,  7., 1.};
  // G4double rOuter[5]  = {21., 6., 15., 15., 38.}; 
  // G4Polyhedra t1("aPhedra",  
  //		53.*deg, 163.*deg,         
  //		8,5,zPlanes,rInner,rOuter); 
  //--------------------------------------------------------
  // Complex Polyhedra from CMS Detector Description via GDML
  //---------------------------------------------------------
  //
  // Example with numSide=1
  //  
  // 
  // G4double zPlanes[6] = {3240., 3704.6, 3750.22, 4460.79,4491.27, 5541.};
  // G4double rInner[6]  = { 1775., 1775.,1775.,2770.71, 2813.42,2813.42};
  // G4double rOuter[6]  = {1866.5,1866.5,1927.03, 2870., 2870., 2870.}; 
  // G4Polyhedra t1("aPhedra",  
  //		350.*deg, 20.*deg,         
  //	1,6,zPlanes,rInner,rOuter); 
  //
  // Example with startPhi<0 and phiTotal=twopi
  //  
     G4double zPlanes[8] = {3893.58, 3980.58, 3980.58, 4461.93,5167.08, 5167.08,5515.08,5541};
     G4double rInner[8]  = { 1712.1,1750.11,399.902,447.946,518.33,518.33,553.065,553.065};
     G4double rOuter[8]  = {1884.78,2000.23,2000.23,2639,2639,2459,2459,2459}; 
     
     G4Polyhedra t1("aPhedra",  
  		-10.*deg,360.*deg,         
  	18,8,zPlanes,rInner,rOuter); 
    
   
  G4Timer time;
  time.Start();
  
  for(i=0; i<N; i++)
  {  G4cout <<"I="<<i<<" stil to check = "<<N-i<<G4endl;
    
    point = t1.GetPointOnSurface();
      
    surf  = t1.Inside(point);
    if(surf != kSurface)
    { 
      n++; what = false;
      G4cout <<"Inconsistence="<<n<<" x "<<point.x()<<" y "<<point.y()<<" z "<<point.z()<<G4endl;
    }  
  }
  
  time.Stop();

  G4cout <<" Check Polyhedra had "<<n<<" inconsistencies"<< G4endl; 
  G4cout <<" Time taken was: "<<time.GetRealElapsed()<<" seconds."<<G4endl<<G4endl;
  return what;
}   

int main() 
{ 
  G4bool what;    
  //G4int N = 10000000;
    G4int N = 50000;
  
  G4cout <<G4endl;
  G4cout <<"********************************************************************************"<<G4endl;
  G4cout <<"**************** TEST GET POINT ON SURFACE METHOD ******************************"<<G4endl;
  G4cout <<"******************** FOR "<<N<<" RANDOM POINTS *********************************"<<G4endl;
  G4cout <<"********************************************************************************"<<G4endl;
   
  G4cout <<G4endl<<G4endl;

   what = checkPolycone(N);           
  // what = checkPolyhedra(N);      
  
  G4cout <<G4endl;
    
  G4cout <<"********************************************************************************"<<G4endl;
  G4cout <<"********************** END OF POINT-ON-SURFACE TEST ****************************"<<G4endl;
  G4cout <<"********************************************************************************"<<G4endl;  
  G4cout <<G4endl;

  return 0;     
}
       
   
 
