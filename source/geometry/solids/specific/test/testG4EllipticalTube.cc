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
// $Id: testG4EllipticalTube.cc
// GEANT4 tag $Name: 
// 
// testG4EllipticalTube
//
// Test file for class G4EllipticalTube derived from G4EllipticalCone
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4GeometryTolerance.hh"

#include "G4ThreeVector.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

const G4String OutputInside(const EInside a)
{
	switch(a) 
        {
		case kInside:  return "Inside"; 
		case kOutside: return "Outside";
		case kSurface: return "Surface";
	}
	return "????";
}


G4bool testG4EllipticalTube()
{

    G4ThreeVector pzero(0,0,0);
    G4ThreeVector pout ( 0*cm, 0*cm, 25*cm ) ;
    G4ThreeVector psurface(2.*mm,0*mm,0*cm);
    G4ThreeVector dir = pzero-pout ;
    dir *= 1/dir.mag();

    //For Normal 
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;//For Normal 
    
    
    G4EllipticalTube t1("Solid EllipticalTube #1",
			2*mm,       // xSemiAxis
			1*mm,       // ySemiAxis
			15*cm);      // zheight   
			
		           
    G4double dist = t1.DistanceToOut(pout,dir) ; 
    G4cout << "*********** Testing DistanceToOut method *************** "<<G4endl;
    G4cout << "Distance = " << dist << G4endl <<G4endl;  
    dist = t1.DistanceToOut(psurface,G4ThreeVector(1,0,0),calcNorm,pgoodNorm,pNorm);
    G4cout << "Distance = " << dist << G4endl ; 
    G4cout << "with Normal = "<< (*pNorm) << G4endl << G4endl ; 
    // Check Inside 
 
    G4cout << "************ Test Inside(p) ****************" << G4endl ;
    G4cout << "pzero :  " << t1.Inside(pzero) << G4endl <<G4endl ; 
    G4cout << "psurface:" << t1.Inside(psurface) << G4endl <<G4endl ;
    //test the name
    G4cout << "The name is : " << t1.GetName() << G4endl ;  

    // testing the volume 
 
    G4double volume = t1.GetCubicVolume() ; 
    G4cout << G4endl ;
    G4cout << "Solid EllipticalTube #1 has Volume = " << volume / cm / cm /cm << " cm^3"  
	   << G4endl << G4endl; 

    return true; 
}

// 
//  This test generates a random point on the surface of the solid and  
//  checks the distance from a point outside in the direction of the 
//  line between the two points 
//
G4bool testDistanceToIn() 
{    
  G4EllipticalTube t1("Solid EllipticalTube #1", 
		      2*cm,       // xSemiAxis
		      1*cm,       // ySemiAxis  
		      15*cm);      // zheight
		       
  
  G4int N = 10000;
  G4int n = 0;

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  for(G4int i=0; i<N; i++)
  {
    G4ThreeVector point = t1.GetPointOnSurface();
    point.setX(std::fabs(point.x()));
    point.setY(std::fabs(point.y()));
    
    G4ThreeVector out (100*cm, 100*cm, 100*cm);
    
    G4ThreeVector dir = point - out;
    G4double dist2 = dir.mag();
    dir /= dist2;
    
    G4double dist1 = t1.DistanceToIn(point,dir);
    G4double diff = std::fabs(dist1 - dist2);
    
    if(diff < 2.*kCarTolerance)
      n++;      
  }

  G4cout <<" ************   For testG4EllipticalTube   ******************"<<G4endl<<G4endl;
  G4cout <<" Number of inconsistencies for testDistanceToIn was: "<< n <<" ..."<<G4endl
	 <<" ... For a total of "<<N<<" trials."<< G4endl <<G4endl;

  return true;
}
G4bool testDistanceToOut() 
{    
  G4EllipticalTube t1("Solid EllipticalTube #1", 
		      2*cm,       // xSemiAxis
		      1*cm,       // ySemiAxis  
		      15*cm);      // zheight 
  
  G4int N = 10000;
  G4int n = 0;
  G4int nn= 0;

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  //For Normal 
    G4ThreeVector *pNorm,norm,SurfNormal;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;
    
    pNorm=&norm;
    pgoodNorm=&goodNorm;
 //

  for(G4int i=0; i<N; i++)
  {
    G4ThreeVector point = t1.GetPointOnSurface();
    point.setX(std::fabs(point.x()));
    point.setY(std::fabs(point.y()));
    
    G4ThreeVector out (0*cm, 0*cm, 0*cm);
    
    G4ThreeVector dir = point - out;
    G4double dist2 = dir.mag();
    dir /= dist2;
    
    G4double dist1 = t1.DistanceToOut(point,dir,calcNorm,pgoodNorm,pNorm);
    G4double diff = std::fabs(dist1 - dist2);
    G4double difNorm=((*pNorm)-t1.SurfaceNormal(point)).mag();
    
    if(difNorm > 2.*kCarTolerance){
    if(std::fabs(std::fabs(point.z())-15.0*cm)>kCarTolerance)  nn++;   
      G4cout<<"point="<<point<<" norm="<< (*pNorm)<<" surfNorm="<<t1.SurfaceNormal(point)<<G4endl;
     
    } 

    if(diff < 2.*kCarTolerance)
      n++;      
  }

  G4cout <<" ************   For testG4EllipticalTube   ******************"<<G4endl<<G4endl;
  G4cout <<" Number of inconsistencies for testDistanceToOut was: "<< n <<" ..."<<G4endl
	 <<" ... For a total of "<<N<<" trials."<< G4endl <<G4endl;
  G4cout <<" Number of differences for test Normal was "<< nn << " ...."<<G4endl <<G4endl;

  return true;
}

 
int main() 
{ 

  G4cout << G4endl;
  G4cout << "*********************************************************************" <<G4endl;
  G4cout << "****************** UNIT TEST FOR ELLIPTICAL TUBE ********************" <<G4endl;
  G4cout << "*********************************************************************" <<G4endl;
  G4cout << G4endl;

  // temporary test
  G4ThreeVector Spoint ;
  G4double dist ;

  G4EllipticalTube t1("Solid EllipticalTube #1",
		      0.5*mm,       // xSemiAxis
		      1*mm,       // ySemiAxis
		      40*cm);      // zheight   
  

  EInside side ;  
  for ( G4int i = 0 ; i < 3 ; i++ ) {
    // G4cout << "Event " << i << G4endl << G4endl ;
    Spoint = t1.GetPointOnSurface() ;
    side = t1.Inside(Spoint) ;
    dist = t1.DistanceToIn(Spoint, -Spoint/Spoint.mag()) ;
    G4cout << "Spoint " << Spoint << " " <<  dist  << " " << side  << G4endl ;
  }

#ifdef NDEBUG 
  G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif 
  assert(testG4EllipticalTube()); 
    
  G4bool what;
  what = testDistanceToIn();
  what = testDistanceToOut();
  G4cout << G4endl;
  G4cout << "*********************************************************************" <<G4endl;
  G4cout << "******************* END OF TEST - THANK YOU!!! **********************" <<G4endl;
  G4cout << "*********************************************************************" <<G4endl;
  G4cout << G4endl;

  return 0;
}
