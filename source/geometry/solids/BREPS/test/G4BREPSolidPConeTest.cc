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
//////////////////////////////////////////////////////////////////////////
// $Id: G4BREPSolidPConeTest.cc,v 1.12 2002-01-28 16:16:23 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//////////////////////////////////////////////////////////////////////////
//
//
// BREP solid test, create by L. Broglia, 20/10/98
// modification of old G4Gerep test
//


#include "G4Timer.hh"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "G4ios.hh" 
#include "G4BREPSolid.hh"
#include "G4BREPSolidPCone.hh"

#include "g4std/fstream"
#include "g4std/iomanip"

void checkNormal( G4BREPSolid* solid, G4ThreeVector& position )
{
  G4ThreeVector normal = solid->SurfaceNormal( position );
  
  G4cout << G4endl
         << "\t\tSurface normal"                 << G4endl
         << "\t\t--------------"                 << G4endl
         << "\t\tposition x="                    << position.x()
         << " y="                                << position.y()
         << " z="                                << position.z() << G4endl
         << "\t\tis normX="                      << normal.x()
         << " normY="                            << normal.y()
         << " normZ="                            << normal.z()   << G4endl;
}

void checkSurfInOut( G4BREPSolid* solid, G4ThreeVector& position, G4ThreeVector& direction )
{
  G4double d;
  
  G4cout << G4endl
         << "\tTest of In & Out " << solid->GetName()   << G4endl
         << "\t--------------------------------------" << G4endl;
  
  d = solid->DistanceToIn( position, direction );
  G4cout << "\tDistance to in="          << d;
  d = solid->DistanceToIn( position );
  G4cout << "\tClosest distance to in="  << d          <<G4endl;
    
  d = solid->DistanceToOut( position, direction );
  G4cout << "\tDistance to out="         << d;
  d = solid->DistanceToOut( position );
  G4cout << "\tClosest distance to out=" << d          << G4endl;
}

void checkSolid( const G4String& where, G4BREPSolid* solid, G4ThreeVector& position, G4ThreeVector& direction, G4double estIn = 0.0, G4double estOut = 0.0 )
{
  G4double d;
  
  EInside in = solid->Inside( position );
  
  G4cout << G4endl
         << where << " of " << solid->GetName()                   << G4endl
         << "--------------------------------" << G4endl
         << "estimation to in  =" << estIn     << G4endl
         << "estimation to out =" << estOut    << G4endl
         << "direction x=" << direction.x()
         << "  y=" << direction.y()
         << "  z=" << direction.z()            << G4endl
         << "position x=" << position.x()
         << "  y=" << position.y()
         << "  z=" << position.z();

  if( in == kInside ) {
    G4cout <<" is inside";
    d = solid->DistanceToOut( position, direction );
    G4cout<<"  distance to out="<<d;
    d = solid->DistanceToOut( position );
    G4cout<<"  closest distance to out="<<d<<G4endl;
  } else if( in == kOutside ) {
    G4cout <<" is outside";
    d = solid->DistanceToIn( position, direction );
    G4cout<<"  distance to in="<<d;
    d = solid->DistanceToIn( position );
    G4cout<<"  closest distance to in="<<d<<G4endl;
  } else {
    G4cout <<" is on the surface"<<G4endl;
    checkSurfInOut( solid, position, direction );
    checkNormal( solid, position );
  }
}

int main(G4int argc, char **argv)
{
  const G4int noZplanes= 8; 
  
  G4double RMINVec[noZplanes]  = {  30,  30,   0,  0,   0,   0,  40,  40 };
  G4double RMAXVec[noZplanes]  = {  70,  70,  70, 40,  40,  80,  80,  60 };
  G4double Z_Values[noZplanes] = { -20, -10, -10,  0,  10,  20,  30,  40 };

  G4double start_angle= 0.0;
  G4double opening_angle= 2*pi;

  G4double zstart= Z_Values[0]; 

  G4cout << "\n=======     PCon test      ========";

  G4BREPSolidPCone *MyPCone = new G4BREPSolidPCone ("MyPCone",
						    start_angle,
						    opening_angle,
						    noZplanes,
						    zstart,
						    Z_Values,
						    RMINVec,
						    RMAXVec );
  
  G4cout << "\n\nPCone created ! "<<G4endl;
  G4cout << "Variety is G4BREPSolidPolycone"<<G4endl;

  G4cout << "Its parameters are: "<<G4endl;

  ///////////////////////////////////////////////////
  // Temporary
  for (G4int x = 0; x < noZplanes; x++)
  {
    G4cout <<    " Z[" << x << "]=" << G4std::setw(5) << Z_Values[x];
    G4cout << " Rmin[" << x << "]=" << G4std::setw(5) << RMINVec[x];
    G4cout << " Rmax[" << x << "]=" << G4std::setw(5) << RMAXVec[x]<<G4endl;
  }

  G4cout<<" start   angle ="<<start_angle<<G4endl;
  G4cout<<" opening angle ="<<opening_angle<<G4endl;
  G4cout<<" zstart =" << zstart << G4endl;


  // -> Check methods :
  //  - Inside
  //  - DistanceToIn
  //  - DistanceToOut

  
  EInside in;
  
  G4cout<<"\n\n==================================================";
  G4ThreeVector  pt(0, -100, 24);
  G4double y; 
  for (y = -100; y<=100; y+=10)
  {
    pt.setY(y);
    in = MyPCone->Inside(pt);
    
    G4cout << "\nx=" << pt.x() << "  y=" << pt.y() << "  z=" << pt.z();
    
    if( in == kInside )
      G4cout <<" is inside";
    else
      if( in == kOutside )
	G4cout <<" is outside";
      else
	G4cout <<" is on the surface";
  }

  G4cout<<"\n\n==================================================";
  G4ThreeVector  start( 0, 0, -30);
  G4ThreeVector dir(1, 1, 0);
  G4double   d;
  
  G4cout<<"\nPdep is (0, 0, z)";
  G4cout<<"\nDir is (1, 1, 0)\n";

  G4double z;
  for(z=-30; z<=50; z+=5)
  {
    start.setZ(z);

    in = MyPCone->Inside(start);
    G4cout<< "x=" << start.x() << "  y=" << start.y() << "  z=" << start.z();
    
    if( in == kInside )
    {
      G4cout <<" is inside";

      d = MyPCone->DistanceToOut(start, dir);
      G4cout<<"  distance to out="<<d;
      d = MyPCone->DistanceToOut(start);
      G4cout<<"  closest distance to out="<<d<<G4endl;
    }
    else if( in == kOutside ) 
    {
      G4cout <<" is outside";

      d = MyPCone->DistanceToIn(start, dir);
      G4cout<<"  distance to in="<<d;
      d = MyPCone->DistanceToIn(start);
      G4cout<<"  closest distance to in="<<d<<G4endl;
    }
    else
      G4cout <<" is on the surface"<<G4endl;

  }

  G4cout<<"\n\n==================================================";
  G4ThreeVector  start2( 0, -100, -30);
  G4ThreeVector dir2(0, 1, 0);
  G4double   d2;

  G4cout<<"\nPdep is (0, -100, z)";
  G4cout<<"\nDir is (0, 1, 0)\n";

  for(z=-30; z<=50; z+=5)
  {
    G4cout<<"  z="<<z;
    start2.setZ(z);
    d2 = MyPCone->DistanceToIn(start2, dir2);
    G4cout<<"  distance to in="<<d2;
    d2 = MyPCone->DistanceToIn(start2);
    G4cout<<"  closest distance to in="<<d2<<G4endl;
  }

  G4cout<<"\n\n==================================================";
  G4ThreeVector  start3( 0, 0, -50);
  G4ThreeVector dir3(0, 0, 1);
  G4double   d3;

  G4cout<<"\nPdep is (0, y, -50)";
  G4cout<<"\nDir is (0, 0, 1)\n";

  for(y=-0; y<=90; y+=5)
  {
    G4cout<<"  y="<<y;
    start3.setY(y);
    d3 = MyPCone->DistanceToIn(start3, dir3);
    G4cout<<"  distance to in="<<d3<<G4endl;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // Test due to the bugfix No. 320
  // Solid setup (only profile of its upper half is shown):
  //
  //     +----+
  //  +--+    +--+
  //  +--+    +--+
  //     +----+
  //
  // -------------- Z axis
  //
 	G4double zValTOB[6]     = { -2800*mm, -1210*mm, -1210*mm, 1210*mm, 1210*mm, 2800*mm };
 	G4double rMinTOB[6]     = { 1160*mm,   1160*mm,  540*mm,  540*mm,  1160*mm, 1160*mm };
 	G4double rMaxTOB[6]     = { 1170*mm,   1170*mm,  1182*mm, 1182*mm, 1170*mm, 1170*mm };
 	G4BREPSolidPCone * cmsTOBSolid = new G4BREPSolidPCone( "BREPPolyconeTOBSolid",
                                                         0.0*rad, twopi*rad,
                                                  			 6, zValTOB[0], zValTOB, rMinTOB, 
                                                         rMaxTOB );
  // The critical gun settings
  G4ThreeVector gunPosition( 904.05833958335438*mm, -761.44764667689935*mm, -382.29360686839397*mm );
  G4ThreeVector gunDirection( -0.32132043849994285, 0.11159991009963287, 0.94037154139624957 );

  checkSolid( "CMS case", cmsTOBSolid, gunPosition, gunDirection );  

  G4ThreeVector gp( 541*mm, 0*mm, 1210*mm );
  G4ThreeVector gd( -0.32132043849994285, 0.11159991009963287, 0.94037154139624957 );
  
  checkSolid( "surface", cmsTOBSolid, gp, gd );  
  
  G4ThreeVector gpedge( 540*mm, 540*mm, -1210*mm );
  G4ThreeVector gdedge( -0.32132043849994285, 0.11159991009963287, 0.94037154139624957 );
  
  checkSolid( "edge", cmsTOBSolid, gpedge, gdedge );  
  
  G4ThreeVector gpin( 540*mm, 540*mm, -1209*mm );
  G4ThreeVector gdin( -0.32132043849994285, 0.11159991009963287, 0.94037154139624957 );
  
  checkSolid( "inside", cmsTOBSolid, gpin, gdin );  

 	G4double zVal1[6]     = { -200*mm,  -100*mm, -100*mm, 100*mm, 100*mm, 200*mm };
 	G4double rMin1[6]     = {   50*mm,    50*mm,   25*mm,  25*mm,  50*mm,  50*mm };
// 	G4double rMin1[6]     = {   50*mm,    50*mm,   50*mm,  50*mm,  50*mm,  50*mm };
 	G4double rMax1[6]     = {   50*mm,   100*mm,  200*mm, 200*mm, 100*mm,  50*mm };
 	G4BREPSolidPCone* s1  = new G4BREPSolidPCone( "s1", 0.0*rad, twopi*rad, 6, zVal1[0], zVal1, rMin1, rMax1 );
    
  G4ThreeVector gps1in( 150*mm, 0*mm, 0*mm );
  G4ThreeVector gds1in( 0, 0, 1 );
  G4ThreeVector gps1out( 0*mm, 0*mm, 0*mm );
  G4ThreeVector gds1out( 1, 0, 0 );
  G4ThreeVector gps1edge( 50*mm, 0*mm, -200*mm );
  G4ThreeVector gds1edge( 0, 0, 1 );
  G4ThreeVector gps1surfout( 200*mm, 0*mm, 0*mm );
  G4ThreeVector gds1surfout( -1, 0, 0 );
  G4ThreeVector gps1surfin( 0*mm, 25*mm, 0*mm );
  G4ThreeVector gds1surfin( 0, 1, 0 );
  G4ThreeVector gps1plansurfleftout( 150*mm, 0*mm, -100*mm );
  G4ThreeVector gps1plansurfleftin( 30*mm, 0*mm, -100*mm );
  G4ThreeVector gps1plansurfleftfict( 0*mm, 0*mm, -100*mm );
  G4ThreeVector gds1plansurfleft( 0, 0, 1 );
  G4ThreeVector gps1plansurfright( 150*mm, 0*mm, 100*mm );
  G4ThreeVector gds1plansurfright( 0, 0, -1 );
  
  checkSolid( "inside", s1, gps1in, gds1in );  
  checkSolid( "outside", s1, gps1out, gds1out );  
  checkSolid( "edge", s1, gps1edge, gds1edge, 0, 400 );  
  checkSolid( "outer surface", s1, gps1surfout, gds1surfout, 0, 100 );  
  checkSolid( "inner surface", s1, gps1surfin, gds1surfin, 0, 175 );  
  checkSolid( "left outer planar surface", s1, gps1plansurfleftout, gds1plansurfleft, 0, 200 );  
  checkSolid( "left outer planar surface inv", s1, gps1plansurfleftout, gds1plansurfright, kInfinity, 0 );  
  checkSolid( "fictious left planar surface", s1, gps1plansurfleftfict, gds1plansurfleft, kInfinity, 25 );  
  checkSolid( "fictious left planar surface inv", s1, gps1plansurfleftfict, gds1plansurfright, kInfinity, 25 );  
  checkSolid( "left inner planar surface", s1, gps1plansurfleftin, gds1plansurfleft, 0, 200 );  
  checkSolid( "left inner planar surface inv", s1, gps1plansurfleftin, gds1plansurfright, kInfinity, 0 );  
  checkSolid( "right planar surface", s1, gps1plansurfright, gds1plansurfright, 0, 200 );  
  checkSolid( "right planar surface inv", s1, gps1plansurfright, gds1plansurfleft, 0 );  

  return EXIT_SUCCESS;
}

