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
// $Id: G4RandomDirectionTest.cc,v 1.1 2008-03-19 16:53:55 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Test for random direction unit vector algorithm
//
// Author: V. Grichine
//
//
// History: 
//
// 19.03.08 - First implementation
// ----------------------------------------------------------------------

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include <iomanip>

#include "globals.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4UnitsTable.hh"
#include "G4Timer.hh"

///////////////////////////////////////////////////////////////////////
//
// Random algorithm from Geant3 RAND3

G4ThreeVector IsotropicCubeRand() 
{
  /* Returns a random isotropic unit vector. */

  G4ThreeVector vect;
  G4double len2;

  do {

    vect.setX(G4UniformRand() - 0.5);
    vect.setY(G4UniformRand() - 0.5);
    vect.setZ(G4UniformRand() - 0.5);

    len2 = vect.mag2();

  } while (len2 < 0.01 || len2 > 0.25);

  return vect.unit();
}

/////////////////////////////////////////////////////////////////////////////
//
// Random distribution over unit radius sphere

G4ThreeVector IsotropicSphereRand() 
{
  G4double cosTheta  = 2*G4UniformRand()-1.;
  G4double sinTheta2 = 1. - cosTheta*cosTheta;
  if( sinTheta2 < 0.)  sinTheta2 = 0.;
  G4double sinTheta  = std::sqrt(sinTheta2); 
  G4double phi       = twopi*G4UniformRand();
  return G4ThreeVector(sinTheta*std::cos(phi),
                       sinTheta*std::sin(phi), cosTheta).unit(); 
}

///////////////////////////////////////////////////////////////////////////////
//
// 8 quadrants algorithm

G4ThreeVector MKRandomDirection()
{
  // Randomization in one of 8 Quadrants (x>0, y>0, z>0)
  //
  G4double x=G4UniformRand(), y=G4UniformRand(), z=G4UniformRand();
  G4double r2= x*x+y*y+z*z;
  while(r2>1.||r2<.000001)
  {
    x = G4UniformRand(); y = G4UniformRand(); z = G4UniformRand();
    r2=x*x+y*y+z*z;
  }
  G4double r=std::sqrt(r2), quad=G4UniformRand();

  if(quad>0.5)
  {
    if(quad>0.75)
    {
      if(quad>0.875)    return G4ThreeVector(-x/r,-y/r,-z/r);
      else              return G4ThreeVector(-x/r,-y/r, z/r);
    }
    else
    {
      if(quad>0.625)    return G4ThreeVector(-x/r, y/r,-z/r);
      else              return G4ThreeVector(-x/r, y/r, z/r);
    }
  }
  else
  {
    if(quad>0.25)
    {
      if(quad>0.375)    return G4ThreeVector( x/r,-y/r,-z/r);
      else              return G4ThreeVector( x/r,-y/r, z/r);
    }
    else if(quad>0.125) return G4ThreeVector( x/r, y/r,-z/r);
  }
  return                       G4ThreeVector( x/r, y/r, z/r);
}


//////////////////////////////////////////////////////////////////////////////
//
// Test main program

int main()
{
  G4int i, iMax = 20;

  G4Timer timer;

  iMax = 1000000;

  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = IsotropicCubeRand();  
  }
  timer.Stop();
  G4cout<<"Total time of volume "<<iMax<<" calls = "
        <<timer.GetUserElapsed()<<" s"<<G4endl<<G4endl;

  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = IsotropicSphereRand();  
  }
  timer.Stop();
  G4cout<<"Total time of surface "<<iMax<<" calls = "
        <<timer.GetUserElapsed()<<" s"<<G4endl<<G4endl;
  
  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = G4RandomDirection();  
  }
  timer.Stop();
  G4cout<<"Total time of G4RandomDirection() "<<iMax
        <<" calls = "<<timer.GetUserElapsed()
        <<" s"<<G4endl<<G4endl;

  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = MKRandomDirection();  
  }
  timer.Stop();
  G4cout<<"Total time of MKRandomDirection() "<<iMax
        <<" calls = "<<timer.GetUserElapsed()
        <<" s"<<G4endl<<G4endl;

  iMax = 1000000;
  G4int j, jMax = 100;
  G4double cosThetaDistr[100], phi[100];

  for( j = 0; j < jMax; j++ )
  {
    cosThetaDistr[j] = 0.;  
    phi[j]           = 0.;  
  }    
  G4double xyPlane, phiNow, cosThetaNow, cosThetaTmp, phiTmp;
  
  for( i = 0; i < iMax; i++ )
  {
    // G4ThreeVector isoVectr = IsotropicCubeRand();
    // G4ThreeVector isoVectr = IsotropicSphereRand();
    G4ThreeVector isoVectr = G4RandomDirection();  

    xyPlane = std::sqrt( isoVectr.x()*isoVectr.x()
                       + isoVectr.y()*isoVectr.y() );
    if ( xyPlane ) 
    {
      phiNow  = std::atan2(isoVectr.y(),isoVectr.x());
      phiNow += pi;                                     // 0-twopi range

      cosThetaNow = isoVectr.z();
    }
    else
    {
      if ( isoVectr.z() >= 0. ) cosThetaNow =  1.;
      else                      cosThetaNow = -1.;
      phiNow = twopi*G4UniformRand();      
    }
    for( j = 0; j < jMax; j++ )
    {
      cosThetaTmp = -1. + 2.*j/jMax;

      if( cosThetaTmp >= cosThetaNow )
      {
        cosThetaDistr[j] += 1.;
        break;
      }
    }    
    for( j = 0; j < jMax; j++ )
    {
      phiTmp      = twopi*j/jMax;
      if( phiTmp >= phiNow )
      {
        phi[j] += 1.;
        break;
      }
    }     
  }
  G4cout << G4endl;
  G4cout <<"cosThetaTmp"<<"\t"<<"cosThetaDistr[j]"<<"\t"
         <<"phi/degree"<<"\t"<<"phi[j]"<<G4endl;
  G4cout << G4endl;

  for( j = 0; j < jMax; j++ )
  {
    cosThetaTmp = -1. + 2.*j/jMax;
    phiTmp      = twopi*j/jMax;
    G4cout <<cosThetaTmp<<"\t"<<cosThetaDistr[j]<<"\t"
           <<phiTmp/degree<<"\t"<<phi[j]<<G4endl;
  }
    
  G4double mean=0., rms2=0., rms, relDelta;

  for( j = 0; j < jMax; j++ )
  {
    mean += cosThetaDistr[j];
  }
  mean /= jMax;

  for( j = 0; j < jMax; j++ )
  {
    rms2 += (mean-cosThetaDistr[j])*(mean-cosThetaDistr[j]);
  }
  rms2 /= jMax-1;

  rms = std::sqrt(rms2);

  relDelta = rms/mean;
  G4cout << G4endl;
  G4cout << "meanCosThetaDistr = "<<mean<<" +- "<<rms
         <<" with relative error = "<<relDelta <<G4endl;

  mean = 0;
  rms2 = 0.;

  for( j = 0; j < jMax; j++ )
  {
    mean += phi[j];
  }
  mean /= jMax;

  for( j = 0; j < jMax; j++ )
  {
    rms2 += (mean-phi[j])*(mean-phi[j]);
  }
  rms2 /= jMax-1;

  rms = std::sqrt(rms2);

  relDelta = rms/mean;
  G4cout << G4endl;
  G4cout << "meanPhiDistr = "<<mean<<" +- "<<rms
         <<" with relative error = "<<relDelta <<G4endl;

  return 1 ;
}
