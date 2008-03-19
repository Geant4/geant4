



// test for random direction unit vector algorithm
// author: V. Grichine
// dased on discussions and suggestions of G. cosmo
//
//
// History: 
//
//  19.03.08 first impementation





#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4UnitsTable.hh"
#include "G4Timer.hh"


#include <iomanip>







///////////////////////////////////////////////////////////////////////
//
// Random algorithm from G4OpDundaryProcess

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
  G4double cosTheta = 2*G4UniformRand()-1.;
  G4double sinTheta2 = 1. - cosTheta*cosTheta;
  if (sinTheta2 < 0.) sinTheta2 = 0.;
  G4double sinTheta = std::sqrt(sinTheta2); 
  G4double phi = twopi*G4UniformRand();
  return G4ThreeVector(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta).unit(); 
}

//////////////////////////////////////////////////////////////////////////////
//
// Test main program

int main()
{
  G4int i, iMax = 20;  // , k, kMax;

  
    
  G4Timer timer;
  iMax = 1000000;

  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = IsotropicCubeRand();  
  }
  timer.Stop();
  G4cout<<"Total time of volume "<<iMax<<" calls = "<<timer.GetUserElapsed()<<" s"<<G4endl<<G4endl;

  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = IsotropicSphereRand();  
  }
  timer.Stop();
  G4cout<<"Total time of surface "<<iMax<<" calls = "<<timer.GetUserElapsed()<<" s"<<G4endl<<G4endl;
  
  timer.Start();
  for( i = 0; i < iMax; i++ )
  {
    G4ThreeVector isoVectr = G4RandomDirection();  
  }
  timer.Stop();
  G4cout<<"Total time of G4RandomDirection() "<<iMax<<" calls = "<<timer.GetUserElapsed()
        <<" s"<<G4endl<<G4endl;

  


  iMax = 1000000;
  G4int j, jMax = 100;
  G4int cosThetaDistr[100], phi[100];

  for( j = 0; j < jMax; j++ )
  {
    cosThetaDistr[j] = 0;  
    phi[j]           = 0;  
  }    
  G4double xyPlane, phiNow, cosThetaNow, cosThetaTmp, phiTmp;
  
  for( i = 0; i < iMax; i++ )
  {
    // G4ThreeVector isoVectr = IsotropicCubeRand();
    G4ThreeVector isoVectr = IsotropicSphereRand();
    // G4ThreeVector isoVectr = G4RandomDirection();  

    xyPlane = std::sqrt( isoVectr.x()*isoVectr.x() + isoVectr.y()*isoVectr.y() );

    if ( xyPlane ) 
    {
      phiNow  = std::atan2(isoVectr.y(),isoVectr.x());
      phiNow += pi;                                     // 0-twopi range

      cosThetaNow = isoVectr.z();
    }
    else
    {
      if ( isoVectr.z() >= 0. ) cosThetaNow = 1.;
      else                      cosThetaNow = -1.;
      phiNow = twopi*G4UniformRand();      
    }
    for( j = 0; j < jMax; j++ )
    {
      cosThetaTmp = -1. + 2.*j/jMax;
      if( cosThetaTmp >= cosThetaNow )
      {
        cosThetaDistr[j]++;
        break;
      }
    }    
    for( j = 0; j < jMax; j++ )
    {
      phiTmp      = twopi*j/jMax;
      if( phiTmp >= phiNow )
      {
        phi[j]++;
        break;
      }
    }     
  }
  G4cout << G4endl;
  G4cout <<"cosThetaTmp"<<"\t"<<"cosThetaDistr[j]"<<"\t"<<"phi/degree"<<"\t"<<"phi[j]"<<G4endl;
  G4cout << G4endl;

  for( j = 0; j < jMax; j++ )
  {
    cosThetaTmp = -1. + 2.*j/jMax;
    phiTmp      = twopi*j/jMax;
    G4cout <<cosThetaTmp<<"\t"<<cosThetaDistr[j]<<"\t"<<phiTmp/degree<<"\t"<<phi[j]<<G4endl;
  }    



  return 1 ;
}



