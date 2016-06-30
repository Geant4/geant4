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
// $Id: G4PolarizationHelper.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizationHelper
//
// Author:        Andreas Schaelicke
//
// Creation date: 12.08.2006
//
// Modifications:
//
// Class Description:
//
// Provides some basic polarization transformation routines.
//
#include "G4PolarizationHelper.hh"
#include "G4PhysicalConstants.hh"
#include "G4StokesVector.hh"


G4ThreeVector G4PolarizationHelper::GetFrame(const G4ThreeVector & mom1, const G4ThreeVector & mom2)
{
  G4ThreeVector normal = (mom1.cross(mom2)).unit();
  return normal;
  //  return 1./normal.mag()*normal;
}

G4ThreeVector G4PolarizationHelper::GetParticleFrameY(const G4ThreeVector &uZ)
{
  // compare also G4ThreeVector::rotateUz()

  if (uZ.x()==0. && uZ.y()==0.) {
    return G4ThreeVector(0.,1.,0.);
  }

  G4double invPerp = 1./std::sqrt(sqr(uZ.x())+sqr(uZ.y()));
  return G4ThreeVector(-uZ.y()*invPerp,uZ.x()*invPerp,0);
}

G4ThreeVector G4PolarizationHelper::GetParticleFrameX(const G4ThreeVector &uZ)
{
  // compare also G4ThreeVector::rotateUz()

  if (uZ.x()==0. && uZ.y()==0.) {
    if (uZ.z()>=0.) return G4ThreeVector(1.,0.,0.);
    return G4ThreeVector(-1.,0.,0.);
  }

  G4double perp    = std::sqrt(sqr(uZ.x())+sqr(uZ.y()));
  G4double invPerp = uZ.z()/perp;
  return G4ThreeVector(uZ.x()*invPerp,uZ.y()*invPerp,-perp);
}

G4ThreeVector G4PolarizationHelper::GetRandomFrame(const G4ThreeVector & mom1)
{
  G4double phi     =2.*pi*G4UniformRand();
  G4ThreeVector normal = std::cos(phi)*GetParticleFrameX(mom1) 
    + std::sin(phi)*G4PolarizationHelper::GetParticleFrameY(mom1);
  return normal;
}


G4ThreeVector G4PolarizationHelper::GetSpinInPRF(const G4ThreeVector &uZ, const G4ThreeVector & spin)
{
  // compare also G4ThreeVector::rotateUz()

  if (uZ.x()==0. && uZ.y()==0.) {
    if (uZ.z()>=0.) return spin;
    return G4ThreeVector(-spin.x(),spin.y(),-spin.z());
  }

  G4double perp    = std::sqrt(sqr(uZ.x())+sqr(uZ.y()));
  G4double invPerp = 1./perp;

  G4ThreeVector uX(uZ.x()*uZ.z()*invPerp,uZ.y()*uZ.z()*invPerp,-perp);
  G4ThreeVector uY(-uZ.y()*invPerp,uZ.x()*invPerp,0); 
  
  return G4ThreeVector(spin*uX,spin*uY,spin*uZ);
}

void G4PolarizationHelper::TestPolarizationTransformations()
{
  G4double theta=0.;
  G4cout<<"========================================\n\n";
  for (G4int i=0; i<=10; ++i) {
    theta=pi*i/10.;
    G4ThreeVector zAxis = G4ThreeVector(std::sin(theta),0.,std::cos(theta));
    if (i==5) zAxis = G4ThreeVector(1.,0.,0.);
    if (i==10) zAxis = G4ThreeVector(0.,0.,-1.);
    G4ThreeVector yAxis = GetParticleFrameY(zAxis);

    G4cout<<zAxis<<" "<<zAxis.mag()<<"\n";
    G4cout<<yAxis<<" "<<yAxis.mag()<<"\n";
    G4ThreeVector xAxis = yAxis.cross(zAxis);
    G4cout<<xAxis<<" "<<xAxis.mag()<<"\n\n";
  }

  G4cout<<"========================================\n\n";

  for (G4int i=0; i<=10; ++i) {
    theta=pi*i/10.;
    G4ThreeVector zAxis = G4ThreeVector(0.,std::sin(theta),std::cos(theta));
    if (i==5) zAxis = G4ThreeVector(0.,1.,0.);
    if (i==10) zAxis = G4ThreeVector(0.,0.,-1.);
    G4ThreeVector yAxis = GetParticleFrameY(zAxis);

    G4cout<<zAxis<<" "<<zAxis.mag()<<"\n";
    G4cout<<yAxis<<" "<<yAxis.mag()<<"\n";
    G4ThreeVector xAxis = yAxis.cross(zAxis);
    G4cout<<xAxis<<" "<<xAxis.mag()<<"\n\n";

    G4cout<<"spat : "<<xAxis*yAxis.cross(zAxis)<<"\n\n";
  }
  G4cout<<"========================================\n\n";
}

void G4PolarizationHelper::TestInteractionFrame()
{
  // check transformation procedure for polarisation transfer 
  // calculation in scattering processes
  //  a) transfer target polarisation in beam particle reference frame (PRF)
  //  b) calc correct asymmetry w.r.t. scattering plane
  //  c) determine incomming polarisation in interaction frame (IF)
  //  d) transfer outgoing polarisation from IF to PRF
  G4cout<<"========================================\n\n";

  G4double theta=0.;

  G4ThreeVector dir0=G4ThreeVector(0.,0.,1.);
  G4ThreeVector dir2=G4ThreeVector(std::sin(theta),0.,std::cos(theta));

  G4StokesVector pol0=G4ThreeVector(0.,0.,1.);
  G4StokesVector pol1=G4ThreeVector(0.,0.,1.); 

  pol1.rotateUz(dir0);

  G4cout<<"========================================\n\n";


}
