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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction2.cc
/// \brief Implementation of the PrimaryGeneratorAction2 class
//
//
// $Id: PrimaryGeneratorAction2.cc 68024 2013-03-13 13:42:01Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PrimaryGeneratorAction2.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction2::PrimaryGeneratorAction2(G4ParticleGun* gun)
: particleGun(gun)
{    
  // energy distribution
  //
  InitFunction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction2::~PrimaryGeneratorAction2()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction2::GeneratePrimaries(G4Event* anEvent)
{
  //cosAlpha uniform in [cos(0), cos(pi)]
  G4double cosAlpha = 1. - 2*G4UniformRand();
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi      = twopi*G4UniformRand();  //psi uniform in [0, 2*pi]  
  G4ThreeVector dir(sinAlpha*std::cos(psi),sinAlpha*std::sin(psi),cosAlpha);

  particleGun->SetParticleMomentumDirection(dir);
  
  //set energy from a tabulated distribution
  //
  //G4double energy = RejectAccept();
  G4double energy = InverseCumul();  
  particleGun->SetParticleEnergy(energy);    

  //create vertex
  //   
  particleGun->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction2::InitFunction()
{
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
  //
  nPoints = 16;
  const G4double xx[] = 
    { 37*keV, 39*keV, 45*keV,  51*keV,  57*keV,  69*keV,  71*keV,  75*keV, 
      83*keV, 91*keV, 97*keV, 107*keV, 125*keV, 145*keV, 159*keV, 160*keV }; 
      
  const G4double ff[] =
    { 0.000,  0.077,  0.380,  2.044, 5.535, 15.077, 12.443, 14.766,
     17.644, 18.518, 17.772, 14.776, 8.372,  3.217,  0.194,  0.000 };
  
  //copy arrays in std::vector and compute fMax
  //
  x.resize(nPoints); f.resize(nPoints);
  fMax = 0.;
  for (G4int j=0; j<nPoints; j++) {
    x[j] = xx[j]; f[j] = ff[j];
    if (fMax < f[j]) fMax = f[j];
  };
     
  //compute slopes
  //
  a.resize(nPoints);
  for (G4int j=0; j<nPoints-1; j++) { 
    a[j] = (f[j+1] - f[j])/(x[j+1] - x[j]);
  };
  
  //compute cumulative function
  //
  Fc.resize(nPoints);  
  Fc[0] = 0.;
  for (G4int j=1; j<nPoints; j++) {
    Fc[j] = Fc[j-1] + 0.5*(f[j] + f[j-1])*(x[j] - x[j-1]);
  };     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction2::RejectAccept()
{
  // tabulated function 
  // f is assumed positive, linear per segment, continuous
  //  
  G4double x_rndm = 0., y_rndm = 0., f_inter = -1.;
  
  while (y_rndm > f_inter) {
    //choose a point randomly
    x_rndm = x[0] + G4UniformRand()*(x[nPoints-1] - x[0]);
    y_rndm = G4UniformRand()*fMax;
    //find bin
    G4int j = nPoints-2;
    while ((x[j] > x_rndm) && (j > 0)) j--;
    //compute f(x_rndm) by linear interpolation
    f_inter = f[j] + a[j]*(x_rndm - x[j]);
  };
  return x_rndm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction2::InverseCumul()
{
  // tabulated function
  // f is assumed positive, linear per segment, continuous 
  // --> cumulative function is second order polynomial
  
  //choose y randomly
  G4double y_rndm = G4UniformRand()*Fc[nPoints-1];
  //find bin
  G4int j = nPoints-2;
  while ((Fc[j] > y_rndm) && (j > 0)) j--;
  //y_rndm --> x_rndm :  Fc(x) is second order polynomial
  G4double x_rndm = x[j];
  G4double aa = a[j];
  if (aa != 0.) {
    G4double b = f[j]/aa, c = 2*(y_rndm - Fc[j])/aa;
    G4double delta = b*b + c;
    G4int sign = 1; if (aa < 0.) sign = -1;
    x_rndm += sign*std::sqrt(delta) - b;    
  } else if (f[j] > 0.) {
    x_rndm += (y_rndm - Fc[j])/f[j];
  };
  return x_rndm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
