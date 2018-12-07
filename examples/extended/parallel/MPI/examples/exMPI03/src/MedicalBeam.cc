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
/// @file MedicalBeam.cc
/// @brief Define beam profile as primary generator

#include <cmath>
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryVertex.hh"
#include "Randomize.hh"
#include "MedicalBeam.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MedicalBeam::MedicalBeam()
  : fparticle(0), fkineticE(1.*MeV), fsourcePosition(G4ThreeVector()),
    fSSD(1.*m), ffieldShape(MedicalBeam::kSQUARE), ffieldR(10.*cm)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  fparticle = particleTable-> FindParticle("proton");

  fkineticE = 200.*MeV;
  fsourcePosition = G4ThreeVector(0.,0.,-125.*cm);
  fSSD = 100.*cm;
  ffieldXY[0] = ffieldXY[1] = 5.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MedicalBeam::~MedicalBeam()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector MedicalBeam::GenerateBeamDirection() const
{
  // uniform distribution in a limitted solid angle
  G4double dr;
  if ( ffieldShape == MedicalBeam::kSQUARE ) {
    dr = std::sqrt(sqr(ffieldXY[0]/2.) + sqr(ffieldXY[1]/2.));
  } else {
    dr = ffieldR;
  }

  G4double sin0 = dr/fSSD;
  G4double cos0 = std::sqrt(1.-sqr(sin0));

  G4double dcos = 0.;
  G4double dsin, dphi, z;

  G4double x = DBL_MAX;
  G4double y = DBL_MAX;

  G4double xmax, ymax;
  if ( ffieldShape == MedicalBeam::kSQUARE ) {
    xmax = ffieldXY[0]/2./fSSD;
    ymax = ffieldXY[1]/2./fSSD;
  } else {
    xmax = ymax = DBL_MAX-1.;
  }

  while(! (std::abs(x) < xmax && std::abs(y) < ymax) ) {
    dcos = G4RandFlat::shoot(cos0, 1.);
    dsin = std::sqrt(1.-sqr(dcos));
    dphi = G4RandFlat::shoot(0., twopi);

    x = std::cos(dphi)*dsin*dcos;
    y = std::sin(dphi)*dsin*dcos;
  }
  z = dcos;

  return G4ThreeVector(x,y,z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MedicalBeam::GeneratePrimaries(G4Event* anEvent)
{
  if ( fparticle == NULL ) return;

  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(fsourcePosition, 0.*ns);

  // momentum
  G4double mass = fparticle-> GetPDGMass();
  G4double p = std::sqrt(sqr(mass+fkineticE)-sqr(mass));
  G4ThreeVector pmon = p*GenerateBeamDirection();
  G4PrimaryParticle* primary =
    new G4PrimaryParticle(fparticle, pmon.x(), pmon.y(), pmon.z());

  // set primary to vertex
  vertex-> SetPrimary(primary);

  // set vertex to event
  anEvent-> AddPrimaryVertex(vertex);
}
