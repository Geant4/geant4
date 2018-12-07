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
/// \file field/field04/src/F04PrimaryGeneratorAction.cc
/// \brief Implementation of the F04PrimaryGeneratorAction class
//

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4GeometryManager.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "F04PrimaryGeneratorAction.hh"

#include "F04DetectorConstruction.hh"
#include "F04PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PrimaryGeneratorAction::
       F04PrimaryGeneratorAction(F04DetectorConstruction* detectorConstruction)
  : fDetector(detectorConstruction), fRndmFlag("off"), fFirst(false),
    fXvertex(0.), fYvertex(0.), fZvertex(0.),
    fVertexdefined(false)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  fGunMessenger = new F04PrimaryGeneratorMessenger(this);

  G4String particleName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  fParticleGun->SetParticleDefinition(particleTable->
                                        FindParticle(particleName="proton"));
  fParticleGun->SetParticleEnergy(500.*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  fZvertex = -0.5*(fDetector->GetTargetThickness());
  fParticleGun->SetParticlePosition(G4ThreeVector(fXvertex,fYvertex,fZvertex));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PrimaryGeneratorAction::~F04PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if (!fFirst) {

     fFirst = true;
     G4ThreeVector direction(0.0,0.0,1.0);

     G4Navigator* theNavigator =
                    G4TransportationManager::GetTransportationManager()->
                                                 GetNavigatorForTracking();
     if ( theNavigator->GetWorldVolume() )
     {
       G4Navigator* aNavigator = new G4Navigator();
       aNavigator->SetWorldVolume(theNavigator->GetWorldVolume());

       G4ThreeVector center(0.,0.,0.);
       aNavigator->LocateGlobalPointAndSetup(center,0,false);

       G4TouchableHistoryHandle touchable = aNavigator->
                                          CreateTouchableHistoryHandle();

       // set Global2local transform
       fGlobal2local = touchable->GetHistory()->GetTopTransform();

       direction = fGlobal2local.Inverse().TransformAxis(direction);
       delete aNavigator;
     }

     fParticleGun->SetParticleMomentumDirection(direction);
  }

  G4double x0,y0,z0 ;

  if(fVertexdefined)
  {
    x0 = fXvertex ;
    y0 = fYvertex ;
    z0 = fZvertex ;
  }
  else
  {
    x0 = 0. ;
    y0 = 0. ;
    z0 = -0.5*(fDetector->GetTargetThickness());
  }

  G4double r0, phi0;

  if (fRndmFlag == "on")
  {
      r0 = (fDetector->GetTargetRadius())*std::sqrt(G4UniformRand());
      phi0 = twopi*G4UniformRand();
      x0 = r0*std::cos(phi0);
      y0 = r0*std::sin(phi0);
  }

  G4ThreeVector localPosition(x0,y0,z0);
  G4ThreeVector globalPosition =
                        fGlobal2local.Inverse().TransformPoint(localPosition);

  fParticleGun->SetParticlePosition(globalPosition);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PrimaryGeneratorAction::SetXvertex(G4double x)
{
  fVertexdefined = true;
  fXvertex = x;
  G4cout << " X coordinate of the primary vertex = " << fXvertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PrimaryGeneratorAction::SetYvertex(G4double y)
{
  fVertexdefined = true;
  fYvertex = y;
  G4cout << " Y coordinate of the primary vertex = " << fYvertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PrimaryGeneratorAction::SetZvertex(G4double z)
{
  fVertexdefined = true;
  fZvertex = z;
  G4cout << " Z coordinate of the primary vertex = " << fZvertex/mm <<
            " mm." << G4endl;
}
