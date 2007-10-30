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
//

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4GeometryManager.hh"

#include "Randomize.hh"

#include "F04PrimaryGeneratorAction.hh"

#include "F04DetectorConstruction.hh"
#include "F04PrimaryGeneratorMessenger.hh"

G4bool F04PrimaryGeneratorAction::first = false;

F04PrimaryGeneratorAction::F04PrimaryGeneratorAction(F04DetectorConstruction* DC)
  : Detector(DC), rndmFlag("off"),
    xvertex(0.), yvertex(0.), zvertex(0.),
    vertexdefined(false)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  gunMessenger = new F04PrimaryGeneratorMessenger(this);

  G4String particleName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  particleGun->SetParticleDefinition(particleTable->
                                        FindParticle(particleName="proton"));
  particleGun->SetParticleEnergy(500.*MeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  zvertex = -0.5*(Detector->GetTargetThickness());
  particleGun->SetParticlePosition(G4ThreeVector(xvertex,yvertex,zvertex));

}

F04PrimaryGeneratorAction::~F04PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void F04PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if (!first) {

     first = true;

     G4Navigator* theNavigator =
                    G4TransportationManager::GetTransportationManager()->
                                                 GetNavigatorForTracking();
     G4Navigator* aNavigator = new G4Navigator();
     if ( theNavigator->GetWorldVolume() )
               aNavigator->SetWorldVolume(theNavigator->GetWorldVolume());

     G4GeometryManager* geomManager = G4GeometryManager::GetInstance();

     if (!geomManager->IsGeometryClosed()) {
        geomManager->OpenGeometry();
        geomManager->CloseGeometry(true);
     }

     G4ThreeVector center(0.,0.,0.);
     aNavigator->LocateGlobalPointAndSetup(center,0,false);

     G4TouchableHistoryHandle fTouchable = aNavigator->
                                         CreateTouchableHistoryHandle();

    // set global2local transform
    global2local = fTouchable->GetHistory()->GetTopTransform();

    G4ThreeVector direction(0.0,0.0,1.0);
    direction = global2local.Inverse().TransformAxis(direction);

    particleGun->SetParticleMomentumDirection(direction);
  }

  G4double x0,y0,z0 ;

  if(vertexdefined)
  {
    x0 = xvertex ;
    y0 = yvertex ;
    z0 = zvertex ;
  }
  else
  {
    x0 = 0. ;
    y0 = 0. ;
    z0 = -0.5*(Detector->GetTargetThickness());
  }

  G4double r0,phi0 ;

  if (rndmFlag == "on")
  {
      r0 = (Detector->GetTargetRadius())*std::sqrt(G4UniformRand());
      phi0 = twopi*G4UniformRand();
      x0 = r0*std::cos(phi0);
      y0 = r0*std::sin(phi0);
  } 

  G4ThreeVector localPosition(x0,y0,z0);
  G4ThreeVector globalPosition =
                        global2local.Inverse().TransformPoint(localPosition);

  particleGun->SetParticlePosition(globalPosition);
  particleGun->GeneratePrimaryVertex(anEvent);
}

void F04PrimaryGeneratorAction::Setxvertex(G4double x)
{
  vertexdefined = true ;
  xvertex = x ;
  G4cout << " X coordinate of the primary vertex = " << xvertex/mm <<
            " mm." << G4endl;
}

void F04PrimaryGeneratorAction::Setyvertex(G4double y)
{
  vertexdefined = true ;
  yvertex = y ;
  G4cout << " Y coordinate of the primary vertex = " << yvertex/mm <<
            " mm." << G4endl;
}

void F04PrimaryGeneratorAction::Setzvertex(G4double z)
{
  vertexdefined = true ;
  zvertex = z ;
  G4cout << " Z coordinate of the primary vertex = " << zvertex/mm <<
            " mm." << G4endl;
}
