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
// $Id: F01PrimaryGeneratorAction.cc,v 1.8 2006-06-29 17:16:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F01PrimaryGeneratorAction.hh"

#include "F01DetectorConstruction.hh"
#include "F01PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
 G4String F01PrimaryGeneratorAction::thePrimaryParticleName = "e-" ; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01PrimaryGeneratorAction::F01PrimaryGeneratorAction(
                                            F01DetectorConstruction* F01DC)
  : F01Detector(F01DC), rndmFlag("off"),
    xvertex(0.), yvertex(0.), zvertex(0.),
    vertexdefined(false)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  // create a messenger for this class
  gunMessenger = new F01PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  
  thePrimaryParticleName = particle->GetParticleName() ;

  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(0.5*GeV);

  zvertex = F01Detector->GetAbsorberZpos() -0.5*(F01Detector->GetAbsorberThickness());
  particleGun->SetParticlePosition(G4ThreeVector(xvertex,yvertex,zvertex));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01PrimaryGeneratorAction::~F01PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  // 
  thePrimaryParticleName = particleGun->GetParticleDefinition()->
                                                GetParticleName() ;
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
    z0 = F01Detector->GetAbsorberZpos()-0.5*(F01Detector->GetAbsorberThickness());
  }
  G4double r0,phi0 ;

  if (rndmFlag == "on")
  {
      r0 = (F01Detector->GetAbsorberRadius())*std::sqrt(G4UniformRand());
      phi0 = twopi*G4UniformRand();
      x0 = r0*std::cos(phi0);
      y0 = r0*std::sin(phi0);
  } 

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
}

///////////////////////////////////////////////////////////////////////
//
//

G4String F01PrimaryGeneratorAction::GetPrimaryName()
{
   return thePrimaryParticleName ;
}

void F01PrimaryGeneratorAction::Setzvertex(G4double z)
{
  vertexdefined = true ;
  zvertex = z ;
  G4cout << " Z coordinate of the primary vertex = " << zvertex/mm <<
            " mm." << G4endl;
}
void F01PrimaryGeneratorAction::Setxvertex(G4double x)
{
  vertexdefined = true ;
  xvertex = x ;
  G4cout << " X coordinate of the primary vertex = " << xvertex/mm <<
            " mm." << G4endl;
}

void F01PrimaryGeneratorAction::Setyvertex(G4double y)
{
  vertexdefined = true ;
  yvertex = y ;
  G4cout << " Y coordinate of the primary vertex = " << yvertex/mm <<
            " mm." << G4endl;
}











