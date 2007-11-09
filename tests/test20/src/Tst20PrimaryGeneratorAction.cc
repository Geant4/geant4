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
// $Id: Tst20PrimaryGeneratorAction.cc,v 1.6 2007-11-09 18:33:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#include "Tst20PrimaryGeneratorAction.hh"

#include "Tst20DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ios.hh"


Tst20PrimaryGeneratorAction::Tst20PrimaryGeneratorAction(Tst20DetectorConstruction* detConstr):detector(detConstr),
											       rndmFlag("off"),
											       xVertex(0.),
											       yVertex(0.),
											       vertexDefined(false)
{
  primaryParticleName = "e-";
  G4int nParticle = 1;
  particleGun  = new G4ParticleGun(nParticle);
  
  //create a messenger for this class
  //  gunMessenger = new Tst20PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(primaryParticleName);
  particleGun->SetParticleDefinition(particle);
  
  // This is nonsense...
  primaryParticleName = particle->GetParticleName();

  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(100.*eV);

  zVertex = -0.5 * detector->GetAbsorberThickness();
  particleGun->SetParticlePosition(G4ThreeVector(xVertex,yVertex,zVertex));

}

Tst20PrimaryGeneratorAction::~Tst20PrimaryGeneratorAction()
{
  delete particleGun;
}

void Tst20PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
 
  primaryParticleName = particleGun->GetParticleDefinition()->GetParticleName() ;
  G4double x0 = 0.;
  G4double y0 = 0.;
  G4double z0 = detector->GetZstartAbs() + 0.5 * detector->GetAbsorberThickness();

  if (vertexDefined)
    {
      x0 = xVertex ;
      y0 = yVertex ;
      z0 = zVertex ;
    }

  G4double r0 = 0.;
  G4double phi0 = 0.;
  
  if (rndmFlag == "on")
    {
      r0 = (detector->GetAbsorberRadius()) * std::sqrt(G4UniformRand());
      phi0 = twopi * G4UniformRand();
      x0 = r0 * std::cos(phi0);
      y0 = r0 * std::sin(phi0);
    } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
}


G4String Tst20PrimaryGeneratorAction::GetPrimaryName()
{
  return primaryParticleName;
}

void Tst20PrimaryGeneratorAction::SetZvertex(G4double z)
{
  vertexDefined = true ;
  zVertex = z ;
  G4cout << " Z coordinate of the primary vertex = " << zVertex/mm << " mm" << G4endl;
}

void Tst20PrimaryGeneratorAction::SetXvertex(G4double x)
{
  vertexDefined = true ;
  xVertex = x ;
  G4cout << " X coordinate of the primary vertex = " << xVertex/mm << " mm" << G4endl;
}

void Tst20PrimaryGeneratorAction::SetYvertex(G4double y)
{
  vertexDefined = true ;
  yVertex = y ;
  G4cout << " Y coordinate of the primary vertex = " << yVertex/mm << " mm" << G4endl;
}
