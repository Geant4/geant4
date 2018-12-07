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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyPrimaryGeneratorMessenger.hh"

#include "HadrontherapyMatrix.hh"
#include "HadrontherapyDetectorSD.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"


#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"

#include "G4Event.hh"
#include "G4Timer.hh"

#include "G4RunManager.hh"



/////////////////////////////////////////////////////////////////////////////
HadrontherapyPrimaryGeneratorAction::HadrontherapyPrimaryGeneratorAction()
{
    PrimaryGeneratorMessenger = new HadrontherapyPrimaryGeneratorMessenger(this);
    SetDefaultPrimaryParticle();
    particleGun = new G4GeneralParticleSource();
    
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPrimaryGeneratorAction::~HadrontherapyPrimaryGeneratorAction()
{
    delete PrimaryGeneratorMessenger;
    delete  particleGun;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPrimaryGeneratorAction::SetDefaultPrimaryParticle()
{
    
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    if(NewSource==true)
    {
        std::ifstream in(calculatedPhaseSpaceFileIN);
        G4double e, xpos, ypos, zpos,dirx,diry,dirz;
        G4int PDG;
        G4ThreeVector pos,dir;
        
        if(in.eof())
        {
            G4Exception("HadrontherapyPrimaryGeneratorAction", "NoParticles", FatalException, "No more particles in the file");
        }
        
        while(!in.eof())
        {
            
            in >> e >> xpos >> ypos >>zpos >>dirx>>diry>>dirz >> PDG;
            dir= G4ThreeVector(dirx,diry,dirz);
            particleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(e);
            
            particleGun->GetCurrentSource()->GetParticlePosition().setX(xpos);
            particleGun->GetCurrentSource()->GetParticlePosition().setY(ypos);
            particleGun->GetCurrentSource()->GetParticlePosition().setZ(zpos);
            particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(dir);
            
            G4ParticleDefinition* particleDef = nullptr;
            if (PDG > 1000000000)
            {
                int a=(PDG-1000000000)-(((PDG-1000000000)/10)*10);
                if(a>0)
                {
                    PDG=PDG-a;
                    particleDef = G4IonTable::GetIonTable()->GetIon(PDG);
                    G4String Nome = particleDef->GetParticleName();
                }
                
                else
                {
                    particleDef = G4IonTable::GetIonTable()->GetIon(PDG);
                    G4String Nome = particleDef->GetParticleName();
                }
            }
            
            else
            {
                particleDef = G4ParticleTable::GetParticleTable()->FindParticle(PDG);
            }
            
            particleGun->GetCurrentSource()->SetParticleDefinition(particleDef);
            particleGun->GeneratePrimaryVertex(anEvent);
            
        }
        
        in.close();
        
    }
    else
    {
        particleGun->GeneratePrimaryVertex(anEvent);
    }
    
}


