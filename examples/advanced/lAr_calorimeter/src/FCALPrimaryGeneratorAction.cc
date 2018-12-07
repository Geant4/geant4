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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <fstream>
#include <cstdlib>

#include "FCALPrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4AutoLock.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


// Migration to MT: there is a single input file that is read by all threads.
// The idea is that the events are read by a single thread and processed
// by all threads. Threads ask for the next ID to be processed. When
// events are all processed we start over from the beginning of the file
namespace {
    G4bool isFileRead = false;
    G4Mutex mFileRead = G4MUTEX_INITIALIZER;
    //Primary kinematics
    G4DataVector fX;
    G4DataVector fY;
    G4DataVector fZ;
    G4DataVector fCosX;
    G4DataVector fCosY;
    G4DataVector fCosZ;
    size_t nextEventId = 0;
    G4Mutex mNextEventId = G4MUTEX_INITIALIZER;

    size_t GetNextId() {
        G4AutoLock l(&mNextEventId);
        if ( nextEventId >= fX.size() ) //file data are over,  restart file
            {
                G4Exception("FCALPrimaryGeneratorAction::GeneratePrimaries","lAr002",
                            JustWarning,"Data file with kinematics is over, restart it");
                nextEventId=0;
            }
        return nextEventId++;
    }
    
    void ReadKinematicFromFile(G4double energy) {
        //Only one thread shoud read input file
        G4AutoLock l(&mFileRead);
        if ( isFileRead ) return;
        // Read Kinematics from file
        G4String file_name = "data-tracks/tracks-80GeV.dat";
        if (energy < 30*GeV)
            file_name = "data-tracks/tracks-20GeV.dat";
        else if (energy < 50*GeV)
            file_name = "data-tracks/tracks-40GeV.dat";
        else if (energy < 70*GeV)
            file_name = "data-tracks/tracks-60GeV.dat";
        else if (energy < 90*GeV)
            file_name = "data-tracks/tracks-80GeV.dat";
        else if (energy < 150*GeV)
            file_name = "data-tracks/tracks-120GeV.dat";
        else
            file_name = "data-tracks/tracks-200GeV.dat";
        std::ifstream Traks_file(file_name);
        if(!Traks_file)
        {
            G4ExceptionDescription ed;
            ed << "Failed to open file " << file_name << G4endl;
            G4Exception("FCALPrimaryGeneratorAction::FCALPrimaryGeneratorAction()",
                        "lAr001",FatalException,ed);
        }
        G4double xx=0,yy=0,zz=0,c1=0,c2=0,c3=0;
        G4int iev = 0;
        while(!(Traks_file.eof())) {
            Traks_file >> iev >> xx >> yy >> zz >> c1 >> c2 >> c3;
            fX.push_back(xx*cm);
            fY.push_back(yy*cm);
            fZ.push_back(zz*cm);
            fCosX.push_back(c1);
            fCosY.push_back(c2);
            fCosZ.push_back(c3);
        }
        G4cout << "Read " << fX.size() << " events from file " << file_name << G4endl;
        isFileRead= true;
        Traks_file.close();
        return;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALPrimaryGeneratorAction::FCALPrimaryGeneratorAction() : 
  fVerbosity(0)
{
  particleGun  = new G4ParticleGun();

  // default Particle
  G4String particleName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();  
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);

  // default Energy
  particleGun->SetParticleEnergy(20*GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALPrimaryGeneratorAction::~FCALPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  ReadKinematicFromFile(particleGun->GetParticleEnergy());
    
    //Get next event to be processed
    size_t nEvent = GetNextId();
    particleGun->SetParticlePosition(G4ThreeVector(fX[nEvent],fY[nEvent],fZ[nEvent]));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(-1.0*fCosX[nEvent],
							  fCosY[nEvent],
							  -1.0*fCosZ[nEvent]));

  particleGun->GeneratePrimaryVertex(anEvent);

  if (fVerbosity)
    {
        G4cout<< " Event  "<<anEvent->GetEventID()<< " Generated Vertex : "
            <<anEvent->GetEventID() <<" (x,y,z)=(" << fX[nEvent] << ","
            <<fY[nEvent] << "," << fZ[nEvent]<< ") (cosX,cosY,cosZ)=("
            << -1.*fCosX[nEvent] << "," << fCosY[nEvent]
            <<"," << -1.*fCosZ[nEvent] << ")"<<G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


