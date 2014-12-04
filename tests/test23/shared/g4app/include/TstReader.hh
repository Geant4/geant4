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
#ifndef TstReader_H
#define TstReader_H 1

#include "G4String.hh"
#include "G4ThreeVector.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>


class TstReader
{

   public:
   
      // ctor & dtor
      TstReader();
      virtual ~TstReader();
      
      void Help();
      void OpenAppConfig( std::string conf );
      void CloseAppConfig();
      
      G4int           GetNEvents()        const { return fNEvents; }
      G4String        GetBeamParticle()   const { return fBeamPart; }
      G4double        GetBeamEnergy()     const { return fBeamEnergy; }
      G4double        GetBeamKineticEnergy() const { return fBeamEKin; } 
      G4double        GetBeamMomentum()   const { return fBeamMomentum; }
      G4ThreeVector   GetDirection()      const { return fDirection; }
      G4ThreeVector   GetPosition()       const { return fPosition; }
      G4double        GetTime()           const { return fTime; }
      G4String        GetTargetMaterial() const { return fTargetMaterial; }
      G4ThreeVector   GetTargetSize()     const { return fTargetSize; }
      G4String        GetTargetShape()    const { return fTargetShape; }
      G4double        GetStep()           const { return fStep; }
      long            GetRndmSeed()       const { return fRndmSeed; }
      G4int           GetJobID()          const { return fJobID; }
      G4int           GetClusterID()      const { return fClusterID; }
      G4int           GetVerbosity()      const { return fVerbose; } 
      G4String        GetExpDataSet()     const { return fExpDataSet; } 
      G4String        GetPhysics()        const { return fPhysics; }  
      G4bool          ForceResDecay()     const { return fForceResDecay; }   
            
      G4bool IsDone() { return fEndConfig; }
      void ProcessConfig();
      void SyncKinematics() const;
            
      virtual G4bool IsDiscreteProcess() const = 0;
      virtual G4bool IsAtRestProcess() const = 0;
      virtual G4bool IsPhysicsList() const = 0; 
      
   protected:

      virtual void ProcessLine( G4String line );
      void SetExpDataSet( G4String dset ) { fExpDataSet=dset; return; }
      void SetPhysics( G4String physname ) { fPhysics=physname; return; }
      
      std::ifstream* fInStream;
      G4bool         fEndConfig;
      
      G4int          fNEvents;
      G4String       fBeamPart;
      mutable G4double       fBeamEnergy;   // MeV
      mutable G4double       fBeamEKin;     // MeV
      G4double       fBeamMomentum; // MeV/c
      G4ThreeVector  fDirection;
      G4ThreeVector  fPosition;     // mm
      G4double       fTime;         // ns
      G4String       fTargetMaterial;
      G4ThreeVector  fTargetSize;   // in mm
      G4String       fTargetShape;
      G4String       fPhysics;
      G4double       fStep;         // mm ...not really used though... 
      long           fRndmSeed;
      G4int          fJobID;
      G4int          fClusterID;
      G4int          fVerbose;
      G4String       fExpDataSet;
      G4bool         fForceResDecay;
      
};



#endif
