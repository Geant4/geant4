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
      G4int           GetVerbosity()      const { return fVerbose; } 
      G4String        GetExpDataSet()     const { return fExpDataSet; } 
      G4String        GetPhysics()        const { return fPhysics; }     
            
      G4bool IsDone() { return fEndConfig; }
      void ProcessConfig();
            
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
      G4double       fBeamEnergy;   // MeV
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
      G4int          fVerbose;
      G4String       fExpDataSet;
      
};



#endif
