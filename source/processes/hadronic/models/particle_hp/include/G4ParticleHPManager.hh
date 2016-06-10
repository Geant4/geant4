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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPManager_h
#define G4ParticleHPManager_h 1

// Class Description
// Manager of NeutronHP 
// Class Description - End

// 121031 First implementation done by T. Koi (SLAC/PPA)
//
#include <map>
#include "globals.hh"

#include "G4ParticleHPReactionWhiteBoard.hh"

class G4ParticleHPManager 
{
   public:
      static G4ParticleHPManager* GetInstance() {
         if ( instance == NULL) instance = new G4ParticleHPManager();
         return instance;
      };

   private: 
      G4ParticleHPManager();
      G4ParticleHPManager( const G4ParticleHPManager& ){};
      ~G4ParticleHPManager();
      static G4ThreadLocal G4ParticleHPManager* instance;

   public:
      G4ParticleHPReactionWhiteBoard* GetReactionWhiteBoard();
      void OpenReactionWhiteBoard();
      void CloseReactionWhiteBoard(){delete RWB; RWB=NULL;};

      void GetDataStream( G4String , std::istringstream& iss );
      void GetDataStream2( G4String , std::istringstream& iss );
      void SetVerboseLevel( G4int i ); 
      G4int GetVerboseLevel() {return verboseLevel; }; 

      void DumpDataSource();
   private:
      void register_data_file( G4String , G4String );
      std::map<G4String,G4String> mDataEvaluation;
      G4ParticleHPReactionWhiteBoard* RWB;

      G4int verboseLevel;
};
#endif
