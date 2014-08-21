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
#ifndef G4NeutronHPManager_h
#define G4NeutronHPManager_h 1

// Class Description
// Manager of NeutronHP 
// Class Description - End

// 121031 First implementation done by T. Koi (SLAC/PPA)
//
#include <map>
#include "globals.hh"

#include "G4NeutronHPReactionWhiteBoard.hh"

class G4NeutronHPMessenger;

class G4NeutronHPManager 
{
   public:
      static G4NeutronHPManager* GetInstance() {
         if ( instance == NULL) instance = new G4NeutronHPManager();
         return instance;
      };

   private: 
      G4NeutronHPManager();
      G4NeutronHPManager( const G4NeutronHPManager& ){};
      ~G4NeutronHPManager();
      static G4ThreadLocal G4NeutronHPManager* instance;

   public:
      G4NeutronHPReactionWhiteBoard* GetReactionWhiteBoard();
      void OpenReactionWhiteBoard();
      void CloseReactionWhiteBoard(){delete RWB; RWB=NULL;};

      void GetDataStream( G4String , std::istringstream& iss );
      void GetDataStream2( G4String , std::istringstream& iss );
      void SetVerboseLevel( G4int i ); 
      G4int GetVerboseLevel() {return verboseLevel; }; 

      void DumpDataSource();

      G4bool GetUseOnlyPhotoEvaporation() { return USE_ONLY_PHOTONEVAPORATION; };
      void SetUseOnlyPhotoEvaporation( G4bool val ) { USE_ONLY_PHOTONEVAPORATION = val; };
      G4bool GetSkipMissingIsotopes() { return SKIP_MISSING_ISOTOPES; };
      G4bool GetNeglectDoppler() { return NEGLECT_DOPPLER; };
      G4bool GetDoNotAdjustFinalState() { return DO_NOT_ADJUST_FINAL_STATE; };
      G4bool GetProduceFissionFragments() { return PRODUCE_FISSION_FRAGMENTS; };

      void SetSkipMissingIsotopes( G4bool val ) { SKIP_MISSING_ISOTOPES = val; };
      void SetNeglectDoppler( G4bool val ) { NEGLECT_DOPPLER = val; };
      void SetDoNotAdjustFinalState( G4bool val ) { DO_NOT_ADJUST_FINAL_STATE = val; };
      void SetProduceFissionFragments( G4bool val ) { PRODUCE_FISSION_FRAGMENTS = val; };

   private:
      void register_data_file( G4String , G4String );
      std::map<G4String,G4String> mDataEvaluation;
      G4NeutronHPReactionWhiteBoard* RWB;

      G4int verboseLevel;

      G4NeutronHPMessenger* messenger;
      G4bool USE_ONLY_PHOTONEVAPORATION;
      G4bool SKIP_MISSING_ISOTOPES;
      G4bool NEGLECT_DOPPLER;
      G4bool DO_NOT_ADJUST_FINAL_STATE;
      G4bool PRODUCE_FISSION_FRAGMENTS;
 
};
#endif
