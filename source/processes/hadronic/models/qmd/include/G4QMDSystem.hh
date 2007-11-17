// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:    G4QMDSystem.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 29 March 2007
// -----------------------------------------------------------------------------

#ifndef G4QMDSystem_hh
#define G4QMDSystem_hh

#include "G4QMDParticipant.hh"

class G4QMDSystem 
{
   public:
      G4QMDSystem();
      ~G4QMDSystem();

      void SetParticipant( G4QMDParticipant* particle ) { participants.push_back ( particle ); };
      void SetSystem ( G4QMDSystem* , G4ThreeVector , G4ThreeVector );

      void SubtractSystem ( G4QMDSystem* );

      void DeleteParticipant( G4int i ) { participants.erase( std::find ( participants.begin() , participants.end() , participants[ i ] ) ); };

      G4int GetTotalNumberOfParticipant() { return participants.size(); };

      G4QMDParticipant* GetParticipant( G4int i ) { return participants[i]; };

      void IncrementCollisionCounter() { numberOfCollision++; };
      G4int GetNOCollision() { return numberOfCollision; };

      void ShowParticipants();

      void Clear();

   protected:
      std::vector< G4QMDParticipant* > participants;

   private:
      G4int numberOfCollision;
};

#endif
