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
//
// 080602 Fix memory leaks by T. Koi 
// 081120 Add EraseParticipant and InsertParticipant Methods by T. Koi

#ifndef G4QMDSystem_hh
#define G4QMDSystem_hh

#include "G4QMDParticipant.hh"

class G4QMDSystem 
{
   public:

      G4QMDSystem();
      virtual ~G4QMDSystem();

      void SetParticipant( G4QMDParticipant* particle )
      {
        participants.push_back ( particle );
      }

      void SetSystem ( G4QMDSystem* , G4ThreeVector , G4ThreeVector );

      void SubtractSystem ( G4QMDSystem* );

      G4QMDParticipant* EraseParticipant( G4int i )
      {
        G4QMDParticipant* particle =  participants[ i ];
        participants.erase(std::find( participants.cbegin(),
                                      participants.cend(), participants[ i ]));
        return particle;
      }

      void DeleteParticipant( G4int i )
      {
        delete participants[ i ];
        participants.erase(std::find ( participants.cbegin(),
                                       participants.cend(), participants[ i ]));
      }

      void InsertParticipant( G4QMDParticipant* particle , G4int j );

      G4int GetTotalNumberOfParticipant()
      {
        return (G4int)participants.size();
      }

      G4QMDParticipant* GetParticipant( G4int i )
      {
        return participants[i];
      }

      void IncrementCollisionCounter()
      {
        ++numberOfCollision;
      }

      G4int GetNOCollision()
      {
        return numberOfCollision;
      }

      void ShowParticipants();

      void Clear();

   protected:
      std::vector< G4QMDParticipant* > participants;

   private:
      G4int numberOfCollision;
};

#endif
