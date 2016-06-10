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
// 081120 Add InsertParticipant Method by T. Koi

#include "G4QMDSystem.hh"
#include <iomanip>

G4QMDSystem::G4QMDSystem()
{
   participants.clear();
   numberOfCollision = 0;
}



G4QMDSystem::~G4QMDSystem()
{
   this->Clear();
}


// Insert nucleus to current system;
void G4QMDSystem::SetSystem ( G4QMDSystem* nucleus , G4ThreeVector dp , G4ThreeVector dr )
{
   std::vector< G4QMDParticipant* >::iterator it; 
   for ( it = nucleus->participants.begin() ; it != nucleus->participants.end() ; it++ ) 
   {
      G4ThreeVector r = (*it)->GetPosition() + dr;
      (*it)->SetPosition ( r );
      G4ThreeVector p = (*it)->GetMomentum() + dp;
      (*it)->SetMomentum ( p );
      this->SetParticipant( *it );
   }
}

void G4QMDSystem::SubtractSystem ( G4QMDSystem* nucleus )
{
   
   for ( G4int i = 0 ; i < nucleus->GetTotalNumberOfParticipant() ; i++ )
   {
      participants.erase ( std::find ( participants.begin() , participants.end() , nucleus->GetParticipant( i ) ) );
   }
}

void G4QMDSystem::Clear ()
{
   for ( G4int i = 0 ; i < this->GetTotalNumberOfParticipant() ; i++ )
   {
      delete participants[i];
   }
   participants.clear();
}



void G4QMDSystem::ShowParticipants()
{
   //store orginal precision
   std::ios::fmtflags oldform = G4cout.flags();

   G4ThreeVector p_sum( 0.0 );
   std::vector< G4QMDParticipant* >::iterator it; 
   G4cout << "Momentum and Position of each participant " << G4endl; 
   G4int i = 0; 
   for ( it = participants.begin() ; it != participants.end() ; it++ ) 
   { 
      G4cout << i
             << " " 
             << (*it)->GetDefinition()->GetParticleName() 
             << " " 
             << std::setprecision( 8 )
             << (*it)->GetMomentum()
             << " " 
             << (*it)->GetPosition() 
             << G4endl;
      p_sum += (*it)->GetMomentum();
      i++;
   }
   G4cout << "Sum upped Momentum and its mag " << p_sum << " " << p_sum.mag() << G4endl;

   //restore orginal precision
   G4cout.flags( oldform );
}



void G4QMDSystem::InsertParticipant ( G4QMDParticipant* particle , G4int n )
{

   if ( (size_t) n > participants.size()+1 )
      G4cout << "G4QMDSystem::InsertParticipant size error" << G4endl;

   std::vector< G4QMDParticipant* >::iterator it; 
   it = participants.begin();

   for ( G4int i = 0; i < n ; i++ )
      it++;

   participants.insert( it, particle );
}
