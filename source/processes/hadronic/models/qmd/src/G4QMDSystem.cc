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
   G4cout << "Sum upped Momentum and mag " << p_sum << " " << p_sum.mag() << G4endl;
}
