//
//   Switch the parameters of this Field Manager between 
//      a 'primary' one (for most particles, energies, ...), and 
//      a 'secondary' one (for particular ones, eg electrons < 10 MeV)
//
//  First implementation: John Apostolakis, 10th June 2003
//

#include "NTSTFieldManager.hh"
#include "G4Track.hh"

NTSTFieldManager::NTSTFieldManager(G4Field         *pCommonField,
				   G4FieldManager  *pfmPrimary, 
				   G4FieldManager  *pfmAlternate,
				   G4double        thresholdKinEnergy) 
  : G4FieldManager( pCommonField, 0, false ),
    pPrimaryFieldMgr(pfmPrimary),
    pAlternativeFieldMgr(pfmAlternate),
    fThresholdKineticEnergy(thresholdKinEnergy),
    fpCurrentFieldMgr(0)
{ 
  G4cout << "Initialising NTST Field Manager with primary FM values." << G4endl;

  // Initialise with the values of the primary Field Manager
  this->CopyValuesAndChordFinder( *pfmPrimary );
  fpCurrentFieldMgr= pfmPrimary;
}
 
void  NTSTFieldManager::ConfigureForTrack( const G4Track *pTrack )
{
   G4double kineticEnergy  = pTrack->GetKineticEnergy() ; 

   if( kineticEnergy < fThresholdKineticEnergy ){
      if( fpCurrentFieldMgr != pAlternativeFieldMgr ){
         this->CopyValuesAndChordFinder( *pAlternativeFieldMgr );
	 fpCurrentFieldMgr= pAlternativeFieldMgr;
	 G4cout << "Choosen alternative Field Manager for NTST FieldManager values." 
		<< G4endl;
      }
   } else { 
      if( fpCurrentFieldMgr != pPrimaryFieldMgr ){
	 this->CopyValuesAndChordFinder( *pPrimaryFieldMgr );
	 fpCurrentFieldMgr= pPrimaryFieldMgr;
	 G4cout << "Choosen primary     Field Manager for NTST FieldManager values." 
		<< G4endl;
      }
   } 
   // G4cout << "Hmin= " 
   //  << G4FieldManager::GetChordFinder()->GetIntegratorDriver()->GetHmin()
   //  << G4endl;
}

const G4FieldManager& NTSTFieldManager::CopyValuesAndChordFinder( const G4FieldManager& newFieldMgr )
  // Copy the accuracy parameters and pointer to ChordFinder; 
  //  do not copy the field pointer -- for now!
{
  G4FieldManager::SetFieldChangesEnergy( newFieldMgr.DoesFieldChangeEnergy() ); 
  
  // this->G4FieldManager::
  SetDeltaIntersection( newFieldMgr.GetDeltaIntersection() ); 

  // G4FieldManager::
  // this->G4FieldManager::
  SetDeltaOneStep( newFieldMgr.GetDeltaOneStep() ); 

  // G4FieldManager nonConst_FieldMgr= const_cast<G4FieldManager&> newFieldMgr; 
  G4FieldManager *nonConst_FieldMgr=  (G4FieldManager*) (&newFieldMgr); 
  G4ChordFinder  *newChordF= nonConst_FieldMgr->GetChordFinder(); 
  G4FieldManager::SetChordFinder( newChordF ); 

  return newFieldMgr;
}


