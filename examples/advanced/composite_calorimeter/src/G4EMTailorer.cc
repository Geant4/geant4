#include "G4EMTailorer.hh"
#include "G4EMBuilder.hh"

void G4EMTailorer::SetNewValue(G4UIcommand* aComm, G4String aS) 
{ 
  if(aComm==theSynch) 
  {
    //cout <<" ############### synched ############## "<<endl;
    theB->Synch(aS);
  }
  if(aComm==theGN) 
  {
    theB->GammaNuclear(aS);
    //cout <<" ############### GN'ed ############## "<<endl;
  }
}

G4EMTailorer::G4EMTailorer(G4EMBuilder * ab)
{
  theB = ab;
  aDir1 = new G4UIdirectory("/physics_engine/");
  aDir1->SetGuidance("commands related to the physics simulation engine.");
  
  // general stuff.
  aDir2 = new G4UIdirectory("/physics_engine/tailor/");
  aDir2->SetGuidance("tailoring the processes");
  
  // command for synchrotron radiation.
  theSynch = new G4UIcmdWithAString("/physics_engine/tailor/SyncRadiation",this);
  theSynch->SetGuidance("Switching on/off synchrotron radiation.");
  theSynch->SetParameterName("status","off");
  theSynch->SetCandidates("on off");
  theSynch->SetDefaultValue("off");
  theSynch->AvailableForStates(G4State_PreInit);    
  
  // command for gamma nuclear physics.  
  theGN = new G4UIcmdWithAString("/physics_engine/tailor/GammaNuclear",this);
  theGN->SetGuidance("Switching on gamma nuclear physics.");
  theGN->SetParameterName("status","off");
  theGN->SetCandidates("on off");
  theGN->SetDefaultValue("off");
  theGN->AvailableForStates(G4State_PreInit);    
}
    
