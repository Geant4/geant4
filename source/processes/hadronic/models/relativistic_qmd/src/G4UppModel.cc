
#include "G4UppModel.hh"
#include "G4UppActionHandler.hh"
#include "G4VUppAction.hh"
#include "G4UppActionAnalyze.hh"
#include "G4UppTrackChange.hh"
#include "G4VUppFieldtransport.hh"
#include "G4Fancy3DNucleus.hh"


void G4UppModel::Initialize (G4Fancy3DNucleus& projectile, 
			     G4Fancy3DNucleus& target)
{
  allTracks.clear();
  allTracks.add(projectile.GetNucleons(),1);
  allTracks.add(target.GetNucleons(),2);
  allTracks.setGlobalTime(0.0);
  //  allTracks.dump();
}


void G4UppModel::Initialize (const G4KineticTrackVector& t)
{
  allTracks.clear();
  for (G4int i=0; i<t.length(); i++) {
    G4UppTrack* aPtr = new G4UppTrack(*(t[i]));
    allTracks.push_back(aPtr);
  }
  allTracks.setGlobalTime(0.0);
  allTracks.dump();
}


G4int G4UppModel::Propagate(const G4VUppFieldtransport& aTransportPtr)
{
  G4UppInteraction actualInteraction;
  const G4VUppAction* actualActionPtr;
  G4UppTrackChange aTrackChange;

  // Fill Collision Array
  myHandler.UpdateActions(allTracks);
  allTracks.resetChanged();
  myHandler.dump(allTracks);
  do {
    // get the first action to be done
    actualActionPtr = &(myHandler.getFirstAction());
    cout << "next action:" << endl;
    actualActionPtr->dump();
    cout << "GlobalTime: " << allTracks.getGlobalTime()/fermi << endl;
    if (actualActionPtr->getActionTime() < allTracks.getGlobalTime()) {
      cout << "WARNING action already gone!" << endl;
      cout << "removing action..." << endl;
    } else {
      // propagate the whole system to the time of the action
      G4double dTime = actualActionPtr->getActionTime() - allTracks.getGlobalTime();
      aTransportPtr.Propagate(allTracks,dTime);
    
      // perform the action
      if (actualActionPtr->Perform(allTracks)) {
	actualActionPtr->Perform(allTracks, actualInteraction);
	aTrackChange = actualInteraction.Perform();
	G4UppTrackChange Save = allTracks.Update(aTrackChange);
	if (allTracks.isPauliBlocked(*(actualInteraction.getOutgoingParticles())))
	  allTracks.Update(Save);
      }
    }
    myHandler.deleteFirstAction();
    myHandler.CleanUp();
    myHandler.UpdateActions(allTracks);
    allTracks.resetChanged();
    
  } while ((allTracks.getGlobalTime() < PropagationTime) && !myHandler.empty());

  return 1;
}


void G4UppModel::addAnalyzer (const G4VUppAnalyzer* aPtr, const G4double time)
{
  G4VUppAction* actionPtr = new G4UppActionAnalyze(time, aPtr, &allTracks);
  myHandler.addAction(actionPtr);
}



