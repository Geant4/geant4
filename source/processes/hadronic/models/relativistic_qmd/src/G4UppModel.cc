
#include "G4UppModel.hh"
#include "G4UppActionHandler.hh"
#include "G4VUppAction.hh"
#include "G4UppActionAnalyze.hh"
#include "G4UppTrackChange.hh"
#include "G4VUppFieldtransport.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4UppGlobal.hh"


void G4UppModel::initialize (G4Fancy3DNucleus& aProjectile, 
			     G4Fancy3DNucleus& aTarget)
{
  allTracks.clear();
  allTracks.add(aProjectile.GetNucleons(),1);
  allTracks.add(aTarget.GetNucleons(),2);
  allTracks.setGlobalTime(0.0);
  // allTracks.dump();
}


void G4UppModel::initialize (const G4KineticTrackVector& tracks)
{
  allTracks.clear();
  for (G4int i=0; i<tracks.length(); i++) {
    G4UppTrack* aTrackPtr = new G4UppTrack( *(tracks[i]) );
    allTracks.push_back(aTrackPtr);
  }
  allTracks.setGlobalTime(0.0);
  // allTracks.dump();
}


G4int G4UppModel::propagate(const G4VUppFieldtransport& aTransport)
{
  const G4VUppAction* actualActionPtr;
  G4UppTrackChange* aTrackChangePtr;

  // Fill Collision Array
  myHandler.updateActions();
  allTracks.resetChangedFlags();
  // myHandler.dump(allTracks);

  do {
    // get the first action to be done
    actualActionPtr = myHandler.getFirstAction();
    cout << "(debug) get action:";
    actualActionPtr->dump();

    G4double dTime = actualActionPtr->getActionTime() - allTracks.getGlobalTime();
    if (dTime < 0) {
      cout << "(warning) action already gone!" << endl;
      cout << "(warning)  GlobalTime: " << allTracks.getGlobalTime()/fermi << endl;
      cout << "(warning)  removing action..." << endl;
    } else {
      // propagate the whole system to the time of the action
      aTransport.propagate(allTracks,dTime);
    
      // now perform the action
      if (aTrackChangePtr = actualActionPtr->perform(allTracks)) {
	G4UppTrackChange undoBackupChange = allTracks.update(*aTrackChangePtr);
	if (allTracks.isPauliBlocked(aTrackChangePtr->newParticles))
	  allTracks.update(undoBackupChange);
      }
    }
    myHandler.deleteFirstAction();
    myHandler.cleanUp();
    myHandler.updateActions();
    allTracks.resetChangedFlags();
    
  } while ((allTracks.getGlobalTime() < propagationTime) && !myHandler.empty());

  return 0;
}


void G4UppModel::addAnalyzer (const G4VUppAnalyzer& anAnalyzer, 
			      const G4double analyzeTime)
{
  G4VUppAction* actionPtr = new G4UppActionAnalyze(analyzeTime, anAnalyzer);
  myHandler.addAction(actionPtr);
}



