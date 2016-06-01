#include "G4DiffractiveHHScatterer.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4ExcitedString.hh"
#include "G4LundStringFragmentation.hh"
#include "G4KineticTrack.hh"
#include "G4DiffractiveSplitableHadron.hh"

G4DiffractiveHHScatterer::G4DiffractiveHHScatterer()
:
  theExcitation(new G4DiffractiveExcitation()),
  theStringFragmentation(new G4LundStringFragmentation())
{}

G4KineticTrackVector * G4DiffractiveHHScatterer::
Scatter(const G4KineticTrack & aTrack, const G4KineticTrack & bTrack)
{
  G4KineticTrackVector * result = new G4KineticTrackVector();

  G4DiffractiveSplitableHadron aHadron(& aTrack);
  G4DiffractiveSplitableHadron bHadron(& bTrack);
  if ( ! theExcitation->ExciteParticipants(& aHadron, & bHadron)) 
  {
	return NULL;
  }
  
  G4ExcitedString * string;
  G4KineticTrackVector * fragments=NULL;

  string = theExcitation->String(& aHadron, true);  // Projectile
  fragments=theStringFragmentation->FragmentString(*string);
  for (G4int aFragment=0; aFragment < fragments->entries(); aFragment++)
  {
  	result->append(fragments->at(aFragment));
  }
  
  string = theExcitation->String(& bHadron, false); // Target
  fragments=theStringFragmentation->FragmentString(*string);
  for (G4int bFragment=0; bFragment < fragments->entries(); bFragment++)
  {
  	result->append(fragments->at(bFragment));
  }
  
  
  return result;
}
