
#ifndef G4UPPMODEL_H
#define G4UPPMODEL_H


#include "G4KineticTrackVector.hh"
#include "G4UppEvent.hh"
#include "G4UppTrackVector.hh"
#include "G4VUppFieldtransport.hh"
#include "G4VUppAnalyzer.hh"
#include "G4UppActionHandler.hh"
#include "G4Fancy3DNucleus.hh"


class G4UppModel
{
public:

  G4UppModel() : PropagationTime(80) {}

  void Initialize (const G4UppEvent& event);
  void Initialize (const G4KineticTrackVector& aState);
  void Initialize (const G4KineticTrackVector& projectile, 
		   const G4KineticTrackVector& target);
  void Initialize (G4Fancy3DNucleus& projectile, 
		   G4Fancy3DNucleus& target);

  void SetModelOptions (const G4double PropTime) 
    { PropagationTime=PropTime; }
  void addAnalyzer (const G4VUppAnalyzer* aPtr, const G4double time);
  void addAnalyzer (const G4VUppAnalyzer* aPtr, const G4double begin,
		    const G4double end, const G4int nSteps);

  G4int Propagate(const G4VUppFieldtransport& aTransport);

  G4KineticTrackVector* GetFinalState();

private:

  G4UppTrackVector allTracks;
  G4double PropagationTime;
  G4UppActionHandler myHandler;

};


#endif // G4UPPMODEL_H





