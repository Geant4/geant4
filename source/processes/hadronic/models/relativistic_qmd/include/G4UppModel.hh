
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

  G4UppModel() 
    : propagationTime(80),myHandler(allTracks) {}

  void initialize(const G4KineticTrackVector& aState);
  void initialize(const G4KineticTrackVector& aProjectile, 
		  const G4KineticTrackVector& aTarget);
  void initialize(G4Fancy3DNucleus& aProjectile, 
		  G4Fancy3DNucleus& aTarget);

  void setModelOptions(const G4double newPropTime) 
    { propagationTime=newPropTime; }

  void addAnalyzer(const G4VUppAnalyzer& anAnalyzer, 
		   const G4double analyzeTime);
  void addAnalyzer(const G4VUppAnalyzer& anAnalyzer, 
		   const G4double beginTime,
		   const G4double endTime, 
		   const G4int nSteps);

  G4int propagate(const G4VUppFieldtransport& aTransport);

  G4KineticTrackVector* getFinalState();

private:

  G4UppTrackVector allTracks;
  G4double propagationTime;
  G4UppActionHandler myHandler;

};


#endif // G4UPPMODEL_H





