//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

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





