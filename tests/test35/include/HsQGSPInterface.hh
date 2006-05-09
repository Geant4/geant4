#ifndef HsQGSPInterface_h
#define HsQGSPInterface_h

//----------------------------------------------------------------------------
//
//  Package   : Simulation 
//
// Description: Algorithm of G4 HARP for Hadron Production in the target
//
// Author:      V.Ivanchenko 05.03.04
//
// Modifications: 
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VIntraNuclearTransportModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

class HsQGSPInterface : public G4VIntraNuclearTransportModel
{
public:
  HsQGSPInterface();

  ~HsQGSPInterface();

  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                                 G4Nucleus& theNucleus);

  virtual G4ReactionProductVector* Propagate(G4KineticTrackVector*,
                                               G4V3DNucleus*) {return 0;}; 
private:
  
  G4TheoFSGenerator * theModel;
  G4ExcitationHandler theHandler;
  G4PreCompoundModel * thePreEquilib;
  G4GeneratorPrecompoundInterface * theCascade;
  G4QGSModel< G4QGSParticipants > theStringModel;
  G4QGSMFragmentation theFragmentation;
  G4ExcitedStringDecay * theStringDecay;
};
#endif
