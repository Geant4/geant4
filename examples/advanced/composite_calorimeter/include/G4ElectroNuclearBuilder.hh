#ifndef G4ElectroNuclearBuilder_h
#define G4ElectroNuclearBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4GammaNuclearReaction.hh"
#include "G4ElectroNuclearReaction.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

class G4ElectroNuclearBuilder 
{
  public: 
    G4ElectroNuclearBuilder();
    virtual ~G4ElectroNuclearBuilder();

  public: 
    virtual void Build();

  protected:
    G4PhotoNuclearProcess thePhotoNuclearProcess;
    G4ElectronNuclearProcess theElectronNuclearProcess;
    G4PositronNuclearProcess thePositronNuclearProcess;
    G4ElectroNuclearReaction * theElectroReaction;
    G4GammaNuclearReaction * theGammaReaction;  
    
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4QGSModel< G4GammaParticipants > theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
};

// 2002 by J.P. Wellisch

#endif





