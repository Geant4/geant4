#ifndef HadronPhysicsLHEP_BIC_h
#define HadronPhysicsLHEP_BIC_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4HadronQEDBuilder.hh"
#include "G4StoppingHadronBuilder.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4LHEPPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LHEPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LHEPNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"

class HadronPhysicsLHEP_BIC : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsLHEP_BIC(const G4String& name ="hadron");
    virtual ~HadronPhysicsLHEP_BIC();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4NeutronBuilder theNeutrons;
    G4LHEPNeutronBuilder theLHEPNeutron;
    G4BinaryNeutronBuilder theBinaryNeutron;
    
    G4PiKBuilder thePiK;
    G4LHEPPiKBuilder theLHEPPiK;
    
    G4ProtonBuilder thePro;
    G4LHEPProtonBuilder theLHEPPro;
    G4BinaryProtonBuilder theBinaryPro;
    
    G4MiscLHEPBuilder theMiscLHEP;
    G4StoppingHadronBuilder theStoppingHadron;
    G4HadronQEDBuilder theHadronQED;
};

// 2002 by J.P. Wellisch

#endif

