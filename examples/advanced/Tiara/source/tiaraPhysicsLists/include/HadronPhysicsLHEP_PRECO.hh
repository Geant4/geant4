#ifndef HadronPhysicsLHEP_PRECO_h
#define HadronPhysicsLHEP_PRECO_h 1

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
#include "G4PrecoProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LHEPNeutronBuilder.hh"
#include "G4PrecoNeutronBuilder.hh"

class HadronPhysicsLHEP_PRECO : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsLHEP_PRECO(const G4String& name ="hadron");
    virtual ~HadronPhysicsLHEP_PRECO();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4NeutronBuilder theNeutrons;
    G4LHEPNeutronBuilder theLHEPNeutron;
    G4PrecoNeutronBuilder thePrecoNeutron;
    
    G4PiKBuilder thePiK;
    G4LHEPPiKBuilder theLHEPPiK;
    
    G4ProtonBuilder thePro;
    G4LHEPProtonBuilder theLHEPPro;
    G4PrecoProtonBuilder thePrecoPro;    
    
    G4MiscLHEPBuilder theMiscLHEP;
    G4StoppingHadronBuilder theStoppingHadron;
    G4HadronQEDBuilder theHadronQED;
};

// 2002 by J.P. Wellisch

#endif

