//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsQGSP_BIC_AllHP_h
#define G4HadronPhysicsQGSP_BIC_AllHP_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4BinaryDeuteronBuilder.hh"
#include "G4BinaryTritonBuilder.hh"
#include "G4BinaryHe3Builder.hh"
#include "G4BinaryAlphaBuilder.hh"
#include "G4ProtonPHPBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4NeutronPHPBuilder.hh"

#include "G4DeuteronBuilder.hh"
#include "G4DeuteronPHPBuilder.hh"

#include "G4TritonBuilder.hh"
#include "G4TritonPHPBuilder.hh"

#include "G4He3Builder.hh"
#include "G4He3PHPBuilder.hh"

#include "G4AlphaBuilder.hh"
#include "G4AlphaPHPBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"


class G4HadronPhysicsQGSP_BIC_AllHP : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsQGSP_BIC_AllHP(G4int verbose =1);
    G4HadronPhysicsQGSP_BIC_AllHP(const G4String& name, G4bool quasiElastic=true);
    virtual ~G4HadronPhysicsQGSP_BIC_AllHP();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();

    struct ThreadPrivate {
      G4NeutronBuilder * theNeutronB;
      G4FTFPNeutronBuilder * theFTFPNeutron;
      G4QGSPNeutronBuilder * theQGSPNeutron;
      G4BinaryNeutronBuilder * theBinaryNeutron;
      G4NeutronPHPBuilder * thePHPNeutron;
    
      G4PiKBuilder * thePiKB;
      G4FTFPPiKBuilder * theFTFPPiK;
      G4QGSPPiKBuilder * theQGSPPiK;
      G4BertiniPiKBuilder * theBertiniPiK;
    
      G4ProtonBuilder * theProtonB;
      G4FTFPProtonBuilder * theFTFPProton;
      G4QGSPProtonBuilder * theQGSPProton;
      G4BinaryProtonBuilder * theBinaryProton;
      G4ProtonPHPBuilder * thePHPProton;

      G4DeuteronBuilder * theDeuteronB;
      G4DeuteronPHPBuilder * thePHPDeuteron; 
      G4BinaryDeuteronBuilder * theBinaryDeuteron;

      G4TritonBuilder * theTritonB;
      G4TritonPHPBuilder * thePHPTriton; 
      G4BinaryTritonBuilder * theBinaryTriton;

      G4He3Builder * theHe3B;
      G4He3PHPBuilder * thePHPHe3; 
      G4BinaryHe3Builder * theBinaryHe3;
      
      G4AlphaBuilder * theAlphaB;
      G4AlphaPHPBuilder * thePHPAlpha; 
      G4BinaryAlphaBuilder * theBinaryAlpha;

      G4HyperonFTFPBuilder * theHyperon;

      G4AntiBarionBuilder * theAntiBaryon;
      G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

      G4VCrossSectionDataSet * xsNeutronCaptureXS;
    };
    static G4ThreadLocal ThreadPrivate* tpdata;

    // G4bool QuasiElastic;
};

#endif

