#ifndef G4FTFCNeutronBuilder_h
#define G4FTFCNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4NeutronInelasticCrossSection.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

class G4FTFCNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4FTFCNeutronBuilder();
    virtual ~G4FTFCNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4FTFModel theStringModel;
    G4LundStringFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;

    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

