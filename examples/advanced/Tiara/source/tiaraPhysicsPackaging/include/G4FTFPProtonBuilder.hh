#ifndef G4FTFPProtonBuilder_h
#define G4FTFPProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4NeutronInelasticCrossSection.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4ProtonInelasticCrossSection.hh"

class G4FTFPProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4FTFPProtonBuilder();
    virtual ~G4FTFPProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4TheoFSGenerator * theModel;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel * thePreEquilib;
    G4GeneratorPrecompoundInterface * theCascade;
    G4FTFModel theStringModel;
    G4LundStringFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

