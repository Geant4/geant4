#ifndef G4FTFPPiKBuilder_h
#define G4FTFPPiKBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4PiNuclearCrossSection.hh"

class G4FTFPPiKBuilder : public G4VPiKBuilder
{
  public: 
    G4FTFPPiKBuilder();
    virtual ~G4FTFPPiKBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4PionPlusInelasticProcess & aP);
    virtual void Build(G4PionMinusInelasticProcess & aP);
    virtual void Build(G4KaonPlusInelasticProcess & aP);
    virtual void Build(G4KaonMinusInelasticProcess & aP);
    virtual void Build(G4KaonZeroLInelasticProcess & aP);
    virtual void Build(G4KaonZeroSInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4PiNuclearCrossSection thePiData;
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

