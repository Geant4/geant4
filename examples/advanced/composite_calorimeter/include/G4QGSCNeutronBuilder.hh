#ifndef G4QGSCNeutronBuilder_h
#define G4QGSCNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4NeutronInelasticCrossSection.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

class G4QGSCNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4QGSCNeutronBuilder();
    virtual ~G4QGSCNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4QGSModel< G4QGSParticipants > theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;

    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

