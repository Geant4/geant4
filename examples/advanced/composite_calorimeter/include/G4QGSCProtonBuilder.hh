#ifndef G4QGSCProtonBuilder_h
#define G4QGSCProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4NeutronInelasticCrossSection.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4ProtonInelasticCrossSection.hh"

class G4QGSCProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4QGSCProtonBuilder();
    virtual ~G4QGSCProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4TheoFSGenerator * theModel;
    G4StringChipsParticleLevelInterface * theCascade;
    G4QGSModel< G4QGSParticipants > theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

