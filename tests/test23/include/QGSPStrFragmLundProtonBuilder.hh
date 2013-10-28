//
//---------------------------------------------------------------------------
//
// ClassName:   QGSPStrFragmLundProtonBuilder
//
// Author: Julia Yarba, FNAL/CD (2013)
//
//----------------------------------------------------------------------------
//
#ifndef QGSPStrFragmLundProtonBuilder_h
#define QGSPStrFragmLundProtonBuilder_h 

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
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"
// #include "G4ProjectileDiffractiveChannel.hh"

class QGSPStrFragmLundProtonBuilder : public G4VProtonBuilder
{
  public: 
  
    // ctor & dtor
    QGSPStrFragmLundProtonBuilder( G4bool quasiElastic=true ); 
    virtual ~QGSPStrFragmLundProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess * aP);
    virtual void Build(G4ProtonInelasticProcess * aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:

    G4TheoFSGenerator*               theModel;
    G4PreCompoundModel*              thePreEquilib;
    G4GeneratorPrecompoundInterface* theCascade;
    G4QGSModel< G4QGSParticipants >* theStringModel;
    G4ExcitedStringDecay*            theStringDecay;
    G4QuasiElasticChannel*           theQuasiElastic;
    G4LundStringFragmentation*       theStrFragm;
    G4ExcitationHandler*             theHandler;
    G4double                         theMin;

};

#endif

