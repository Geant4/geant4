#ifndef G4NeutronHPBuilder_h
#define G4NeutronHPBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"

class G4NeutronHPBuilder : public G4VNeutronBuilder
{
  public: 
    G4NeutronHPBuilder();
    virtual ~G4NeutronHPBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);

  private:
    G4NeutronHPElastic * theHPElastic;
    G4NeutronHPElasticData * theHPElasticData;
    G4NeutronHPInelastic * theHPInelastic;
    G4NeutronHPInelasticData * theHPInelasticData;
    G4NeutronHPFission * theHPFission;
    G4NeutronHPFissionData * theHPFissionData;
    G4NeutronHPCapture * theHPCapture;
    G4NeutronHPCaptureData * theHPCaptureData;

};

// 2002 by J.P. Wellisch

#endif

