#ifndef G4LHEPNeutronBuilder_h
#define G4LHEPNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4LElastic.hh"   
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LENeutronInelastic.hh"
#include "G4HENeutronInelastic.hh"

class G4LHEPNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4LHEPNeutronBuilder();
    virtual ~G4LHEPNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) 
    {
      theMin = aM;
      theIMin = aM;
    }
    void SetMinInelasticEnergy(G4double aM) 
    {
      theIMin = aM;
    }

  private:
    G4LElastic * theElasticModel;
    G4LENeutronInelastic * theLENeutronModel;
    G4HENeutronInelastic * theHENeutronModel;
    G4LFission * theNeutronFissionModel;
    G4LCapture * theNeutronCaptureModel;
    
    G4double theMin;
    G4double theIMin;

};

// 2002 by J.P. Wellisch

#endif

