#ifndef G4LEPNeutronBuilder_h
#define G4LEPNeutronBuilder_h 1

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

class G4LEPNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4LEPNeutronBuilder();
    virtual ~G4LEPNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) 
    {
      theMin=aM;
      theIMin = theMin;
    }
    void SetMinInelasticEnergy(G4double aM) 
    {
      theIMin=aM;
    }
    void SetMaxEnergy(G4double aM) 
    {
      theMax=aM;
    }
    
  private:
    G4double theMin;
    G4double theIMin;
    G4double theMax;
    G4LElastic * theElasticModel;
    G4LENeutronInelastic * theLENeutronModel;
    G4LFission * theNeutronFissionModel;
    G4LCapture * theNeutronCaptureModel;

};

// 2002 by J.P. Wellisch

#endif

