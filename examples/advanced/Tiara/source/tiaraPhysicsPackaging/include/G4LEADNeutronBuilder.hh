#ifndef G4LEADNeutronBuilder_h
#define G4LEADNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4Mars5GeV.hh"   
#include "G4NeutronInelasticCrossSection.hh"

class G4LEADNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4LEADNeutronBuilder();
    virtual ~G4LEADNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4Mars5GeV * theModel;    
    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

