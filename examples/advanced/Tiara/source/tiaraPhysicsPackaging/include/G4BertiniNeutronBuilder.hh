#ifndef G4BertiniNeutronBuilder_h
#define G4BertiniNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4CascadeInterface.hh"   
#include "G4NeutronInelasticCrossSection.hh"

class G4BertiniNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4BertiniNeutronBuilder();
    virtual ~G4BertiniNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}

  private:
    G4CascadeInterface * theModel;    
    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;
    G4double theMax;

};

// 2002 by J.P. Wellisch

#endif

