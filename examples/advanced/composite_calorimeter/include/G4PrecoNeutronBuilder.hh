#ifndef G4PrecoNeutronBuilder_h
#define G4PrecoNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4PreCompoundModel.hh" 
#include "G4ExcitationHandler.hh"  
#include "G4NeutronInelasticCrossSection.hh"

class G4PrecoNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4PrecoNeutronBuilder();
    virtual ~G4PrecoNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4HadronFissionProcess & aP);
    virtual void Build(G4HadronCaptureProcess & aP);
    virtual void Build(G4NeutronInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4PreCompoundModel * theModel;   
    G4ExcitationHandler theHandler;
    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;
    G4double theMax;

};

// 2002 by J.P. Wellisch

#endif

