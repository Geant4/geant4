#ifndef G4BertiniProtonBuilder_h
#define G4BertiniProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4CascadeInterface.hh"   
#include "G4ProtonInelasticCrossSection.hh"

class G4BertiniProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4BertiniProtonBuilder();
    virtual ~G4BertiniProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4CascadeInterface * theModel;    
    G4double theMin;
    G4double theMax;

};

// 2002 by J.P. Wellisch

#endif

