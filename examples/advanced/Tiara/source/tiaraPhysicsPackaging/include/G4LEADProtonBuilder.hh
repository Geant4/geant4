#ifndef G4LEADProtonBuilder_h
#define G4LEADProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4Mars5GeV.hh"   
#include "G4ProtonInelasticCrossSection.hh"

class G4LEADProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4LEADProtonBuilder();
    virtual ~G4LEADProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4Mars5GeV * theModel;    
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

