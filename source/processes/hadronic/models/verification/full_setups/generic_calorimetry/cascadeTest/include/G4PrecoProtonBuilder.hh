#ifndef G4PrecoProtonBuilder_h
#define G4PrecoProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4ExcitationHandler.hh"  
#include "G4PreCompoundModel.hh"   
#include "G4ProtonInelasticCrossSection.hh"

class G4PrecoProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4PrecoProtonBuilder();
    virtual ~G4PrecoProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4ProtonInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4ProtonInelasticCrossSection theXSec;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel * theModel;    
    G4double theMin;
    G4double theMax;

};

// 2002 by J.P. Wellisch

#endif

