#ifndef G4LHEPProtonBuilder_h
#define G4LHEPProtonBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4LElastic.hh"   
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"

class G4LHEPProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4LHEPProtonBuilder();
    virtual ~G4LHEPProtonBuilder();

  public: 
    virtual void Build(G4ProtonInelasticProcess & aP);
    virtual void Build(G4HadronElasticProcess & aP);
    
    void SetMinEnergy(G4double aM) 
    {
      theMin = aM;
    }

  private:
    G4LElastic * theElasticModel;
    G4LEProtonInelastic * theLEProtonModel;
    G4HEProtonInelastic * theHEProtonModel;
    
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

