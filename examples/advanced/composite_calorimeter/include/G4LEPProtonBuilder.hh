#ifndef G4LEPProtonBuilder_h
#define G4LEPProtonBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4LElastic.hh"   
#include "G4LEProtonInelastic.hh"

class G4LEPProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4LEPProtonBuilder();
    virtual ~G4LEPProtonBuilder();

  public: 
    virtual void Build(G4ProtonInelasticProcess & aP);
    virtual void Build(G4HadronElasticProcess & aP);
    
    void SetMinEnergy(G4double aM) 
    {
      theMin=aM;
    }
    void SetMaxEnergy(G4double aM) 
    {
      theMax=aM;
    }
    
  private:
    G4double theMin;
    G4double theMax;
    G4LElastic * theElasticModel;
    G4LEProtonInelastic * theLEProtonModel;

};

// 2002 by J.P. Wellisch

#endif

