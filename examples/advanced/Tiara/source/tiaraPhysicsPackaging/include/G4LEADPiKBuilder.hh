#ifndef G4LEADPiKBuilder_h
#define G4LEADPiKBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4Mars5GeV.hh"   

class G4LEADPiKBuilder : public G4VPiKBuilder
{
  public: 
    G4LEADPiKBuilder();
    virtual ~G4LEADPiKBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4PionPlusInelasticProcess & aP);
    virtual void Build(G4PionMinusInelasticProcess & aP);
    virtual void Build(G4KaonPlusInelasticProcess & aP);
    virtual void Build(G4KaonMinusInelasticProcess & aP);
    virtual void Build(G4KaonZeroLInelasticProcess & aP);
    virtual void Build(G4KaonZeroSInelasticProcess & aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4PiNuclearCrossSection thePiData;
    G4Mars5GeV * theModel;    
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

