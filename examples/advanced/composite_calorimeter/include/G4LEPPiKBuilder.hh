#ifndef G4LEPPiKBuilder_h
#define G4LEPPiKBuilder_h

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPiKBuilder.hh"

#include "G4LElastic.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"

class G4LEPPiKBuilder : public G4VPiKBuilder
{
  public: 
    G4LEPPiKBuilder();
    virtual ~G4LEPPiKBuilder();
    
  public: 
    virtual void Build(G4HadronElasticProcess & aP);
    virtual void Build(G4PionPlusInelasticProcess & aP);
    virtual void Build(G4PionMinusInelasticProcess & aP);
    virtual void Build(G4KaonPlusInelasticProcess & aP);
    virtual void Build(G4KaonMinusInelasticProcess & aP);
    virtual void Build(G4KaonZeroLInelasticProcess & aP);
    virtual void Build(G4KaonZeroSInelasticProcess & aP);
    
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
    
    G4LEPionPlusInelastic* theLEPionPlusModel;
    G4LEPionMinusInelastic* theLEPionMinusModel;
    G4LEKaonPlusInelastic* theLEKaonPlusModel;
    G4LEKaonMinusInelastic* theLEKaonMinusModel;
    G4LEKaonZeroLInelastic* theLEKaonZeroLModel;
    G4LEKaonZeroSInelastic* theLEKaonZeroSModel;

};
// 2002 by J.P. Wellisch

#endif
