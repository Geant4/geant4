#ifndef G4VPiKBuilder_h
#define G4VPiKBuilder_h

#include "G4HadronElasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"

class G4VPiKBuilder
{
  public:
    G4VPiKBuilder() {}
    virtual ~G4VPiKBuilder() {}
    virtual void Build(G4HadronElasticProcess & aP) = 0;
    virtual void Build(G4PionPlusInelasticProcess & aP) = 0;
    virtual void Build(G4PionMinusInelasticProcess & aP) = 0;
    virtual void Build(G4KaonPlusInelasticProcess & aP) = 0;
    virtual void Build(G4KaonMinusInelasticProcess & aP) = 0;
    virtual void Build(G4KaonZeroLInelasticProcess & aP) = 0;
    virtual void Build(G4KaonZeroSInelasticProcess & aP) = 0;
};
// 2002 by J.P. Wellisch

#endif
