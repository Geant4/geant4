#ifndef G4VNeutronBuilder_h
#define G4VNeutronBuilder_h

class G4HadronElasticProcess;
class G4HadronFissionProcess;
class G4HadronCaptureProcess;
class G4NeutronInelasticProcess;

class G4VNeutronBuilder
{
  public:
    G4VNeutronBuilder() {}
    virtual ~G4VNeutronBuilder() {}
    virtual void Build(G4HadronElasticProcess & aP) = 0;
    virtual void Build(G4HadronFissionProcess & aP) = 0;
    virtual void Build(G4HadronCaptureProcess & aP) = 0;
    virtual void Build(G4NeutronInelasticProcess & aP) = 0;
};
// 2002 by J.P. Wellisch

#endif
