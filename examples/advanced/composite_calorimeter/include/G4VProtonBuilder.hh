#ifndef G4VProtonBuilder_h
#define G4VProtonBuilder_h

class G4ProtonInelasticProcess;

class G4VProtonBuilder
{
  public:
    G4VProtonBuilder() {}
    virtual ~G4VProtonBuilder() {}
    virtual void Build(G4HadronElasticProcess & aP) = 0;
    virtual void Build(G4ProtonInelasticProcess & aP) = 0;
};
// 2002 by J.P. Wellisch

#endif
