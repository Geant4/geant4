#ifndef G4BertiniElasticCollision_h
#define G4BertiniElasticCollision_h 1

#include "G4BertiniModel.hh"

class G4BertiniElasticCollision : public G4BertiniModel {

public:
  G4BertiniElasticCollision();
  ~G4BertiniElasticCollision();
  void interpolateElasticNeutronData(G4int medium, G4int kdd, G4double e);
  void geti(G4double es, G4int npts, G4double e, G4int i);
  void scatteringWithHydrogen(); 

private:
  G4double pt[16];
  G4double col[23];
};

#endif

  












