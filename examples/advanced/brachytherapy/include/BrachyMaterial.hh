#ifndef BrachyMaterial_H
#define BrachyMaterial_H 1
#include "globals.hh"
class G4Material;


class BrachyMaterial
{ public:
  BrachyMaterial();
  ~ BrachyMaterial();

public:
  void  DefineMaterials();

public:
  G4Material* GetMat(G4String);
};

#endif
