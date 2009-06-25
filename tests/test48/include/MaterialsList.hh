#ifndef MaterialsList_h
#define MaterialsList_h 1

#include "globals.hh"

class G4Material;

class MaterialsList {

public:
  
  MaterialsList();
  ~MaterialsList();
     
  G4Material* GetMaterial(const G4String&);     
	     
};

#endif

 


