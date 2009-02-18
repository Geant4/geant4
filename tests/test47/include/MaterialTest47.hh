#ifndef MaterialTest47_h
#define MaterialTest47_h 1

#include "globals.hh"

class G4Material;

class MaterialTest47 {

public:
  
  MaterialTest47();
  ~MaterialTest47();
     
  G4Material* GetMaterial(const G4String&);     
	     
};

#endif

 


