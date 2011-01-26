#ifndef MATERIALS_HH
#define MATERIALS_HH

#include "globals.hh"

class G4Material;


class Materials {

 public:
   virtual ~Materials();
   static Materials* Instance();
   virtual void Destroy();
   
   G4Material* GetMaterial(G4String matName);
   
 protected:
   // virtual ~Materials();
   Materials();

   Materials(const Materials& only);
   const Materials& operator=(const Materials& only);

 private:
   static Materials* instance;

   G4Material* hydrogen;
   G4Material* beryllium;
   G4Material* graphite; 
   G4Material* magnesium;
   G4Material* aluminium;
   G4Material* silicon;
   G4Material* liquidArgon;  
   G4Material* titanium;
   G4Material* iron; 
   G4Material* cobalt;
   G4Material* nickel;
   G4Material* indium;
   G4Material* tin; 
   G4Material* copper; 
   G4Material* zinc;  
   G4Material* gallium;
   G4Material* germanium;
   G4Material* zirconium;
   G4Material* molybdenium;
   G4Material* silver;
   G4Material* cadmium;
   G4Material* cesium; 
   G4Material* samarium;
   G4Material* ytterbium; 
   G4Material* tantalum;
   G4Material* tungsten;
   G4Material* gold; 
   G4Material* lead;
   G4Material* uranium;
   G4Material* water; 
   G4Material* quartz; 
   G4Material* ossigeno;
   G4Material* air; 
   G4Material* vacuum;
   G4Material* nytrogen;
};

#endif // MATERIALS_HH
