#ifndef G4NeutronHPThermalScatteringNames_h
#define G4NeutronHPThermalScatteringNames_h 1

// Class Description
// Name list of Elements for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV; 
// Based on Thermal neutron scattering files 
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

// 15-Nov-06 First implementation is done by T. Koi (SLAC/SCCS)

#include "globals.hh"
#include <map> 

class G4NeutronHPThermalScatteringNames 
{
   public:
   
      G4NeutronHPThermalScatteringNames();
   
      ~G4NeutronHPThermalScatteringNames();

      G4bool IsThisThermalElement ( G4String ); 
      size_t GetSize() { return names.size(); };
      G4String GetTS_NDL_Name( G4String nameG4Element ) { return  names.find ( nameG4Element )->second; };
      //G4String GetTS_G4E_Name( G4int i ) { return  names[i]->first; };
   
   private:

//              G4Element  NDL name 
      std::map< G4String , G4String > names;

};

#endif
