#ifndef exrdmEnergyDeposition_h
#define exrdmEnergyDeposition_h 1

#include "globals.hh"

class exrdmEnergyDeposition
{
  public:   // with description

    exrdmEnergyDeposition();
    exrdmEnergyDeposition( const exrdmEnergyDeposition &right );
    exrdmEnergyDeposition( G4double, G4double, G4double );
    virtual ~exrdmEnergyDeposition();
         // Constructor and virtual destructor

    G4bool operator==(const exrdmEnergyDeposition &right) const ;
    G4bool operator< (const exrdmEnergyDeposition &right) const ;
    G4bool operator<=(const exrdmEnergyDeposition &right) const ;
  // Operators  

  G4double GetEnergy() {return Energy;};
  G4double GetTime() {return Time;};
  G4double GetWeight() {return Weight;};
  // Accessors

  private:

    G4double Energy;  
    G4double Time;    
    G4double Weight;
};
#endif



