#ifndef MLColour_HH
#define MLColour_HH
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4Colour.hh"
#include <vector>

class MLColourMessenger;
////////////////////////////////////////////////////////////////////////////////
//
class MLColour
{
public:
  
  MLColour();
  ~MLColour();

public:
  
  void AddColour (G4String, G4double, G4double, G4double );
  G4Colour GetColour (G4int i)  {return Colour[i];};
  G4int GetNbOfColours () {return Colour.size();};
  void ListColour();

private:              

  MLColourMessenger* colourMessenger;

  std::vector<G4Colour> Colour;
  std::vector<G4String> CName;

};
////////////////////////////////////////////////////////////////////////////////
#endif
