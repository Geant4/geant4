////////////////////////////////////////////////////////////////////////////////
//
#include "MLColour.hh"
#include "MLColourMessenger.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include <vector>
#include <iomanip>  
////////////////////////////////////////////////////////////////////////////////
//
MLColour::MLColour()
{
  Colour.clear();
  
  static G4bool bcol = false;
  if (!bcol) {
    // some default coulours
    G4Colour colour = G4Colour (1.0, 1.0, 1.0) ; // white
    Colour.push_back(colour); 
    CName.push_back("white");
    colour = G4Colour (0.5, 0.5, 0.5) ;  //grey 
    Colour.push_back(colour); 
    CName.push_back("grey");
    colour = G4Colour (.75, .75, .75) ; // dark grey
    Colour.push_back(colour); 
    CName.push_back("dark_grey");
    colour = G4Colour (1.0, 0.0, 0.0) ; // red
    Colour.push_back(colour); 
    CName.push_back("red");
    colour = G4Colour (0.0, 0.0, 1.0) ; // blue
    Colour.push_back(colour); 
    CName.push_back("blue");
    colour = G4Colour (0.0, 1.0, 1.0) ; // cyan
    Colour.push_back(colour); 
    CName.push_back("cyan");
    colour = G4Colour (1.0, 0.0, 1.0) ; // magenta 
    Colour.push_back(colour); 
    CName.push_back("magenta");
    colour = G4Colour (1.0, 1.0, 0.0) ; // yellow
    Colour.push_back(colour); 
    CName.push_back("yellow");
    colour = G4Colour (0.0, 0.0, .75) ; // dark blue
    Colour.push_back(colour); 
    CName.push_back("dark_blue");
    colour = G4Colour (0.0, 1.0, 0.0) ; // green
    Colour.push_back(colour); 
    CName.push_back("green");
    colour = G4Colour (0.0, .75, 0.0) ; // dark green
    Colour.push_back(colour); 
    CName.push_back("dark_green");
    colour = G4Colour (0.75, 0.0, 0.0) ; // dark red
    Colour.push_back(colour); 
    CName.push_back("dark_red");
  
    colourMessenger = new MLColourMessenger(this);
    
    bcol = true;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
MLColour::~MLColour()
{ 
  delete colourMessenger;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLColour::AddColour(G4String name, G4double red, G4double green,
  G4double blue)
{
  G4Colour colour = G4Colour (red, green, blue);
  Colour.push_back(colour); 
  CName.push_back(name);
}
////////////////////////////////////////////////////////////////////////////////
//
void MLColour::ListColour ()
{
  G4cout <<" There are " <<std::setw(3) <<Colour.size() <<" colour defined "
         <<G4endl;
  G4cout <<G4endl;
  for (size_t i = 0; i< Colour.size(); i++) 
    G4cout <<"     Colour Index " <<std::setw(3) <<i+1 <<" "
           <<std::setw(8) <<CName[i] <<G4endl;
  G4cout <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
