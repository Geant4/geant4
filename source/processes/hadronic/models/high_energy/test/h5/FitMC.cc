#include "G4HEPlot.hh"
#include "globals.hh"

int main()
{
  G4int i;

  G4HEPlot* plot = new G4HEPlot[5];

  for(i=1; i<=4; i++) 
    { 
      plot[i].GetFromFile(i, "ProtonAu.plot");      
      plot[i].Print("hist",100+i);
    }
  delete plot;
}


