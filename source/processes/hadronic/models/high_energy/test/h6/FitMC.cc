#include "G4HEPlot.hh"
#include "globals.hh"

int main()
{
  G4int i;

  G4HEPlot* plot = new G4HEPlot[39];

  plot[2].GetFromFile(2, "Proton.plot");
  plot[2].Print("hist",102);

  for(i=10; i<=15; i++) 
    { 
      plot[i].GetFromFile(i, "Proton.plot");      
      plot[i].Print("hist",100+i);
    }
  for(i=24; i<=29; i++) 
    { 
      plot[i].GetFromFile(i, "Proton.plot");      
      plot[i].Print("hist",100+i);
    }
  for(i=32; i<=38; i++) 
    { 
      plot[i].GetFromFile(i, "Proton.plot");      
      plot[i].Print("hist",100+i);
    }
  delete plot;
}


