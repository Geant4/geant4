//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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


