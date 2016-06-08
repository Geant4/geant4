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
#include "G4MesonSplitter.hh"
#include "Randomize.hh"

G4bool G4MesonSplitter::SplitMeson(G4int PDGcode, G4int* aEnd, G4int* bEnd)
{
  G4int absPDGcode = abs(PDGcode);
  if(absPDGcode == 22)
  {
    G4int it=1;
    if(G4UniformRand()<.5) it++;
    *aEnd = it;
    *bEnd = -it;
    return TRUE; 	 
  }
  if (absPDGcode >= 1000) return FALSE;

  G4int heavy =  absPDGcode/100;
  G4int light = (absPDGcode%100)/10;
  G4int anti  = 1 - 2*(G4std::max(heavy, light)%2);
  if (PDGcode < 0 ) anti = -anti;
  heavy *=  anti;
  light *= -anti;
  if ( anti < 0) 
     G4SwapObj(&heavy, &light);
  *aEnd = heavy;
  *bEnd = light;
  return TRUE; 	 
}
