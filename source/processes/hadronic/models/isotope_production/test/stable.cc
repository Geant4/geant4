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
#include "G4StableIsotopes.hh"
#include "globals.hh"
#include "g4std/iostream"

int main()
{
// To loop over all isotopes for element with protonount Z
 G4int Z;
 G4StableIsotopes theIso;
 for(Z=1; Z<93; Z++)
 {
   for (G4int i=0; i<theIso.GetNumberOfIsotopes(Z); i++)
   {
      cout <<theIso.GetName(Z)<<" ";
      cout <<Z<<" ";
      cout <<theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+i)<<" ";
      cout << G4endl;
    }
    cout << G4endl;
  }
}
