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
//
// $Id: G4StableIsotopesTest.cc,v 1.3 2001-07-11 10:08:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "G4StableIsotopes.hh"

main()
{
   G4StableIsotopes theIso;
   for (int Z=1; Z<92; Z++)
   {
     G4int protonCount = theIso.GetProtonCount(Z);
     G4String theName = theIso.GetName(Z);
     for (G4int i=0; i<theIso.GetNumberOfIsotopes(Z); i++)
     {
       G4int nucleons = theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+i);
       G4double fracInPercent=theIso.GetAbundance(theIso.GetFirstIsotope(Z)+i);
       G4cout << nucleons << " " << fracInPercent << " " << protonCount << " " << theName << G4endl;
     }
     G4cout << G4endl;
   }
   G4cout << G4endl;
}
