//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
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
