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
#ifndef G4IsoResult_h
#define G4IsoResult_h 

#include "globals.hh"
#include "G4Nucleus.hh"
// Class Description
// Communication class for isotope production; 
// Only interesting, if you want to implement your own isotope production model.
// Class Description - End


class G4IsoResult
{

public:

  G4IsoResult(const G4String & anIso, const G4Nucleus & aMother)
  {
    theIsoName = anIso;
    theMotherNucleus = aMother;
  }
  
  ~G4IsoResult();
  
  G4String GetIsotope() { return theIsoName; }
  G4Nucleus GetMotherNucleus() { return theMotherNucleus; }
  
private:
  G4IsoResult() {}
  G4IsoResult & operator = (const G4IsoResult & ) { return *this;}

private:
  G4String theIsoName;
  G4Nucleus theMotherNucleus;

};

#endif
