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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
  G4IsoResult & operator = (const G4IsoResult & aResult) { return *this;}

private:
  G4String theIsoName;
  G4Nucleus theMotherNucleus;

};

#endif
