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
#ifndef G4ASCCrossSection_h
#define G4ASCCrossSection_h

#include "globals.hh"
#include "G4VAnnihilationCrossSection.hh"

class G4ASCCrossSection : public G4VAnnihilationCrossSection
{
  public:
    G4ASCCrossSection(G4int, G4int, G4double,  G4double, G4double, G4double);
    G4bool InCharge(G4int aCode, G4int bCode);
    G4double GetXsec(G4double s);
  private:
    
    G4int theCode1;
    G4int theCode2;
    G4double theX;
    G4double theY;
    G4double theEta;
    G4double theEps;
};


inline G4bool G4ASCCrossSection::
InCharge(G4int aCode, G4int bCode)
{
  G4bool result;
  result = (aCode==theCode1&&bCode==theCode2)||(aCode==theCode2&&bCode==theCode1);
  return result;
}

inline G4double G4ASCCrossSection::
GetXsec(G4double s)
{
   G4double result = theX*pow(s, theEps) + theY*pow(s, -theEta);
   return result;
}

#endif
