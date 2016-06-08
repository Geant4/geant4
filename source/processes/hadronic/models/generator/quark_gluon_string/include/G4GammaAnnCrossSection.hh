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
#ifndef G4GammaAnnCrossSection_h
#define G4GammaAnnCrossSection_h

#include "globals.hh"
#include "g4std/vector"
#include "G4VAnnihilationCrossSection.hh"
#include "G4ASCCrossSection.hh"

class G4GammaAnnCrossSection : public G4VAnnihilationCrossSection
{
  public:
    G4GammaAnnCrossSection();
    G4bool InCharge(G4int aCode, G4int bCode);
    G4double GetXsec(G4double s);
    virtual ~G4GammaAnnCrossSection(){}
    
  private:
  
    G4std::vector<G4ASCCrossSection*> theGammaNucXSections;
};

#endif
