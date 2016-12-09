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
#ifndef G4AnnihilationCrossSection_h
#define G4AnnihilationCrossSection_h

#include "globals.hh"
#include <vector>  
#include "G4VAnnihilationCrossSection.hh"

class G4AnnihilationCrossSection
{
  public:
    G4AnnihilationCrossSection();

    G4double GetCrossSection(int aCode, int bCode, G4double S);
  private:

    std::vector<G4VAnnihilationCrossSection*> theDataSets;
};

inline G4double G4AnnihilationCrossSection::GetCrossSection(G4int aCode, G4int bCode, G4double S)
{
  G4double result = 0.;
  typedef std::vector<G4VAnnihilationCrossSection*>::iterator iter;
  iter i;
  for(i=theDataSets.begin(); i!=theDataSets.end(); i++)
  {
    if((*i)->InCharge(aCode, bCode))
    {
      result = (*i)->GetXsec(S);
      break;
    }
  }
  return result;
}

#endif

