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
#ifndef G4ping_h
#define G4ping_h

#include <vector>

#include "globals.hh"
#include "G4LorentzVector.hh"

  class G4ping
  {
    public:
    G4ping(G4String aName) : theName(aName) {};
    
    void push_back(G4String aS) {theS.push_back(aS);}
    void push_back(G4double aD) {theD.push_back(aD);}
    void push_back(G4LorentzVector aV) {theV.push_back(aV);}
    
    void dump()
    {
      if(getenv(theName))
      {
        size_t i(0);
	for(i=0; i<theS.size(); i++)
	{
	  G4cout << theS[i]<<", ";
	}
	for(i=0; i<theD.size(); i++)
	{
	  G4cout << theD[i]<<", ";
	}
	for(i=0; i<theV.size(); i++)
	{
	  G4cout << theV[i]<<", ";
	}
	G4cout << G4endl;
      }
      theS.clear();
      theD.clear();
      theV.clear();
    }
    private:
    std::vector<G4String> theS;
    std::vector<G4double> theD;
    std::vector<G4LorentzVector> theV;
    
    G4String theName;
  };
#endif
