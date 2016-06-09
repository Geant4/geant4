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
#ifndef G4ping_h
#define G4ping_h

#include "G4String.hh"
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
