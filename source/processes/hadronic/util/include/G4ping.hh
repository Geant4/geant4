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
    vector<G4String> theS;
    vector<G4double> theD;
    vector<G4LorentzVector> theV;
    
    G4String theName;
  };
#endif
