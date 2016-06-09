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
#ifndef G4NeutronHPHash_h
#define G4NeutronHPHash_h

#include "globals.hh"
#include <vector>

class G4NeutronHPHash
{
public:
  G4NeutronHPHash() 
  {
    theUpper = 0;
    prepared = false;
  }
  
  ~G4NeutronHPHash()
  {
    if(theUpper) delete theUpper;
  }
  
  G4NeutronHPHash(const G4NeutronHPHash & aHash)
  {
    theIndex = aHash.theIndex;
    theData = aHash.theData;
    if(aHash.theUpper != 0)
    {
      theUpper = new G4NeutronHPHash(*(aHash.theUpper));
    }
    else
    {
      theUpper = 0;
    }
  }
  
  G4NeutronHPHash & operator = (const G4NeutronHPHash & aHash)
  {
    if(&aHash != this)
    {
      theIndex = aHash.theIndex;
      theData = aHash.theData;
      if(aHash.theUpper != 0)
      {
        theUpper = new G4NeutronHPHash(*(aHash.theUpper));
      }
      else
      {
        theUpper = 0;
      }
    }
    return *this;
  }
  
  void Clear()
  {
    if(theUpper) 
    {
      delete theUpper;
      theUpper = 0;
    }
    theIndex.clear();
    theData.clear();
    prepared = false;
  }
  
  G4bool Prepared() const {return prepared;}
  inline void SetData(G4int index, G4double x, G4double y) 
  { 
    prepared = true;
    G4NeutronHPDataPoint aPoint;
    aPoint.SetData(x, y);
    theData.push_back(aPoint);
    theIndex.push_back(index);
    if(0 == theData.size()%10 && 0!=theData.size())
    {
      if(0 == theUpper) theUpper = new G4NeutronHPHash();
      theUpper->SetData(theData.size()-1, x, y);
    }
  }
    
  G4int GetMinIndex(G4double e) const
  {
    G4int result=-1;
    if(theData.size() == 0) return 0;
    if(theData[0].GetX()>e) return 0;
    
    G4int lower=0;
    if(theUpper != 0)
    {
      lower = theUpper->GetMinIndex(e);
    }
    unsigned int i;
    for(i=lower; i<theData.size(); i++)
    {
      if(theData[i].GetX()>e)
      {
        result = theIndex[i-1];
	break;
      }
    }
    if(result == -1) result = theIndex[theIndex.size()-1];
    return result;
  }
  
private:

  G4bool prepared;
  G4NeutronHPHash * theUpper;
  std::vector<int> theIndex;
  std::vector<G4NeutronHPDataPoint> theData; // the data
};
#endif
