#ifndef ANADataPoint_h
#define ANADataPoint_h

#include <iostream>

class ANADataPoint
{
  public:
  ANADataPoint(G4double aBin, G4double aMass) : 
         highBin(aBin), theMass(aMass), meanKinetic(0), entries(0), xSec(0) {}
  
  G4bool InsertAt(G4double anEnergy, G4double aXsec)
  {
//    G4cout << "Entering a data point "<< anEnergy << " " << aXsec 
//           << " " << theMass << " " << highBin << G4endl;
    G4bool result = (anEnergy-theMass<highBin);
    if (result)
    {
      G4cout << "Entering a data point "<< anEnergy << " " << aXsec 
             << " " << theMass << " " << highBin << G4endl;
      xSec += aXsec;
      meanKinetic = (entries*meanKinetic+(anEnergy-theMass))/G4double(++entries);
    }
    return result;
  }
  
  void DumpInfo(ostream & aOStream)
  {
    // G4cout << "DataPoint Info "<<meanKinetic/GeV<<" "<<xSec<<" "<<highBin<<" "<<G4endl;
    aOStream <<meanKinetic/GeV<<" "<<xSec<<G4endl;
  }
  G4double GetMeanKinetic(){return meanKinetic;}
  G4double GetXSec(){return xSec;}
  
  private:

  G4double highBin;
  G4double theMass;

  G4double meanKinetic;
  G4int entries;
  G4double xSec;
};
#endif
