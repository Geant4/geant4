#ifndef ANADataPoint_h
#define ANADataPoint_h

#include <iostream>
#include "g4std/vector"

class ANADataPoint
{
  public:
  ANADataPoint(G4double aBin, G4double aWidth, G4double aMass) : 
         highBin(aBin), theMass(aMass), meanKinetic(0), 
	 entries(0), xSec(0), theWidth(aWidth){}
  
  G4bool InsertAt(G4double anEnergy, G4double aXsec)
  {
//    G4cout << "Entering a data point "<< anEnergy << " " << aXsec 
//           << " " << theMass << " " << highBin << G4endl;
    G4bool result = (anEnergy-theMass<highBin);
    if (result)
    {
//      G4cout << "Entering a data point "<< anEnergy << " " << aXsec 
//             << " " << theMass << " " << highBin << G4endl;
//      input.push_back(aXsec);
      xSec += aXsec;
      meanKinetic = (entries*meanKinetic+(anEnergy-theMass))/G4double(++entries);
    }
    return result;
  }
  
  void DumpInfo(ostream & aOStream, G4double aWeight)
  {
    // G4cout << "DataPoint Info "<<meanKinetic/GeV<<" "<<xSec<<" "<<highBin<<" "<<G4endl;
    G4double err=0;
    if(entries) err=xSec/sqrt(G4double(entries));
    G4double energyInGeV = meanKinetic/GeV;
    G4double doubleDiffCrossSectionWithUnits = aWeight*(xSec/theWidth)/(millibarn/GeV/GeV);
    G4double errorWithUnits = aWeight*(err/theWidth)/(millibarn/GeV/GeV);
    aOStream <<energyInGeV<<" "
             << 0.<<" "
             << doubleDiffCrossSectionWithUnits <<" "
	     << errorWithUnits <<G4endl;
//    aOStream <<"points contributing to the above"<<G4endl;
//    G4int i;
//    aOStream << "Total sum = "<<xSec<<G4endl;
//    for(i=0; i<input.size(); i++)
//    {
//      aOStream << "1./(4*pi*p) = "<<input[i]<<G4endl;
//    }
  }
  G4double GetMeanKinetic(){return meanKinetic;}
  G4double GetXSec(){return xSec;}
  
  private:

  G4double highBin;
  G4double theMass;
  G4double theWidth;
//  G4std::vector<G4double> input;

  G4double meanKinetic;
  G4int entries;
  G4double xSec;
};
#endif
