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
#ifndef ANADataPoint_h
#define ANADataPoint_h

#include <iostream>
#include "g4std/vector"

class ANADataPoint
{
public:
//           bin highEdge, WidthForBin Norm, MassOfOutParticle
  ANADataPoint(G4double aBin, G4double aWidth, G4double aMass) : 
    highBin(aBin), theMass(aMass), meanK(0), entries(0), xSec(0), theWidth(aWidth){}

  //              TotEnergyOfPart  , WeightOfEvent
  G4bool InsertAt(G4double anEnergy, G4double aXsec)
  {
    G4double kinE=anEnergy-theMass;
    G4double k=(kinE+sqrt(kinE*(anEnergy+theMass)))/2.; // protoQuark's energy/momentum
//
//    G4cout << "Check the data point E="<< anEnergy << ", Xs=" << aXsec 
//           << ", M=" << theMass << ", k=" << k << ", hB=" << highBin << G4endl;
//
    G4bool result = (k<highBin);
    if (result)
    {
//
//    G4cout << "Fill the data point E="<< anEnergy << ", Xs=" << aXsec 
//           << ", M=" << theMass << ", k=" << k << ", hB=" << highBin << G4endl;
//      input.push_back(aXsec);
      xSec += aXsec; // Collect total weight for the bin
      meanK = (entries*meanK+k)/G4double(++entries);
    }
    return result;
  }
  
  void DumpInfo(ostream & aOStream, G4double aWeight)
  {
    // G4cout << "DataPoint Info "<<meanK/GeV<<" "<<xSec<<" "<<highBin<<" "<<G4endl;
    G4double err=0;
    if(entries) err=xSec/sqrt(G4double(entries));
    G4double energyInGeV = meanK/GeV;
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
  G4double GetMeanK(){return meanK;}
  G4double GetXSec(){return xSec;}
  
  private:

  G4double highBin;
  G4double theMass;
//  G4std::vector<G4double> input;
  G4double meanK;
  G4int entries;
  G4double xSec;
  G4double theWidth;
};
#endif
