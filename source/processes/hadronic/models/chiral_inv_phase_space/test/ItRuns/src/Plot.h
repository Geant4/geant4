#ifndef ANAPlot_h
#define ANAPlot_h

#include "g4std/vector"
#include "Analysis/src/DataPoint.h"
#include "Analysis/src/TVPlot.h"

template <class DataPointType>
class ANAPlot : public TVANAPlot<DataPointType>
{
  public:
    ANAPlot(G4int aPDG, G4double aMass, G4double aTotalXsec, G4String fn = "");
  
    G4bool Insert(G4int aPDG, G4double anEnergy, G4double aXsec)
    {
      if(aPDG==thePDG)
      {
        for(G4int i=0; i<theDataPoints.size(); i++)
	{
	  if(theDataPoints[i].InsertAt(anEnergy, aXsec)) break;
	}
      }
      return aPDG==thePDG;
    }
    
    void DumpInfo(ostream & aOStream)
    {
      G4cout <<" Dumping info for PDG = "<<thePDG<<G4endl;
      ofstream theOutput(theOutputFile);
      for(G4int i=0; i<theDataPoints.size(); i++)
      {
        theDataPoints[i].DumpInfo(aOStream);
        if(theOutputFile!="") theDataPoints[i].DumpInfo(theOutput);
      }
    }
    
  private:
    G4int thePDG;
    G4double theTotalXsec;
    
    G4std::vector<DataPointType> theDataPoints;
    G4String theOutputFile;
};

template<class DataPointType>
ANAPlot<DataPointType>::ANAPlot(G4int aPDG, G4double aMass, G4double aTotalXsec, G4String fn) : 
         thePDG(aPDG), theTotalXsec(aTotalXsec)
{
  theOutputFile = fn;
  theDataPoints.push_back(DataPointType(0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.075*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.125*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.15*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.175*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.2*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.25*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.3*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.35*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.4*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.45*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.5*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.55*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.6*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.7*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.8*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.9*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.0*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.2*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.3*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.4*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.5*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.6*GeV, aMass));
  theDataPoints.push_back(DataPointType(2.0*GeV, aMass));
  theDataPoints.push_back(DataPointType(2.5*GeV, aMass));
}

#endif
