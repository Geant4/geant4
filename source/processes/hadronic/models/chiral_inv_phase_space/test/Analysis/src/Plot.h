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
#ifndef ANAPlot_h
#define ANAPlot_h

#include "g4std/vector"
#include "Analysis/src/DataPoint.h"
#include "Analysis/src/TVPlot.h"

template <class DataPointType, class FilterType>
class ANAPlot : public TVANAPlot<DataPointType>
{
  public:
    ANAPlot(G4int aPDG, G4double aMass, G4double aTotalXsec, 
           G4String fn = "", FilterType * aFilter = 0);
  
    void SetNevents(G4int aNumber) {theStatistics = aNumber;}
    
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
    
    void DumpInfo(ostream & aOStream, G4String aPreFix)
    {
      G4cout <<" Dumping info for PDG = "<<thePDG<<G4endl;
      G4String dot = ".";
      G4String it = theOutputFile+theFilter->GetName()+dot+aPreFix;
      G4cout<< it<<G4endl;
      ofstream theOutput(it);
      G4double aWeight = theTotalXsec/theStatistics/theFilter->RelativeGeometricalAcceptance();
      for(G4int i=0; i<theDataPoints.size(); i++)
      {
//        theDataPoints[i].DumpInfo(aOStream);
        if(theOutputFile!="") 
	{
	  theDataPoints[i].DumpInfo(theOutput, aWeight);
	}
      }
    }
    
    G4bool Filter(ANAParticle * aPart) 
    {
      G4double cosTheta = aPart->GetCosTheta();
      return theFilter->Accept(cosTheta);
    }
    
  private:
    G4int thePDG;
    G4double theTotalXsec;
    G4int theStatistics;
    FilterType * theFilter;
    
    G4std::vector<DataPointType> theDataPoints;
    G4String theOutputFile;
};

template<class DataPointType, class FilterType>
ANAPlot<DataPointType, FilterType>::ANAPlot(G4int aPDG, G4double aMass, G4double aTotalXsec,
                                            G4String fn, FilterType * aFilter) : 
         thePDG(aPDG), theTotalXsec(aTotalXsec), theFilter(aFilter)
{
  theStatistics = 10000;
  theOutputFile = fn;
  theDataPoints.push_back(DataPointType(0.05*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.075*GeV, 0.025*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.1*GeV, 0.025*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.125*GeV, 0.025*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.15*GeV, 0.025*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.175*GeV, 0.025*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.2*GeV, 0.025*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.25*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.3*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.35*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.4*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.45*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.5*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.55*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.6*GeV, 0.05*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.7*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.8*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(0.9*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.0*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.1*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.2*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.3*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.4*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.5*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(1.6*GeV, 0.1*GeV, aMass));
  theDataPoints.push_back(DataPointType(2.0*GeV, 0.4*GeV, aMass));
  theDataPoints.push_back(DataPointType(2.5*GeV, 0.5*GeV, aMass));
}

#endif
