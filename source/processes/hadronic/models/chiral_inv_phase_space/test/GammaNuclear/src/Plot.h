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

//#define debug

#include "g4std/vector"
#include "GammaNuclear/src/DataPoint.h"
#include "GammaNuclear/src/TVPlot.h"

template <class DataPointType, class FilterType>
class ANAPlot : public TVANAPlot<DataPointType>
{
  public:
    ANAPlot(G4int aPDG, G4double aMass, G4double aTotalXsec, 
           G4String fn = "", FilterType * aFilter = 0);
  
    void SetNevents(G4int aNumber) {theStatistics = aNumber;}
    
    //            ParticlePDGF, ParticleEnergy, WeightOfEvent
    G4bool Insert(G4int aPDG, G4double anEnergy, G4double aXsec)
    {
      if(aPDG==thePDG)
        for(unsigned i=0; i<theDataPoints.size(); i++)
	      if(theDataPoints[i].InsertAt(anEnergy, aXsec)) break;
      return aPDG==thePDG;
    }
    
    void DumpInfo(ostream & aOStream, G4String aPreFix)
    {
#ifdef debug
      G4cout <<" Dumping info for PDG = "<<thePDG<<G4endl;
#endif
      G4String dot = ".";
      G4String it = theOutputFile+theFilter->GetName()+dot+aPreFix;
#ifdef debug
      G4cout<< it<<G4endl;
#endif
      ofstream theOutput(it);
      G4double aWeight = theTotalXsec/theStatistics/theFilter->RelativeGeometricalAcceptance();
      for(unsigned i=0; i<theDataPoints.size(); i++)
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
      G4double theta = aPart->GetCosTheta();
      return theFilter->Accept(theta);
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
  static const G4int nBins = 33;
  theStatistics = 10000;
  theOutputFile = fn;

  G4double hB=0.155;                   // Instead of this one can use the hBins array
  G4double lowE=0.;                    // Value of the first low bin for Energy
  for (G4int i=0; i<nBins; i++)
  {
    hB+=0.005;                         // Increment in case of the equidistant bins
    G4double k = hB*GeV;               // protoQuark Energy
    G4double EpP = k+k+aMass;          // E+p
    G4double highE = (aMass*aMass+EpP*EpP)/(EpP+EpP); // HighBin for the energy
    theDataPoints.push_back(DataPointType(k, highE-lowE, aMass));
    lowE=highE;                        // HighBin becomes LowBin for the next Bin
  }
}

#endif
