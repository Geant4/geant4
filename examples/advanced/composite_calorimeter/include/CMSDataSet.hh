#ifndef CMSDataSet_h
#define CMSDataSet_h

#include "G4NeutronHPInterpolator.hh"
#include <vector>
#include <utility>

// by JPW, working, but to be cleaned up. @@@@

class CMSDataSet 
{
public:

  void insert(double anEnergy, double aXsection)
  {
    pair<double, double> aPoint(anEnergy, aXsection);
    theData.push_back(aPoint);
  }
  
  double getCrossSection(double anEnergy)
  {
    double result;
    if(anEnergy < theData[0].first)
    {
      result = 0;
    }
    else
    {
      double x1,x2;
      double y1,y2;
      if(anEnergy > theData[theData.size()-1].second)
      {
        // extrapolation
        int n=theData.size();
        x1 = theData[n-1].first;
        y1 = theData[n-1].second;
        x2 = theData[n-2].first;
        y2 = theData[n-2].second;
      }
      else
      {
        // linear search and interpolation
        unsigned int i;
        for(i=0; i<theData.size(); i++)
        {
          if(theData[i].first>anEnergy) break;
        }
        x1 = theData[i-1].first;
        y1 = theData[i-1].second;
        x2 = theData[i].first;
        y2 = theData[i].second;
      }
      // interpolate logarithmic in energy and linear in cross-section
      result = theInterpolator.Interpolate(LOGLIN, anEnergy,  x1,x2,y1,y2);
    }
    return result*barn;
  }

private:

  G4NeutronHPInterpolator theInterpolator;
  
private:

  vector<pair<double, double> > theData;
};
#endif
