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
#ifndef CCalDataSet_h
#define CCalDataSet_h

#include "G4NeutronHPInterpolator.hh"
#include "g4std/vector"

// by JPW, working, but to be cleaned up. @@@@

class CCalDataSet 
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
