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
#ifndef CCalDataSet_h
#define CCalDataSet_h

#include "G4ParticleHPInterpolator.hh"
#include <vector>

// by JPW, working, but to be cleaned up. @@@@

class CCalDataSet 
{
public:

  void insert(G4double anEnergy, G4double aXsection)
  {
    pair<G4double, G4double> aPoint(anEnergy, aXsection);
    theData.push_back(aPoint);
  }
  
  double getCrossSection(G4double anEnergy)
  {
    G4double result;
    if(anEnergy < theData[0].first)
    {
      result = 0;
    }
    else
    {
      G4double x1,x2;
      G4double y1,y2;
      if(anEnergy > theData[theData.size()-1].second)
      {
        // extrapolation
        G4int n=theData.size();
        x1 = theData[n-1].first;
        y1 = theData[n-1].second;
        x2 = theData[n-2].first;
        y2 = theData[n-2].second;
      }
      else
      {
        // linear search and interpolation
        unsigned int i;
        for(i=0; i<theData.size(); ++i)
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

  G4ParticleHPInterpolator theInterpolator;
  
private:

  vector<pair<G4double, G4double> > theData;
};
#endif
