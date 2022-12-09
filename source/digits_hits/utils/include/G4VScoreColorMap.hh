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
//
//

#ifndef G4VScoreColorMap_h
#define G4VScoreColorMap_h 1

#include "globals.hh"

class G4VVisManager;

class G4VScoreColorMap
{
 public:
  G4VScoreColorMap(G4String mName);
  virtual ~G4VScoreColorMap() = default;

 public:
  virtual void GetMapColor(G4double val, G4double color[4]) = 0;

 public:
  inline G4String GetName() const { return fName; }
  inline void SetFloatingMinMax(G4bool vl = true) { ifFloat = vl; }
  inline G4bool IfFloatMinMax() const { return ifFloat; }
  inline void SetMinMax(G4double minVal, G4double maxVal)
  {
    if(minVal >= maxVal)
    {
      G4cerr << "WARNING: G4VScoreColoMap::SetMinMax() : minimum is larger "
                "than or equal to maximum. Verify values you set, ["
             << minVal << ", " << maxVal << "]" << G4endl;
      fMinVal = maxVal;
      fMaxVal = minVal;
    }
    else
    {
      fMinVal = minVal;
      fMaxVal = maxVal;
    }
  }
  inline G4double GetMin() const { return fMinVal; }
  inline G4double GetMax() const { return fMaxVal; }

  // draw a color chart
  virtual void DrawColorChart(G4int nPoint = 5);

  virtual void DrawColorChartBar(G4int nPoint);

  virtual void DrawColorChartText(G4int nPoint);

  void SetPSUnit(G4String& unit) { fPSUnit = unit; }
  void SetPSName(G4String& psName) { fPSName = psName; }

 protected:
  G4String fName;
  G4bool ifFloat = true;
  G4double fMinVal = 0.0;
  G4double fMaxVal = DBL_MAX;
  G4VVisManager* fVisManager = nullptr;
  G4String fPSUnit = "";
  G4String fPSName = "";
};

#endif
