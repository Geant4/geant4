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
// $Id: G4VScoreColorMap.hh,v 1.3 2007/11/06 17:17:14 asaim Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#ifndef G4VScoreColorMap_h
#define G4VScoreColorMap_h 1

#include "globals.hh"

class G4VScoreColorMap
{
  public:
      G4VScoreColorMap(G4String mName);
      virtual ~G4VScoreColorMap();

  public:
      virtual void GetMapColor(G4double val, G4double color[4]) = 0;

  public:
      inline G4String GetName() const
      { return fName; }
      inline void SetFloatingMinMax(G4bool vl=true)
      { ifFloat = vl; }
      inline G4bool IfFloatMinMax() const 
      { return ifFloat; }
      inline void SetMinMax(G4double minVal, G4double maxVal)
      {
	if(fMinVal >= fMaxVal)
	{ G4cerr << "G4VScoreColoMap::SetMinMax() : minimum is larger than maximum. Verify values you set, ["
		 << fMinVal << ", " << fMaxVal << "]" << G4endl; }
        else
        { fMinVal = minVal; fMaxVal = maxVal; }
      }
      inline G4double GetMin() const
      { return fMinVal; }
      inline G4double GetMax() const
      { return fMaxVal; }

  protected:
      G4String fName;
      G4bool ifFloat;
      G4double fMinVal;
      G4double fMaxVal;
};

#endif

