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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRBiasingRegionInfo.hh
//   Header file of Region Information class to hold parameters of
//   geometry importance biasing.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRBiasingRegionInfo_h
#define GRBiasingRegionInfo_h 1

#include "G4VUserRegionInformation.hh"
#include "globals.hh"

class GRBiasingRegionInfo : public G4VUserRegionInformation
{
  public:
    GRBiasingRegionInfo()
    {;}
    virtual ~GRBiasingRegionInfo()
    {;}

  public:
    virtual void Print() const
    {
      G4cout << "GRBiasingRegionInfo : biasing factor " << factor
             << "  bias probability : " << probability << G4endl;
    }

 public:
    void SetBiasingFactor(G4int val)
    {
      if(val>0)
      {
        factor = val;
        if(val>4)
        {
          G4Exception("G4GRBiasingRegionInfo::SetBiasingFactor","Gorad0005",JustWarning,
                      "Biasing factor maybe too large. Results should be carefully checked.");
        }
      }
      else
      { G4Exception("G4GRBiasingRegionInfo::SetBiasingFactor","Gorad0005",JustWarning,"Invalid value"); }
    }
    inline G4int GetBiasingFactor() const
    { return factor; }
    void SetProbability(G4double val)
    {
      if(val>=0.&&val<=1.)
      { probability = val; }
      else
      { G4Exception("G4GRBiasingRegionInfo::SetProbability","Gorad0006",JustWarning,"Invalid value"); }
    }
    inline G4double GetProbability() const
    { return probability; }

 private:
    G4int factor = 2;
    G4double probability = 1.;
};

#endif

