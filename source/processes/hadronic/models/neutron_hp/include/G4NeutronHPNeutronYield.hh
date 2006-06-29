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
// $Id: G4NeutronHPNeutronYield.hh,v 1.10 2006-06-29 20:49:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPNeutronYield_h
#define G4NeutronHPNeutronYield_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPPolynomExpansion.hh"
#include "G4NeutronHPList.hh"

class G4NeutronHPNeutronYield
{
  public:
  G4NeutronHPNeutronYield()
  {
    simpleMean = true;
    spontPrompt = true;
    hasPromptData = false;
    hasDelayedData = false;
  }
  ~G4NeutronHPNeutronYield(){}
  
  G4double GetTargetMass() { return targetMass; }
  
  void InitMean(std::ifstream & aDataFile)
  {
    G4int iflag;
    aDataFile >> targetMass >>iflag;
    if(iflag == 1) simpleMean=false;
    if(simpleMean)
    {
      theSimpleMean.Init(aDataFile, eV);
    }
    else
    {
      theMean.Init(aDataFile);
    }
  }

  void InitPrompt(std::ifstream & aDataFile)
  { 
    hasPromptData = true;
    G4int iflag;
    aDataFile >> targetMass >>iflag;
    if(iflag == 2) spontPrompt = false;
    if(spontPrompt)
    {
      aDataFile >> theSpontPrompt;
    }
    else
    {
      thePrompt.Init(aDataFile, eV);
    }
  }
 
  void InitDelayed(std::ifstream & aDataFile)
  {
    hasDelayedData = true;
    G4int iflag;
    aDataFile >> targetMass >>iflag;
    thePrecursorDecayConstants.Init(aDataFile, 1./s); // s is the CLHEP unit second
    if(iflag == 2) spontDelayed = false;
    if(spontDelayed)
    {
      aDataFile >> theSpontDelayed;
    }
    else
    {
      theDelayed.Init(aDataFile, eV);
    }
  }
  
  G4double GetMean(G4double anEnergy)
  {
    if(simpleMean)
    {
    return theSimpleMean.GetY(anEnergy);
    }
    return theMean.GetValue(anEnergy);
  }
  
  G4double GetPrompt(G4double anEnergy)
  {
    if(!hasPromptData) return 0;
    if(spontPrompt)
    {
      return theSpontPrompt;
    }
    return thePrompt.GetY(anEnergy);
  }
  
  G4double GetDelayed(G4double anEnergy)
  {
    if(!hasDelayedData) return 0;
    if(spontDelayed)
    {
      return theSpontDelayed;
    }
    return theDelayed.GetY(anEnergy);
  }
  
  inline G4double GetDecayConstant(G4int i)
  {
    return thePrecursorDecayConstants.GetValue(i);
  }
  
  private:
  
  G4double targetMass;
  // total mean
  G4bool simpleMean;
  G4NeutronHPPolynomExpansion theMean;
  G4NeutronHPVector theSimpleMean;
  
  // Prompt neutrons
  G4bool hasPromptData;
  G4bool spontPrompt;
  G4NeutronHPVector thePrompt;
  G4double theSpontPrompt;
  
  // delayed neutrons
  G4bool hasDelayedData;
  G4bool spontDelayed;
  G4NeutronHPList thePrecursorDecayConstants;
  G4NeutronHPVector theDelayed;
  G4double theSpontDelayed;

};
#endif
