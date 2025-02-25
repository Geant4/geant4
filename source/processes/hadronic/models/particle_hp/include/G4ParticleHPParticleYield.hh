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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPParticleYield_h
#define G4ParticleHPParticleYield_h 1

#include "G4ParticleHPList.hh"
#include "G4ParticleHPPolynomExpansion.hh"
#include "G4ParticleHPVector.hh"
#include "globals.hh"

#include <CLHEP/Units/SystemOfUnits.h>

class G4ParticleHPParticleYield
{
  public:

    G4ParticleHPParticleYield()
    {
      simpleMean = true;
      spontPrompt = true;
      hasPromptData = false;
      hasDelayedData = false;

      targetMass = 0.0;
      theSpontPrompt = 0.0;
      spontDelayed = true;
      theSpontDelayed = 0.0;
    }

    ~G4ParticleHPParticleYield() = default;

    inline G4double GetTargetMass() const { return targetMass; }

    inline void InitMean(std::istream& aDataFile)
    {
      G4int iflag;
      aDataFile >> targetMass >> iflag;
      if (iflag == 1) simpleMean = false;
      if (simpleMean) {
        theSimpleMean.Init(aDataFile, CLHEP::eV);
      }
      else {
        theMean.Init(aDataFile);
      }
    }

    inline void InitPrompt(std::istream& aDataFile)
    {
      hasPromptData = true;
      G4int iflag;
      aDataFile >> targetMass >> iflag;
      if (iflag == 2) spontPrompt = false;
      if (spontPrompt) {
        aDataFile >> theSpontPrompt;
      }
      else {
        thePrompt.Init(aDataFile, CLHEP::eV);
      }
    }

    inline void InitDelayed(std::istream& aDataFile)
    {
      hasDelayedData = true;
      G4int iflag;
      aDataFile >> targetMass >> iflag;
      thePrecursorDecayConstants.Init(aDataFile, 1. / CLHEP::s);  // s is the CLHEP unit second
      if (iflag == 2) spontDelayed = false;
      if (spontDelayed) {
        aDataFile >> theSpontDelayed;
      }
      else {
        theDelayed.Init(aDataFile, CLHEP::eV);
      }
    }

    inline G4double GetMean(G4double anEnergy) const
    {
      if (simpleMean) {
        return theSimpleMean.GetY(anEnergy);
      }
      return theMean.GetValue(anEnergy);
    }

    inline G4double GetPrompt(G4double anEnergy) const
    {
      if (!hasPromptData) return 0;
      if (spontPrompt) {
        return theSpontPrompt;
      }
      return thePrompt.GetY(anEnergy);
    }

    inline G4double GetDelayed(G4double anEnergy) const
    {
      if (!hasDelayedData) return 0;
      if (spontDelayed) {
        return theSpontDelayed;
      }
      return theDelayed.GetY(anEnergy);
    }

    inline G4double GetDecayConstant(G4int i) const
    {
      return thePrecursorDecayConstants.GetValue(i);
    }

  private:

    G4double targetMass;
    // total mean
    G4bool simpleMean;
    G4ParticleHPPolynomExpansion theMean;
    G4ParticleHPVector theSimpleMean;

    // Prompt neutrons
    G4bool hasPromptData;
    G4bool spontPrompt;
    G4ParticleHPVector thePrompt;
    G4double theSpontPrompt;

    // delayed neutrons
    G4bool hasDelayedData;
    G4bool spontDelayed;
    G4ParticleHPList thePrecursorDecayConstants;
    G4ParticleHPVector theDelayed;
    G4double theSpontDelayed;
};

#endif
