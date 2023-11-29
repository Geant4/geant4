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
// G4VRangeToEnergyConverter
//
// Class description:
//
// Base class for Range to Energy Converters.
// Cut in energy corresponding to given cut value in range
// is calculated for a material by using Convert() method.

// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------
#ifndef G4VRangeToEnergyConverter_hh
#define G4VRangeToEnergyConverter_hh 1

#include <vector>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"

class G4VRangeToEnergyConverter
{
public:

  explicit G4VRangeToEnergyConverter();

  virtual ~G4VRangeToEnergyConverter();

  // operators are not used
  G4VRangeToEnergyConverter(const G4VRangeToEnergyConverter& r) = delete;
  G4VRangeToEnergyConverter& operator=
  (const G4VRangeToEnergyConverter &r) = delete;
  G4bool operator==(const G4VRangeToEnergyConverter& r) const = delete;
  G4bool operator!=(const G4VRangeToEnergyConverter& r) const = delete;

  // Calculate energy cut from given range cut for the material
  virtual G4double Convert(const G4double rangeCut, const G4Material* material);

  // Set energy range for all particle type
  // if highedge > 10 GeV, highedge value is not changed
  static void SetEnergyRange(const G4double lowedge, const G4double highedge);

  // Get energy range for all particle type
  static G4double GetLowEdgeEnergy();
  static G4double GetHighEdgeEnergy();

  // Get/set max cut energy for all particle type
  // No check on the value
  static G4double GetMaxEnergyCut();
  static void SetMaxEnergyCut(const G4double value);
  
  // Return pointer to the particle type which this converter takes care of
  inline const G4ParticleDefinition* GetParticleType() const;

  inline void SetVerboseLevel(G4int value);
  inline G4int GetVerboseLevel() const;
      // control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

protected:

  virtual G4double ComputeValue(const G4int Z, const G4double kinEnergy) = 0;

private:

  static void FillEnergyVector(const G4double emin, const G4double emax);

  G4double ConvertForGamma(const G4double rangeCut, const G4Material* material);

  G4double ConvertForElectron(const G4double rangeCut, 
                              const G4Material* material);

  inline G4double LiniearInterpolation(const G4double e1, const G4double e2, 
                                       const G4double r1, const G4double r2, 
                                       const G4double r);

protected:

  const G4ParticleDefinition* theParticle = nullptr;
  G4int fPDG = 0;

private:

  static G4double sEmin;
  static G4double sEmax; 
  static std::vector<G4double>* sEnergy;
  static G4int sNbinPerDecade;
  static G4int sNbin;

  G4int verboseLevel = 1;
  G4bool isFirstInstance = false;
};

// ------------------
// Inline methods
// ------------------

inline 
void G4VRangeToEnergyConverter::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline 
G4int G4VRangeToEnergyConverter::GetVerboseLevel() const
{
  return verboseLevel;
}

inline 
const G4ParticleDefinition* G4VRangeToEnergyConverter::GetParticleType() const
{
  return theParticle;
}

inline G4double G4VRangeToEnergyConverter::LiniearInterpolation(
                const G4double e1, const G4double e2, 
                const G4double r1, const G4double r2, const G4double r)
{
  return (r1 == r2) ? e1 : e1 + (e2 - e1)*(r - r1)/(r2 - r1);
}

#endif
