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
// G4VelocityTable
//
// Class description:
//
// This class keeps a table of velocity as a function of the ratio
// kinetic erngy and mass. G4VelocityTable is used by
// G4Track::CalculateVelocity().

// Author: Hisaya Kurashige,  17 August 2011
// --------------------------------------------------------------------
#ifndef G4VelocityTable_hh
#define G4VelocityTable_hh 1

#include <vector>
#include <iostream>

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreadLocalSingleton.hh"

class G4VelocityTable
{
  friend class G4ThreadLocalSingleton<G4VelocityTable>;

  using G4VTDataVector = std::vector<G4double>;

  public:

    G4double Value(G4double theEnergy);
      // Get the cross-section/energy-loss value corresponding to the
      // given energy. An appropriate interpolation is used to calculate
      // the value

    static G4VelocityTable* GetVelocityTable();

    static void SetVelocityTableProperties(G4double t_max, G4double t_min,
                                           G4int nbin);
    static G4double GetMaxTOfVelocityTable();
    static G4double GetMinTOfVelocityTable();
    static G4int GetNbinOfVelocityTable();

  private:

    G4VelocityTable();
    ~G4VelocityTable();

    void PrepareVelocityTable();

    std::size_t FindBinLocation(G4double theEnergy) const;
      // Find the bin# in which theEnergy belongs - pure virtual function

    inline G4double Interpolation() const;

    // --------------------------------------------------------------------

    G4double edgeMin = 0.0;  // Energy of first point
    G4double edgeMax = 0.0;  // Energy of the last point

    std::size_t numberOfNodes = 0;

    G4VTDataVector dataVector;     // Vector to keep the crossection/energyloss
    G4VTDataVector binVector;      // Vector to keep energy
    G4VTDataVector secDerivative;  // Vector to keep second derivatives

    G4double dBin = 0.0;     // Bin width - useful only for fixed binning
    G4double baseBin = 0.0;  // Set this in constructor for performance

    G4double lastEnergy = -DBL_MAX;  // Cache the last input value
    G4double lastValue = 0.0;        // Cache the last output value
    std::size_t lastBin = 0;         // Cache the last bin location

    static G4ThreadLocal G4VelocityTable* theInstance;

    G4double maxT = 1000.0;
    G4double minT = 0.0001;
    G4int NbinT = 500;
};

// ----------------------
// Inline methods
// ----------------------

inline G4double G4VelocityTable::Interpolation() const
{
  // Linear interpolation is used to get the value. If the give energy
  // is in the highest bin, no interpolation will be Done. Because
  // there is an extra bin hidden from a user at locBin=numberOfBin,
  // the following interpolation is valid even the current locBin=
  // numberOfBin-1.

  G4double intplFactor =
    (lastEnergy - binVector[lastBin]) /
    (binVector[lastBin + 1] - binVector[lastBin]);  // Interpol. factor

  return dataVector[lastBin] +
         (dataVector[lastBin + 1] - dataVector[lastBin]) * intplFactor;
}

#endif
