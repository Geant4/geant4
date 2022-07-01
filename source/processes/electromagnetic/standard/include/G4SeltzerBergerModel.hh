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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4SeltzerBergerModel
//
// Author:        Andreas Schaelicke & Vladimir Ivantchenko
//
// Creation date: 04.10.2011
//
// Modifications:
//
// 24.07.2018 Introduced possibility to use sampling tables to sample the
//            emitted photon energy (instead of using rejectio) from the 
//            Seltzer-Berger scalled DCS for bremsstrahlung photon emission. 
//            Using these sampling tables option gives faster(30-70%) final 
//            state generation than the original rejection but takes some 
//            extra memory (+ ~6MB in the case of the full CMS detector). 
//            (M Novak)
//
// Class Description:
//
// Implementation of the bremssrahlung energy spectrum using
// 1. S.M. Seltzer and M.J. Berger Nucl. Instr. Meth. B12 (1985) 95
// 2. S.M. Seltzer and M.J. Berger Atomic data and Nuclear Data
//    Tables 35 (1986) 345
// Cross section computation in the base class G4eBremsstrahlungRelModel

// -------------------------------------------------------------------
//

#ifndef G4SeltzerBergerModel_h
#define G4SeltzerBergerModel_h 1

#include "G4eBremsstrahlungRelModel.hh"
#include "globals.hh"

class G4Physics2DVector;
class G4SBBremTable;

class G4SeltzerBergerModel : public G4eBremsstrahlungRelModel
{

public:

  explicit G4SeltzerBergerModel(const G4ParticleDefinition* p = nullptr,
                                const G4String& nam = "eBremSB");

  ~G4SeltzerBergerModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*,
                         const G4DynamicParticle*,
			 G4double cutEnergy,
                         G4double maxEnergy) override;

  void SetupForMaterial(const G4ParticleDefinition*,
                        const G4Material*, G4double) override;

  inline void SetBicubicInterpolationFlag(G4bool val) 
  { fIsUseBicubicInterpolation = val; };

  // hide assignment operator and cctr
  G4SeltzerBergerModel & operator=(const G4SeltzerBergerModel &right) = delete;
  G4SeltzerBergerModel(const G4SeltzerBergerModel&) = delete;

protected:

  G4double ComputeDXSectionPerAtom(G4double gammaEnergy) override;

private:

  void ReadData(G4int Z);

  const G4String& FindDirectoryPath();

  G4double SampleEnergyTransfer(const G4double kineticEnergy, 
                                const G4double logKineticEnergy, 
                                const G4double cut,
                                const G4double emax);

  static constexpr G4int    gMaxZet       = 101;
  static constexpr G4double gExpNumLimit  = -12.;
  static G4double           gYLimitData[gMaxZet];
  static G4Physics2DVector* gSBDCSData[gMaxZet];
  static G4SBBremTable*     gSBSamplingTable;
  static G4String           gDataDirectory;

  G4bool                    fIsUseBicubicInterpolation;
  G4bool                    fIsUseSamplingTables;

  G4int                     fNumWarnings;

  size_t                    fIndx;
  size_t                    fIndy;
};

#endif
