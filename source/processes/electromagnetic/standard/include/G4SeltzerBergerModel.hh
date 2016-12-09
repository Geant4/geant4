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
// $Id: G4SeltzerBergerModel.hh 98737 2016-08-09 12:51:38Z gcosmo $
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

class G4SeltzerBergerModel : public G4eBremsstrahlungRelModel
{

public:

  explicit G4SeltzerBergerModel(const G4ParticleDefinition* p = nullptr, 
				const G4String& nam = "eBremSB");

  virtual ~G4SeltzerBergerModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double cutEnergy,
				 G4double maxEnergy) override;

  inline void SetBicubicInterpolationFlag(G4bool);

protected:

  virtual G4double ComputeDXSectionPerAtom(G4double gammaEnergy) override;

  virtual G4String DirectoryPath() const;

private:

  void ReadData(G4int Z, const char* path = 0);

  // hide assignment operator
  G4SeltzerBergerModel & operator=(const  G4SeltzerBergerModel &right) = delete;
  G4SeltzerBergerModel(const  G4SeltzerBergerModel&) = delete;

  static G4Physics2DVector* dataSB[101];
  static G4double ylimit[101];
  static G4double expnumlim;
  G4int  nwarn;
  size_t idx;
  size_t idy;
  G4bool useBicubicInterpolation;
};

inline void G4SeltzerBergerModel::SetBicubicInterpolationFlag(G4bool val)
{
  useBicubicInterpolation = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif
