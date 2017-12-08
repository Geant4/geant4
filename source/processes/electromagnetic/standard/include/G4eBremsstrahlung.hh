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
// $Id: G4eBremsstrahlung.hh 106717 2017-10-20 09:41:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eBremsstrahlung
//
// Author:        Laszlo Urban
//
// Creation date: 24.06.1996
//
// Modifications:
//
// 01-10-96 new type G4OrderedTable;  ComputePartialSumSigma()
// 20-03-97 new energy loss+ionisation+brems scheme, L.Urban
// 01-09-98 new method  PrintInfo()
// 10-02-00 modifications , new e.m. structure, L.Urban
// 07-08-00 new cross section/en.loss parametrisation, LPM flag , L.Urban
// 09-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 19-09-01 come back to previous process name "eBrem"
// 29-10-01 all static functions no more inlined (mma)
// 07-01-02 new design of em processes (V.Ivanchenko)
// 26-12-02 secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 17-10-03 PrintInfoDefinition - virtual (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 21-01-04 Migrade to G4ParticleChangeForLoss (V.Ivanchenko)
// 04-11-04 add gamma threshold (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 11-04-04 Move MaxSecondaryEnergy to models (V.Ivanchenko)
//
// Class Description:
//
// This class manages the bremsstrahlung for e-/e+
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLoss.
//

// -------------------------------------------------------------------
//

#ifndef G4eBremsstrahlung_h
#define G4eBremsstrahlung_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;

class G4eBremsstrahlung : public G4VEnergyLossProcess
{

public:

  explicit G4eBremsstrahlung(const G4String& name = "eBrem");

  virtual ~G4eBremsstrahlung();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) final;
  
  // print documentation in html format
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  // Print out of the class parameters
  virtual void StreamProcessInfo(std::ostream& outFile,
                             G4String endOfLine=G4String("\n")) const override;

  virtual void 
  InitialiseEnergyLossProcess(const G4ParticleDefinition*,
			      const G4ParticleDefinition*) override;

  G4bool   isInitialised;

private:

  // hide assignment operator
  G4eBremsstrahlung & operator=(const G4eBremsstrahlung &right) = delete;
  G4eBremsstrahlung(const G4eBremsstrahlung&) = delete;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

