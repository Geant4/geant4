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
/// \file electromagnetic/TestEm16/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
: G4Run(),
  f_n_gam_sync(0),
  f_e_gam_sync(0.),
  f_e_gam_sync2(0.),
  f_e_gam_sync_max(0.),
  f_lam_gam_sync(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // accumulate sums
  //
  f_n_gam_sync += localRun->f_n_gam_sync;
  f_e_gam_sync += localRun->f_e_gam_sync;
  f_e_gam_sync2 += localRun->f_e_gam_sync2;
  f_e_gam_sync_max = std::max(f_e_gam_sync_max, localRun->f_e_gam_sync_max);
  f_lam_gam_sync += localRun->f_lam_gam_sync;
 
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
    if (f_n_gam_sync > 0)
    {
      G4double Emean = f_e_gam_sync/f_n_gam_sync;
      G4double E_rms = std::sqrt(f_e_gam_sync2/f_n_gam_sync - Emean*Emean);
      G4cout
      << "Summary for synchrotron radiation :" << '\n' << std::setprecision(4)
      << "  Number of photons = " << f_n_gam_sync << '\n'
      << "  Emean             = " << Emean/keV << " +/- "
      << E_rms/(keV * std::sqrt((G4double) f_n_gam_sync)) << " keV" << '\n'
      << "  E_rms             = " << G4BestUnit(E_rms,"Energy") << '\n'
      << "  Energy Max / Mean = " << f_e_gam_sync_max / Emean << '\n'
      << "  MeanFreePath      = " << G4BestUnit(f_lam_gam_sync/f_n_gam_sync,
                                            "Length")
      << G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
