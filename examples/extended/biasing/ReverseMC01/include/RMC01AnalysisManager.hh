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
/// \file biasing/ReverseMC01/include/RMC01AnalysisManager.hh
/// \brief Definition of the RMC01AnalysisManager class
//
// $Id$
//
//////////////////////////////////////////////////////////////
//  Class Name:        RMC01AnalysisManager
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//                 17-11-2009 creation by L. Desorgher
//                24-11-2009 L.Desorgher,
//                   -registering in  Conv* ASCII files every 5000 events the computed
//                        edep with  precision.
//                   -Correction of the adjoint computed current and answer matrices
//                        by a factor n_asked/n_processed for the case where a run is aborted
//                        because the user expected precision on e_dep has been reached.
//
//-------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RMC01AnalysisManager_HH
#define RMC01AnalysisManager_HH

#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"G4Run.hh"

class Histo1DVar;
class Histo2DVar;
class G4Timer;
class RMC01AnalysisManagerMessenger;

enum PRIM_SPECTRUM_TYPE{EXPO,POWER}; 

class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RMC01AnalysisManager
{
public:
  
  ~RMC01AnalysisManager();
  static RMC01AnalysisManager* GetInstance();
  
  public:
   
   void BeginOfRun(const G4Run*); 
   void EndOfRun(const G4Run*); 
   void BeginOfEvent(const G4Event*); 
   void EndOfEvent(const G4Event*); 
   
   void SetPrimaryExpSpectrumForAdjointSim(
                                       const G4String& particle_name, G4double fluence,
                                              G4double E0, G4double Emin,G4double Emax);
   void SetPrimaryPowerLawSpectrumForAdjointSim(
                                      const  G4String& particle_name,G4double fluence,
                                            G4double alpha, G4double Emin,G4double Emax);
   //precision of the simulation results is given in % by the user
   inline void SetPrecision(G4double precision){
                                         fPrecision_to_reach =precision/100.;};

private:

  static RMC01AnalysisManager* fInstance;

private:

  RMC01AnalysisManager(); 
private:
  
  void EndOfEventForForwardSimulation(const G4Event* anEvent);
  void EndOfEventForAdjointSimulation(const G4Event* anEvent);
  G4double PrimDiffAndDirFluxForAdjointSim(G4double prim_energy);
  void WriteHisto(Histo1DVar* anHisto, G4double scaling_factor,
                                              G4String fileName, G4String header_lines);
  void WriteHisto(Histo2DVar* anHisto, G4double scaling_factor,
                                              G4String fileName, G4String header_lines);
  void ResetHistograms();
  void ComputeMeanEdepAndError(const G4Event* anEvent,
                                                        G4double& mean,G4double& error);
  
private:
  
RMC01AnalysisManagerMessenger*  fMsg;
  
  //Histos for  fwd simulation
  //--------------
  Histo1DVar* fEdep_vs_prim_ekin;
  Histo1DVar* fElectron_current;
  Histo1DVar* fProton_current;
  Histo1DVar* fGamma_current;
  
  //Fluence
  //------------
  G4double fOmni_fluence_for_fwd_sim;
  
  //Variable to check the convergence of the energy deposited
  //            for forward and adjoint simulations
  //---------------------------------------------------------
  G4double fAccumulated_edep;
  G4double fAccumulated_edep2;
  G4double fMean_edep;
  G4double fError_mean_edep;
  G4double fRelative_error;
  G4double fElapsed_time;
  G4double fPrecision_to_reach;
  G4bool fStop_run_if_precision_reached;
  G4int fNb_evt_modulo_for_convergence_test;
  
  
  //Histos for adjoint simulation
  //-----------------------------
  Histo1DVar* fEdep_rmatrix_vs_electron_prim_energy;
  Histo2DVar* fElectron_current_rmatrix_vs_electron_prim_energy;
  Histo2DVar* fGamma_current_rmatrix_vs_electron_prim_energy;
  
  Histo1DVar* fEdep_rmatrix_vs_gamma_prim_energy;
  Histo2DVar* fElectron_current_rmatrix_vs_gamma_prim_energy;
  Histo2DVar* fGamma_current_rmatrix_vs_gamma_prim_energy;
  
  Histo1DVar* fEdep_rmatrix_vs_proton_prim_energy;
  Histo2DVar* fElectron_current_rmatrix_vs_proton_prim_energy;
  Histo2DVar* fProton_current_rmatrix_vs_proton_prim_energy;
  Histo2DVar* fGamma_current_rmatrix_vs_proton_prim_energy;
  
  
  //Prim spectrum to which the adjoint simulation will be normalised
  //Answer matrices will be also registered for post processing normalisation
  //--------------------------------------------------------
  PRIM_SPECTRUM_TYPE fPrimSpectrumType;
  G4int fPrimPDG_ID;
  G4double fAlpha_or_E0;
  G4double fAmplitude_prim_spectrum;
  G4double fEmin_prim_spectrum;
  G4double fEmax_prim_spectrum;
  G4bool fAdjoint_sim_mode;
  G4int fNb_evt_per_adj_evt;

  //Timer
  //------
  G4Timer* fTimer;
  std::fstream fConvergenceFileOutput;
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
