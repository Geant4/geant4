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
/// \file biasing/ReverseMC01/src/RMC01AnalysisManager.cc
/// \brief Implementation of the RMC01AnalysisManager class
//
// $Id$
//
//////////////////////////////////////////////////////////////
//      Class Name:        RMC01AnalysisManager
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RMC01AnalysisManager.hh"
#include "G4AdjointSimManager.hh"
#include "G4SDManager.hh"
#include "RMC01SD.hh"
#include "G4THitsCollection.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "Histo1DVar.hh"
#include "Histo2DVar.hh"
#include "G4Timer.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "RMC01AnalysisManagerMessenger.hh"

RMC01AnalysisManager* RMC01AnalysisManager::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManager::RMC01AnalysisManager()
 :fOmni_fluence_for_fwd_sim(1./cm2),
  fAccumulated_edep(0.), fAccumulated_edep2(0.), fMean_edep(0.),
  fError_mean_edep(0.), fRelative_error(0.), fElapsed_time(0.),
  fPrecision_to_reach(0.),fStop_run_if_precision_reached(true),
  fNb_evt_modulo_for_convergence_test(5000),
  fPrimSpectrumType(EXPO),
  fAlpha_or_E0(.5*MeV),fAmplitude_prim_spectrum (1.),
  fEmin_prim_spectrum(1.*keV),fEmax_prim_spectrum (20.*MeV),
  fAdjoint_sim_mode(true),fNb_evt_per_adj_evt(2)
{ 
  
  fMsg = new RMC01AnalysisManagerMessenger(this);
  
  //----------------------
  //Creation of histograms
  //----------------------
  
  //Energy binning of the histograms : 60 log bins over [1keV-1GeV]
  
  G4double bins[61];
  size_t nbin=60;
  G4double emin=1.*keV;
  G4double emax=1.*GeV;
  for ( size_t i=0; i <= nbin; i++) {
           G4double val_bin;
        val_bin=emin * std::pow(10., i * std::log10(emax/emin)/nbin);
        G4double exp_10=4.-int(std::log10(val_bin));
        G4double factor =std::pow(10., exp_10);
        val_bin=int(factor*val_bin)/factor;
        bins[i] = val_bin;
  
  }
  
  //Histograms for :
  //        1)the forward simulation results 
  //        2)the Reverse MC  simulation results normalised to a user spectrum
  //------------------------------------------------------------------------
   
  fEdep_vs_prim_ekin = new Histo1DVar("Edep_vs_prim_ekin",bins,  nbin+1, LEFT);
  fElectron_current  = new Histo1DVar("Electron_current",bins,  nbin+1, LEFT);
  fProton_current= new Histo1DVar("Proton_current",bins,  nbin+1, LEFT);
  fGamma_current= new Histo1DVar("Gamma_current",bins,  nbin+1, LEFT);
  
  //Response matrices for the adjoint simulation only
  //-----------------------------------------------
  
  //Response matrices for external isotropic e- source
  
  fEdep_rmatrix_vs_electron_prim_energy =
     new Histo1DVar("Edep_rmatrix_vs_electron_prim_energy",
                                                                                                                            bins, nbin+1, LEFT);
          
  fElectron_current_rmatrix_vs_electron_prim_energy =
         new Histo2DVar("Electron_current_rmatrix_vs_electron_prim_energy",
                                                bins, nbin+1, LEFT, bins, nbin+1, LEFT);
          
  fGamma_current_rmatrix_vs_electron_prim_energy =
         new Histo2DVar("Gamma_current_rmatrix_vs_electron_prim_energy",
                                            bins, nbin+1, LEFT, bins, nbin+1, LEFT);
          
  
  //Response matrices for external isotropic gamma source
  
  fEdep_rmatrix_vs_gamma_prim_energy =
     new Histo1DVar("Edep_rmatrix_vs_gamma_prim_energy",
                                       bins, nbin+1, LEFT);
          
  fElectron_current_rmatrix_vs_gamma_prim_energy =
     new Histo2DVar("Electron_current_rmatrix_vs_gamma_prim_energy",
                                                                        bins, nbin+1, LEFT, bins, nbin+1, LEFT);
          
  fGamma_current_rmatrix_vs_gamma_prim_energy =
     new Histo2DVar("Gamma_current_rmatrix_vs_gamma_prim_energy",
                                                                        bins, nbin+1, LEFT, bins, nbin+1, LEFT);
          
  //Response matrices for external isotropic proton source
  
  fEdep_rmatrix_vs_proton_prim_energy =
     new Histo1DVar("Edep_rmatrix_vs_proton_prim_energy",
                                                        bins, nbin+1, LEFT);
          
  fElectron_current_rmatrix_vs_proton_prim_energy =
          new Histo2DVar("Electron_current_rmatrix_vs_proton_prim_energy",
                                                  bins, nbin+1, LEFT, bins, nbin+1, LEFT);
          
  fGamma_current_rmatrix_vs_proton_prim_energy =
          new Histo2DVar("Gamma_current_rmatrix_vs_proton_prim_energy",
                                                  bins, nbin+1, LEFT, bins, nbin+1, LEFT);
          
  fProton_current_rmatrix_vs_proton_prim_energy =
          new Histo2DVar("Proton_current_rmatrix_vs_proton_prim_energy",
                                                  bins, nbin+1, LEFT, bins, nbin+1, LEFT);
  
  //-------------
  //Timer for convergence vector
  //-------------
  
  fTimer = new G4Timer();

  //---------------------------------
  //Primary particle ID for normalisation of adjoint results
  //---------------------------------
  
  fPrimPDG_ID = G4Electron::Electron()->GetPDGEncoding();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManager::~RMC01AnalysisManager() 
{ 
  delete fEdep_vs_prim_ekin;
  delete fElectron_current;
  delete fProton_current;
  delete fGamma_current;
  
  delete fEdep_rmatrix_vs_electron_prim_energy;
  delete fElectron_current_rmatrix_vs_electron_prim_energy;
  delete fGamma_current_rmatrix_vs_electron_prim_energy;
  
  delete fEdep_rmatrix_vs_gamma_prim_energy;
  delete fElectron_current_rmatrix_vs_gamma_prim_energy;
  delete fGamma_current_rmatrix_vs_gamma_prim_energy;
  
  delete fEdep_rmatrix_vs_proton_prim_energy;
  delete fElectron_current_rmatrix_vs_proton_prim_energy;
  delete fProton_current_rmatrix_vs_proton_prim_energy;
  delete fGamma_current_rmatrix_vs_proton_prim_energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManager* RMC01AnalysisManager::GetInstance()
{
  if (fInstance == 0) fInstance = new RMC01AnalysisManager;
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::BeginOfRun(const G4Run* aRun)
{  fAccumulated_edep =0.;
   fAccumulated_edep2 =0.;
   fRelative_error=1.;
   fMean_edep=0.;
   fError_mean_edep=0.;
   fAdjoint_sim_mode =G4AdjointSimManager::GetInstance()->GetAdjointSimMode();

   if (fAdjoint_sim_mode){
           fNb_evt_per_adj_evt=aRun->GetNumberOfEventToBeProcessed()/
                                                   G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
        fConvergenceFileOutput.open("ConvergenceOfAdjointSimulationResults.txt",
                                                                             std::ios::out);
          fConvergenceFileOutput<<"Normalised Edep[MeV]\terror[MeV]\tcomputing_time[s]"
                                                                                  <<std::endl;
   }
   else {
           fConvergenceFileOutput.open("ConvergenceOfForwardSimulationResults.txt",
                                                                                std::ios::out);
          fConvergenceFileOutput<<"Edep per event [MeV]\terror[MeV]\tcomputing_time[s]"
                                                                                  <<std::endl;
   }
   fConvergenceFileOutput.setf(std::ios::scientific);
   fConvergenceFileOutput.precision(6);         
   ResetHistograms();
   fTimer->Start();
   fElapsed_time=0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::EndOfRun(const G4Run* aRun)
{ fTimer->Stop();
 
  if (!fAdjoint_sim_mode){
    G4cout<<"Results of forward simulation!"<<std::endl;
        G4cout<<"edep per event [MeV] = "<<fMean_edep<<std::endl;
        G4cout<<"error[MeV] = "<<fError_mean_edep<<std::endl;
        G4int nb_evt=aRun->GetNumberOfEvent();
         WriteHisto(fEdep_vs_prim_ekin,1./nb_evt,G4String("Fwd_Edep_vs_EkinPrim.txt"),
                               G4String("E1[MeV]\t\tE2[MeV]\t\tEdep[MeV]\terr_Edep[MeV]\n"));
         WriteHisto(fElectron_current,1./nb_evt, G4String("Fwd_ElectronCurrent.txt"),
                              G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering electron\terr\n"));
         WriteHisto(fProton_current,1./nb_evt, G4String("Fwd_ProtonCurrent.txt"),
                              G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering proton\terr[]\n"));
          WriteHisto(fGamma_current,1./nb_evt, G4String("Fwd_GammaCurrent.txt"),
                                G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering gamma\terr[]\n"));
  }
  
  else {
          G4cout<<"Results of reverse/adjoint simulation!"<<std::endl;
        G4cout<<"normalised edep [MeV] = "<<fMean_edep<<std::endl;
        G4cout<<"error[MeV] = "<<fError_mean_edep<<std::endl;

        
        G4double factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun()
                                                              *fNb_evt_per_adj_evt/aRun->GetNumberOfEvent();
                
        WriteHisto(fEdep_vs_prim_ekin,factor, G4String("Adj_Edep_vs_EkinPrim.txt"),
                              G4String("E1[MeV]\t\tE2[MeV]\t\tEdep[MeV]\terr_Edep[MeV]\n"));
         WriteHisto(fElectron_current,factor, G4String("Adj_ElectronCurrent.txt"),
                              G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering electron\terr\n"));
         WriteHisto(fProton_current,factor, G4String("Adj_ProtonCurrent.txt"),
                              G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering proton\terr[]\n"));
          WriteHisto(fGamma_current,factor, G4String("Adj_GammaCurrent.txt"),
                                G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering gamma\terr[]\n"));
  
          
        WriteHisto(fEdep_rmatrix_vs_electron_prim_energy,factor,
                                           G4String("Adj_Edep_vs_EkinPrimElectron_Answer.txt"),
                                           G4String("E1[MeV]\t\tE2[MeV]\t\tEdep Efficiency[MeV*cm2*MeV*str]\terr_Edep[MeV*cm2*MeV*str]\n"));
         
          WriteHisto(fElectron_current_rmatrix_vs_electron_prim_energy,factor,
                                           G4String("Adj_ElectronCurrent_vs_EkinPrimElectron_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
         
        WriteHisto(fGamma_current_rmatrix_vs_electron_prim_energy,factor,
                                           G4String("Adj_GammaCurrent_vs_EkinPrimElectron_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
        
        WriteHisto(fEdep_rmatrix_vs_gamma_prim_energy,factor,
                                           G4String("Adj_Edep_vs_EkinPrimGamma_Answer.txt"),
                                           G4String("E1[MeV]\t\tE2[MeV]\t\tEdep Efficiency[MeV*cm2*MeV*str]\terr_Edep[MeV*cm2*MeV*str]\n"));
         
          WriteHisto(fElectron_current_rmatrix_vs_gamma_prim_energy,factor,
                                           G4String("Adj_ElectronCurrent_vs_EkinPrimGamma_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
         
        WriteHisto(fGamma_current_rmatrix_vs_gamma_prim_energy,factor,
                                           G4String("Adj_GammaCurrent_vs_EkinPrimGamma_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
         
        
        
        WriteHisto(fEdep_rmatrix_vs_proton_prim_energy,factor,
                                           G4String("Adj_Edep_vs_EkinPrimProton_Answer.txt"),
                                           G4String("E1[MeV]\t\tE2[MeV]\t\tEdep Efficiency[MeV*cm2*MeV*str]\terr_Edep[MeV*cm2*MeV*str]\n"));
         
          WriteHisto(fElectron_current_rmatrix_vs_proton_prim_energy,factor,
                                           G4String("Adj_ElectronCurrent_vs_EkinPrimProton_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
         
        WriteHisto(fGamma_current_rmatrix_vs_proton_prim_energy,factor,
                                           G4String("Adj_GammaCurrent_vs_EkinPrimProton_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
         
        WriteHisto(fProton_current_rmatrix_vs_proton_prim_energy,factor,
                                           G4String("Adj_ProtonCurrent_vs_EkinPrimProton_Answer.txt"),
                                           G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
  }
  fConvergenceFileOutput.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::BeginOfEvent(const G4Event* )
{;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::EndOfEvent(const G4Event* anEvent)
{  

   if (fAdjoint_sim_mode) EndOfEventForAdjointSimulation(anEvent);
   else EndOfEventForForwardSimulation(anEvent);

   //Test convergence. The error is already computed
   //--------------------------------------
   G4int nb_event=anEvent->GetEventID()+1;
   //G4double factor=1.;
   if (fAdjoint_sim_mode) {
           G4double  n_adj_evt= nb_event/fNb_evt_per_adj_evt;
        // nb_event/fNb_evt_per_adj_evt;
        if (n_adj_evt*fNb_evt_per_adj_evt == nb_event) {
                nb_event =static_cast<G4int>(n_adj_evt);
                //factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
        }        
        else nb_event=0;
        
   }
   
   if (nb_event>100 && fStop_run_if_precision_reached &&
                                                fPrecision_to_reach >fRelative_error) {
                G4cout<<fPrecision_to_reach*100.<<"%  Precision reached!"<<std::endl;
                fTimer->Stop();
                fElapsed_time+=fTimer->GetRealElapsed();
                fConvergenceFileOutput<<fMean_edep<<'\t'<<fError_mean_edep
                                                               <<'\t'<<fElapsed_time<<std::endl;
                G4RunManager::GetRunManager()->AbortRun(true);
   }
   
   if (nb_event>0 && nb_event % fNb_evt_modulo_for_convergence_test == 0) {
           fTimer->Stop();
        fElapsed_time+=fTimer->GetRealElapsed();
        fTimer->Start();
        fConvergenceFileOutput<<fMean_edep<<'\t'<<fError_mean_edep<<'\t'
                                                                 <<fElapsed_time<<std::endl;
   }
}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RMC01AnalysisManager::EndOfEventForForwardSimulation(
                                                                const G4Event* anEvent)
{  
   
   G4SDManager* SDman = G4SDManager::GetSDMpointer();
   G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
   RMC01DoubleWithWeightHitsCollection* edepCollection =
                   (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("edep")));

   RMC01DoubleWithWeightHitsCollection* electronCurrentCollection =
                   (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_electron")));

   RMC01DoubleWithWeightHitsCollection* protonCurrentCollection =
                      (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_proton")));

   RMC01DoubleWithWeightHitsCollection* gammaCurrentCollection =
                      (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_gamma")));
   
   //Total energy deposited in Event
   //-------------------------------
   G4double totEdep=0; 
   G4int i;
   for (i=0;i<edepCollection->entries();i++)
        totEdep+=(*edepCollection)[i]->GetValue()*(*edepCollection)[i]->GetWeight();
   
   if (totEdep>0.){
           fAccumulated_edep +=totEdep ;
           fAccumulated_edep2 +=totEdep*totEdep;
    G4PrimaryParticle* thePrimary=anEvent->GetPrimaryVertex()->GetPrimary();
           G4double E0= thePrimary->GetG4code()->GetPDGMass();
           G4double P=thePrimary->GetMomentum().mag();
           G4double prim_ekin =std::sqrt(E0*E0+P*P)-E0;
           fEdep_vs_prim_ekin->fill(prim_ekin,totEdep);
   } 
   ComputeMeanEdepAndError(anEvent,fMean_edep,fError_mean_edep);
   if (fError_mean_edep>0) fRelative_error= fError_mean_edep/fMean_edep;
                   
   //Particle current on sensitive cylinder
   //-------------------------------------
   
   for (i=0;i<electronCurrentCollection->entries();i++) {
           G4double ekin =(*electronCurrentCollection)[i]->GetValue();
        G4double weight=(*electronCurrentCollection)[i]->GetWeight();
        fElectron_current->fill(ekin,weight);
   }
   
   for (i=0;i<protonCurrentCollection->entries();i++) {
           G4double ekin =(*protonCurrentCollection)[i]->GetValue();
        G4double weight=(*protonCurrentCollection)[i]->GetWeight();
        fProton_current->fill(ekin,weight);
   }        
   
   for (i=0;i<gammaCurrentCollection->entries();i++) {
           G4double ekin =(*gammaCurrentCollection)[i]->GetValue();
        G4double weight=(*gammaCurrentCollection)[i]->GetWeight();
        fGamma_current->fill(ekin,weight);
   }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RMC01AnalysisManager::EndOfEventForAdjointSimulation(
                                                                const G4Event* anEvent)
{  
  //Output from Sensitive volume computed during the forward tracking phase
  //-----------------------------------------------------------------------
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
  RMC01DoubleWithWeightHitsCollection* edepCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("edep")));

  RMC01DoubleWithWeightHitsCollection* electronCurrentCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_electron")));

  RMC01DoubleWithWeightHitsCollection* protonCurrentCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_proton")));

  RMC01DoubleWithWeightHitsCollection* gammaCurrentCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_gamma")));
  
  //Output from adjoint tracking phase
  //----------------------------------------------------------------------------
  
  G4AdjointSimManager* theAdjointSimManager = G4AdjointSimManager::GetInstance();
  G4int pdg_nb =theAdjointSimManager->GetFwdParticlePDGEncodingAtEndOfLastAdjointTrack();
  G4double prim_ekin=theAdjointSimManager->GetEkinAtEndOfLastAdjointTrack();
  G4double adj_weight=theAdjointSimManager->GetWeightAtEndOfLastAdjointTrack();
 
  
  //Factor of normalisation to user selected primary spectrum (power law or exponential) 
  //------------------------------------------------------------------------------------
  
  G4double normalised_weight = 0.;
  if (pdg_nb== fPrimPDG_ID && prim_ekin>= fEmin_prim_spectrum
                                                     && prim_ekin<= fEmax_prim_spectrum)
                  normalised_weight =
                                  adj_weight*PrimDiffAndDirFluxForAdjointSim(prim_ekin);
  
  //Answer matrices
  //-------------
   Histo1DVar* edep_rmatrix =0;
   Histo2DVar* electron_current_rmatrix =0;
   Histo2DVar* gamma_current_rmatrix =0;
   Histo2DVar* proton_current_rmatrix =0;

   if (pdg_nb == G4Electron::Electron()->GetPDGEncoding()){ //e- answer matrices
        edep_rmatrix = fEdep_rmatrix_vs_electron_prim_energy;

        electron_current_rmatrix = fElectron_current_rmatrix_vs_electron_prim_energy;

        gamma_current_rmatrix = fGamma_current_rmatrix_vs_electron_prim_energy;
   }
   else if (pdg_nb == G4Gamma::Gamma()->GetPDGEncoding()){ //gammma answer matrices

        edep_rmatrix = fEdep_rmatrix_vs_gamma_prim_energy;
        electron_current_rmatrix = fElectron_current_rmatrix_vs_gamma_prim_energy;
        gamma_current_rmatrix = fGamma_current_rmatrix_vs_gamma_prim_energy;
   }
   else if (pdg_nb == G4Proton::Proton()->GetPDGEncoding()){ //proton answer matrices
        edep_rmatrix = fEdep_rmatrix_vs_proton_prim_energy;
        electron_current_rmatrix = fElectron_current_rmatrix_vs_proton_prim_energy;
        gamma_current_rmatrix = fGamma_current_rmatrix_vs_proton_prim_energy;
        proton_current_rmatrix = fProton_current_rmatrix_vs_proton_prim_energy;
   }
  
  //Registering of total energy deposited in Event
  //-------------------------------
   G4double totEdep=0; 
   G4int i;
   for (i=0;i<edepCollection->entries();i++)
           totEdep+=(*edepCollection)[i]->GetValue()*
                                            (*edepCollection)[i]->GetWeight();
   
   G4bool new_mean_computed=false;
   if (totEdep>0.){
           if (normalised_weight>0.){
                G4double edep=totEdep* normalised_weight;
                
                //Check if the edep is not wrongly too high
                //----------------------------------------- 
                G4double new_mean , new_error;
                
                fAccumulated_edep +=edep;
                   fAccumulated_edep2 +=edep*edep;
                
                ComputeMeanEdepAndError(anEvent,new_mean,new_error);
                G4double new_relative_error = 1.;
                if ( new_error >0) new_relative_error = new_error/ new_mean;
                if (fRelative_error <0.10 && new_relative_error>1.5*fRelative_error) {
                        G4cout<<"Potential wrong adjoint weight!"<<std::endl;
                        G4cout<<"The results of this event will not be registered!"
                                                                                     <<std::endl;
                        G4cout<<"previous mean edep [MeV] "<< fMean_edep<<std::endl;
                        G4cout<<"previous relative error "<< fRelative_error<<std::endl;
                        G4cout<<"new rejected mean edep [MeV] "<< new_mean<<std::endl;
                        G4cout<<"new rejected relative error "<< new_relative_error
                                                                                      <<std::endl;
                        fAccumulated_edep -=edep;
                        fAccumulated_edep2 -=edep*edep;
                        return; 
                }
                else { //accepted
                        fMean_edep = new_mean;
                        fError_mean_edep = new_error;
                        fRelative_error =new_relative_error;
                        new_mean_computed=true;
                }         
                fEdep_vs_prim_ekin->fill(prim_ekin,edep); 
        }
        
        // Registering answer matrix
        //---------------------------
        
        edep_rmatrix->fill(prim_ekin,totEdep*adj_weight/cm2);
   }
   if (!new_mean_computed){
            ComputeMeanEdepAndError(anEvent,fMean_edep,fError_mean_edep);
            if (fError_mean_edep>0) fRelative_error= fError_mean_edep/fMean_edep;
   }         
    
   
  //Registering of current of particles on the sensitive volume
  //------------------------------------------------------------
   
   for (i=0;i<electronCurrentCollection->entries();i++) {
           G4double ekin =(*electronCurrentCollection)[i]->GetValue();
        G4double weight=(*electronCurrentCollection)[i]->GetWeight();
        fElectron_current->fill(ekin,weight*normalised_weight);
        electron_current_rmatrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
   }
   
   for (i=0;i<protonCurrentCollection->entries();i++) {
           G4double ekin =(*protonCurrentCollection)[i]->GetValue();
        G4double weight=(*protonCurrentCollection)[i]->GetWeight();
        fProton_current->fill(ekin,weight*normalised_weight);
        proton_current_rmatrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
  }        
   
   for (i=0;i<gammaCurrentCollection->entries();i++) {
           G4double ekin =(*gammaCurrentCollection)[i]->GetValue();
        G4double weight=(*gammaCurrentCollection)[i]->GetWeight();
        fGamma_current->fill(ekin,weight*normalised_weight);
        gamma_current_rmatrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);

   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RMC01AnalysisManager::PrimDiffAndDirFluxForAdjointSim(
                                                                  G4double prim_energy)
{ 
  G4double flux=fAmplitude_prim_spectrum;
  if ( fPrimSpectrumType ==EXPO)        flux*=std::exp(-prim_energy/fAlpha_or_E0);
  else flux*=std::pow(prim_energy, -fAlpha_or_E0);
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RMC01AnalysisManager::WriteHisto(Histo1DVar* anHisto,
                      G4double scaling_factor, G4String fileName, G4String header_lines)
{ std::fstream FileOutput(fileName, std::ios::out);
  FileOutput<<header_lines;
  FileOutput.setf(std::ios::scientific);
  FileOutput.precision(6);
  size_t nxbins = anHisto->part.total_bins();
  for (size_t i =0;i<nxbins;i++) {
        FileOutput<<anHisto->part.get_bin_position(i)
              <<'\t'<<anHisto->part.get_bin_position(i+1)
              <<'\t'<<anHisto->get_bin_value(int(i))*scaling_factor
              <<'\t'<<anHisto->get_bin_error(int(i))*scaling_factor<<std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RMC01AnalysisManager::WriteHisto(Histo2DVar* anHisto,
                     G4double scaling_factor, G4String fileName, G4String header_lines)
{ std::fstream FileOutput(fileName, std::ios::out);
  FileOutput<<header_lines;
  
  FileOutput.setf(std::ios::scientific);
  FileOutput.precision(6);
  
  size_t nxbins = anHisto->x_part.total_bins();
  size_t nybins = anHisto->y_part.total_bins();
  for (size_t i =0;i<nxbins;i++) {
    for (size_t j =0;j<nybins;j++) {
          FileOutput<<anHisto->x_part.get_bin_position(i)
                        <<'\t'<<anHisto->x_part.get_bin_position(i+1)
                        <<'\t'<<anHisto->y_part.get_bin_position(j)
                        <<'\t'<<anHisto->y_part.get_bin_position(j+1)
                        <<'\t'<<anHisto->get_bin_value(int(i),int(j))*scaling_factor
                      <<'\t'<<anHisto->get_bin_error(int(i),int(j))*scaling_factor
                                                                            <<std::endl;
           }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::ResetHistograms()
{ fEdep_vs_prim_ekin->reset();
  fElectron_current->reset();
  fProton_current->reset();
  fGamma_current->reset();
  
  fEdep_rmatrix_vs_electron_prim_energy->reset();
  fElectron_current_rmatrix_vs_electron_prim_energy->reset();
  fGamma_current_rmatrix_vs_electron_prim_energy->reset();
  
  fEdep_rmatrix_vs_gamma_prim_energy->reset();
  fElectron_current_rmatrix_vs_gamma_prim_energy->reset();
  fGamma_current_rmatrix_vs_gamma_prim_energy->reset();
  
  fEdep_rmatrix_vs_proton_prim_energy->reset();
  fElectron_current_rmatrix_vs_proton_prim_energy->reset();
  fProton_current_rmatrix_vs_proton_prim_energy->reset();
  fGamma_current_rmatrix_vs_proton_prim_energy->reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::ComputeMeanEdepAndError(
                                         const G4Event* anEvent,G4double& mean,G4double& error)
{  
   G4int nb_event=anEvent->GetEventID()+1;
   G4double factor=1.;
   if (fAdjoint_sim_mode) {
           nb_event /=fNb_evt_per_adj_evt;
        factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
   }
   
   //error computation
   if (nb_event>1) {
             mean = fAccumulated_edep/nb_event;
          G4double mean_x2 =fAccumulated_edep2/nb_event;
            error = factor*std::sqrt(mean_x2-mean*mean)/std::sqrt(G4double(nb_event));
          mean *=factor;
   } else {
          mean=0;
          error=0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::SetPrimaryExpSpectrumForAdjointSim(
                     const G4String& particle_name, G4double omni_fluence,
                                              G4double E0, G4double Emin,G4double Emax)
{ fPrimSpectrumType = EXPO;
  if (particle_name == "e-" ) fPrimPDG_ID =
                      G4Electron::Electron()->GetPDGEncoding();
  else if (particle_name == "gamma") fPrimPDG_ID =
                                       G4Gamma::Gamma()->GetPDGEncoding();
  else if (particle_name == "proton") fPrimPDG_ID =
                                      G4Proton::Proton()->GetPDGEncoding();
  else {
          G4cout<<"The particle that you did select is not in the candidate list for primary [e-, gamma, proton]!"<<G4endl;
          return;
  }        
  fAlpha_or_E0 = E0 ;
  fAmplitude_prim_spectrum = omni_fluence/E0/
                                         (std::exp(-Emin/E0)-std::exp(-Emax/E0))/4./pi;
  fEmin_prim_spectrum = Emin ;
  fEmax_prim_spectrum = Emax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::SetPrimaryPowerLawSpectrumForAdjointSim(
                                  const G4String& particle_name, G4double omni_fluence,
                                            G4double alpha, G4double Emin,G4double Emax)
{ fPrimSpectrumType  =POWER;
  if (particle_name == "e-" ) fPrimPDG_ID =
                              G4Electron::Electron()->GetPDGEncoding();
  else if (particle_name == "gamma") fPrimPDG_ID =
                              G4Gamma::Gamma()->GetPDGEncoding();
  else if (particle_name == "proton") fPrimPDG_ID =
                              G4Proton::Proton()->GetPDGEncoding();
  else {
          G4cout<<"The particle that you did select is not in the candidate list for primary [e-, gamma, proton]!"<<G4endl;
          return;
  }        
  

 if (alpha ==1.) {
         fAmplitude_prim_spectrum = omni_fluence/std::log(Emax/Emin)/4./pi;
 }
 else {
         G4double p=1.-alpha;
         fAmplitude_prim_spectrum = omni_fluence/p/(std::pow(Emax,p)
                                                   -std::pow(Emin,p))/4./pi;
 }

  fAlpha_or_E0 = alpha ;
  fEmin_prim_spectrum = Emin ;
  fEmax_prim_spectrum = Emax;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

