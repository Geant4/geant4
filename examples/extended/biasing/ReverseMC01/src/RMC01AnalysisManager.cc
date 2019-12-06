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
#include "G4Timer.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "RMC01AnalysisManagerMessenger.hh"

RMC01AnalysisManager* RMC01AnalysisManager::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManager::RMC01AnalysisManager()
 :fAccumulated_edep(0.), fAccumulated_edep2(0.), fMean_edep(0.),
  fError_mean_edep(0.), fRelative_error(0.), fElapsed_time(0.),
  fPrecision_to_reach(0.),fStop_run_if_precision_reached(true),
  fNb_evt_modulo_for_convergence_test(5000),
  fEdep_rmatrix_vs_electron_prim_energy(0),
  fElectron_current_rmatrix_vs_electron_prim_energy(0),
  fGamma_current_rmatrix_vs_electron_prim_energy(0),
  fEdep_rmatrix_vs_gamma_prim_energy(0),
  fElectron_current_rmatrix_vs_gamma_prim_energy(0),
  fGamma_current_rmatrix_vs_gamma_prim_energy(0),
  fEdep_rmatrix_vs_proton_prim_energy(0),
  fElectron_current_rmatrix_vs_proton_prim_energy(0),
  fProton_current_rmatrix_vs_proton_prim_energy(0),
  fGamma_current_rmatrix_vs_proton_prim_energy(0),
  fFactoryOn(false),
  fPrimSpectrumType(EXPO),
  fAlpha_or_E0(.5*MeV),fAmplitude_prim_spectrum (1.),
  fEmin_prim_spectrum(1.*keV),fEmax_prim_spectrum (20.*MeV),
  fAdjoint_sim_mode(true),fNb_evt_per_adj_evt(2)
{ 
  
  fMsg = new RMC01AnalysisManagerMessenger(this);

  //-------------
  //Timer for convergence vector
  //-------------
  
  fTimer = new G4Timer();

  //---------------------------------
  //Primary particle ID for normalisation of adjoint results
  //---------------------------------
  
  fPrimPDG_ID = G4Electron::Electron()->GetPDGEncoding();
  
  fFileName[0] = "sim";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManager::~RMC01AnalysisManager() 
{;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManager* RMC01AnalysisManager::GetInstance()
{
  if (fInstance == 0) fInstance = new RMC01AnalysisManager;
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::BeginOfRun(const G4Run* aRun)
{

   fAccumulated_edep =0.;
   fAccumulated_edep2 =0.;
   fNentry = 0.0;
   fRelative_error=1.;
   fMean_edep=0.;
   fError_mean_edep=0.;
   fAdjoint_sim_mode =G4AdjointSimManager::GetInstance()->GetAdjointSimMode();

   if (fAdjoint_sim_mode){
           fNb_evt_per_adj_evt=aRun->GetNumberOfEventToBeProcessed()/
                       G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
    fConvergenceFileOutput.open("ConvergenceOfAdjointSimulationResults.txt",
                                                                std::ios::out);
    fConvergenceFileOutput<<
           "Normalised Edep[MeV]\terror[MeV]\tcomputing_time[s]"<<std::endl;
   }
   else {
     fConvergenceFileOutput.open("ConvergenceOfForwardSimulationResults.txt",
                                                               std::ios::out);
     fConvergenceFileOutput<<
         "Edep per event [MeV]\terror[MeV]\tcomputing_time[s]"
                                                                 <<std::endl;
   }
   fConvergenceFileOutput.setf(std::ios::scientific);
   fConvergenceFileOutput.precision(6);         

   fTimer->Start();
   fElapsed_time=0.;

   Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::EndOfRun(const G4Run* aRun)
{ fTimer->Stop();
  G4int nb_evt=aRun->GetNumberOfEvent();
  G4double factor =1./ nb_evt;
  if (!fAdjoint_sim_mode){
   G4cout<<"Results of forward simulation!"<<std::endl;
   G4cout<<"edep per event [MeV] = "<<fMean_edep<<std::endl;
   G4cout<<"error[MeV] = "<<fError_mean_edep<<std::endl;
  }
  
  else {
   G4cout<<"Results of reverse/adjoint simulation!"<<std::endl;
   G4cout<<"normalised edep [MeV] = "<<fMean_edep<<std::endl;
   G4cout<<"error[MeV] = "<<fError_mean_edep<<std::endl;
   factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun()
                                 *fNb_evt_per_adj_evt/aRun->GetNumberOfEvent();
  }
  Save(factor);
  fConvergenceFileOutput.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::BeginOfEvent(const G4Event* )
{ ;
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
         (RMC01DoubleWithWeightHitsCollection*)
              (HCE->GetHC(SDman->GetCollectionID("edep")));

   RMC01DoubleWithWeightHitsCollection* electronCurrentCollection =
             (RMC01DoubleWithWeightHitsCollection*)
                (HCE->GetHC(SDman->GetCollectionID("current_electron")));

   RMC01DoubleWithWeightHitsCollection* protonCurrentCollection =
             (RMC01DoubleWithWeightHitsCollection*)
                   (HCE->GetHC(SDman->GetCollectionID("current_proton")));

   RMC01DoubleWithWeightHitsCollection* gammaCurrentCollection =
             (RMC01DoubleWithWeightHitsCollection*)
                     (HCE->GetHC(SDman->GetCollectionID("current_gamma")));
   
   //Total energy deposited in Event
   //-------------------------------
   G4double totEdep=0; 
   std::size_t i;
   for (i=0;i<edepCollection->entries();++i)
        totEdep+=(*edepCollection)[i]->GetValue()
                  *(*edepCollection)[i]->GetWeight();
   
   if (totEdep>0.){
           fAccumulated_edep +=totEdep ;
           fAccumulated_edep2 +=totEdep*totEdep;
           fNentry += 1.0;
           G4PrimaryParticle* thePrimary=
                                 anEvent->GetPrimaryVertex()->GetPrimary();
           G4double E0= thePrimary->GetG4code()->GetPDGMass();
           G4double P=thePrimary->GetMomentum().mag();
           G4double prim_ekin =std::sqrt(E0*E0+P*P)-E0;
           fEdep_vs_prim_ekin->fill(prim_ekin,totEdep);
   } 
   ComputeMeanEdepAndError(anEvent,fMean_edep,fError_mean_edep);
   if (fError_mean_edep>0) fRelative_error= fError_mean_edep/fMean_edep;
                   
   //Particle current on sensitive cylinder
   //-------------------------------------
   
   for (i=0;i<electronCurrentCollection->entries();++i) {
           G4double ekin =(*electronCurrentCollection)[i]->GetValue();
        G4double weight=(*electronCurrentCollection)[i]->GetWeight();
        fElectron_current->fill(ekin,weight);
   }
   
   for (i=0;i<protonCurrentCollection->entries();++i) {
           G4double ekin =(*protonCurrentCollection)[i]->GetValue();
        G4double weight=(*protonCurrentCollection)[i]->GetWeight();
        fProton_current->fill(ekin,weight);
   }        
   
   for (i=0;i<gammaCurrentCollection->entries();++i) {
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
                  (RMC01DoubleWithWeightHitsCollection*)(
                    HCE->GetHC(SDman->GetCollectionID("edep")));

  RMC01DoubleWithWeightHitsCollection* electronCurrentCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(
                 HCE->GetHC(SDman->GetCollectionID("current_electron")));

  RMC01DoubleWithWeightHitsCollection* protonCurrentCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(
                 HCE->GetHC(SDman->GetCollectionID("current_proton")));

  RMC01DoubleWithWeightHitsCollection* gammaCurrentCollection =
                  (RMC01DoubleWithWeightHitsCollection*)(
                 HCE->GetHC(SDman->GetCollectionID("current_gamma")));
  
  //Computation of total energy deposited in fwd tracking phase
  //-------------------------------
   G4double totEdep=0;
   std::size_t i;
   for (i=0;i<edepCollection->entries();++i)
       totEdep+=(*edepCollection)[i]->GetValue()*
                                            (*edepCollection)[i]->GetWeight();

  //Output from adjoint tracking phase
  //----------------------------------------------------------------------------
  
  G4AdjointSimManager* theAdjointSimManager =
                   G4AdjointSimManager::GetInstance();

  size_t nb_adj_track =
       theAdjointSimManager->GetNbOfAdointTracksReachingTheExternalSurface();
  G4double total_normalised_weight = 0.;

  //We need to loop over the adjoint tracks that have reached the external
  //surface.
  for (std::size_t j=0;j<nb_adj_track;++j) {
    G4int pdg_nb =theAdjointSimManager
         ->GetFwdParticlePDGEncodingAtEndOfLastAdjointTrack(j);
    G4double prim_ekin=theAdjointSimManager
                              ->GetEkinAtEndOfLastAdjointTrack(j);
    G4double adj_weight=theAdjointSimManager
                             ->GetWeightAtEndOfLastAdjointTrack(j);
 
  
    //Factor of normalisation to user defined prim spectrum (power law or exp)
    //------------------------------------------------------------------------
    G4double normalised_weight = 0.;
    if (pdg_nb== fPrimPDG_ID && prim_ekin>= fEmin_prim_spectrum
                                         && prim_ekin<= fEmax_prim_spectrum)
      normalised_weight =
                adj_weight*PrimDiffAndDirFluxForAdjointSim(prim_ekin);
    total_normalised_weight += normalised_weight;
  
    //Answer matrices
    //-------------
    G4H1* edep_rmatrix =0;
    G4H2* electron_current_rmatrix =0;
    G4H2* gamma_current_rmatrix =0;
    G4H2* proton_current_rmatrix =0;

    if (pdg_nb == G4Electron::Electron()->GetPDGEncoding()){ //e- matrices
      edep_rmatrix = fEdep_rmatrix_vs_electron_prim_energy;
      electron_current_rmatrix =
                          fElectron_current_rmatrix_vs_electron_prim_energy;
      gamma_current_rmatrix = fGamma_current_rmatrix_vs_electron_prim_energy;
    }
    else if (pdg_nb == G4Gamma::Gamma()->GetPDGEncoding()){
      //gammma answer matrices
      edep_rmatrix = fEdep_rmatrix_vs_gamma_prim_energy;
      electron_current_rmatrix = fElectron_current_rmatrix_vs_gamma_prim_energy;
      gamma_current_rmatrix = fGamma_current_rmatrix_vs_gamma_prim_energy;
    }
    else if (pdg_nb == G4Proton::Proton()->GetPDGEncoding()){
      //proton answer matrices
      edep_rmatrix = fEdep_rmatrix_vs_proton_prim_energy;
      electron_current_rmatrix =
                         fElectron_current_rmatrix_vs_proton_prim_energy;
      gamma_current_rmatrix = fGamma_current_rmatrix_vs_proton_prim_energy;
      proton_current_rmatrix = fProton_current_rmatrix_vs_proton_prim_energy;
    }
    //Register histo edep vs prim ekin
    //----------------------------------
    if (normalised_weight>0) fEdep_vs_prim_ekin
                        ->fill(prim_ekin,totEdep*normalised_weight);
    // Registering answer matrix
    //---------------------------
    edep_rmatrix->fill(prim_ekin,totEdep*adj_weight/cm2);
    
  //Registering of current of particles on the sensitive volume
  //------------------------------------------------------------
   
   for (i=0;i<electronCurrentCollection->entries();++i) {
     G4double ekin =(*electronCurrentCollection)[i]->GetValue();
     G4double weight=(*electronCurrentCollection)[i]->GetWeight();
     fElectron_current->fill(ekin,weight*normalised_weight);
     electron_current_rmatrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
   }
   for (i=0;i<protonCurrentCollection->entries();++i) {
     G4double ekin =(*protonCurrentCollection)[i]->GetValue();
     G4double weight=(*protonCurrentCollection)[i]->GetWeight();
     fProton_current->fill(ekin,weight*normalised_weight);
     proton_current_rmatrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
   }
   for (i=0;i<gammaCurrentCollection->entries();++i) {
     G4double ekin =(*gammaCurrentCollection)[i]->GetValue();
     G4double weight=(*gammaCurrentCollection)[i]->GetWeight();
     fGamma_current->fill(ekin,weight*normalised_weight);
     gamma_current_rmatrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
   }
  }

  //Registering of total energy deposited in Event
  //-------------------------------
     G4bool new_mean_computed=false;
     if (totEdep>0.){
       if (total_normalised_weight>0.){
         G4double edep=totEdep* total_normalised_weight;

         //Check if the edep is not wrongly too high
         //-----------------------------------------
         G4double new_mean , new_error;
         fAccumulated_edep +=edep;
         fAccumulated_edep2 +=edep*edep;
         fNentry += 1.0;
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
           fNentry -= 1.0;
           return;
         }
         else { //accepted
           fMean_edep = new_mean;
           fError_mean_edep = new_error;
           fRelative_error =new_relative_error;
           new_mean_computed=true;
         }

     }
    if (!new_mean_computed){
         ComputeMeanEdepAndError(anEvent,fMean_edep,fError_mean_edep);
         if (fError_mean_edep>0) fRelative_error= fError_mean_edep/fMean_edep;
    }
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RMC01AnalysisManager::PrimDiffAndDirFluxForAdjointSim(
                                                          G4double prim_energy)
{ 
  G4double flux=fAmplitude_prim_spectrum;
  if ( fPrimSpectrumType ==EXPO)      flux*=std::exp(-prim_energy/fAlpha_or_E0);
  else flux*=std::pow(prim_energy, -fAlpha_or_E0);
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void  RMC01AnalysisManager::WriteHisto(G4H1* anHisto,
            G4double scaling_factor, G4String fileName, G4String header_lines)
{ std::fstream FileOutput(fileName, std::ios::out);
  FileOutput<<header_lines;
  FileOutput.setf(std::ios::scientific);
  FileOutput.precision(6);

  for (G4int i =0;i<G4int(anHisto->axis().bins());++i) {
        FileOutput<<anHisto->axis().bin_lower_edge(i)
              <<'\t'<<anHisto->axis().bin_upper_edge(i)
              <<'\t'<<anHisto->bin_height(i)*scaling_factor
              <<'\t'<<anHisto->bin_error(i)*scaling_factor<<std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RMC01AnalysisManager::WriteHisto(G4H2* anHisto,
            G4double scaling_factor, G4String fileName, G4String header_lines)
{ std::fstream FileOutput(fileName, std::ios::out);
  FileOutput<<header_lines;
  
  FileOutput.setf(std::ios::scientific);
  FileOutput.precision(6);

  for (G4int i =0;i<G4int(anHisto->axis_x().bins());++i) {
    for (G4int j =0;j<G4int(anHisto->axis_y().bins());++j) {
       FileOutput<<anHisto->axis_x().bin_lower_edge(i)
                     <<'\t'<<anHisto->axis_x().bin_upper_edge(i)
                   <<'\t'<<anHisto->axis_y().bin_lower_edge(i)
                   <<'\t'<<anHisto->axis_y().bin_upper_edge(i)
                   <<'\t'<<anHisto->bin_height(i,j)*scaling_factor
                <<'\t'<<anHisto->bin_error(i,j)*scaling_factor
                                                                  <<std::endl;
        }
  }
}
*/
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
   
   // VI: error computation now is based on number of entries and not 
   //     number of events
   // LD: This is wrong! With the use of fNentry the results were no longer
   //     correctly normalised. The mean and the error should be computed
   //     with nb_event. The old computation has been reset.
   G4float nb_event_float = G4float(nb_event);
   if (nb_event_float >1.) {
      mean = fAccumulated_edep/nb_event_float;
      G4double mean_x2 = fAccumulated_edep2/nb_event_float;
      /*
      G4cout << "Nevt= " << nb_event <<  " mean= " << mean 
             << "  mean_x2= " <<  mean_x2 << " x2 - x*x= " 
             << mean_x2-mean*mean << G4endl;
      */
      error = factor*std::sqrt(mean_x2-mean*mean)/std::sqrt(nb_event_float);
      mean *=factor;
   }
   else {
      mean=0;
      error=0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::SetPrimaryExpSpectrumForAdjointSim(
                     const G4String& particle_name, G4double omni_fluence,
                                              G4double E0, G4double Emin,
                                                G4double Emax)
{ fPrimSpectrumType = EXPO;
  if (particle_name == "e-" ) fPrimPDG_ID =
                      G4Electron::Electron()->GetPDGEncoding();
  else if (particle_name == "gamma") fPrimPDG_ID =
                                       G4Gamma::Gamma()->GetPDGEncoding();
  else if (particle_name == "proton") fPrimPDG_ID =
                                      G4Proton::Proton()->GetPDGEncoding();
  else {
   G4cout<<"The particle that you did select is not in the candidate "<<
       "list for primary [e-, gamma, proton]!"<<G4endl;
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
          G4cout<<"The particle that you did select is not in the candidate"<<
          " list for primary [e-, gamma, proton]!"<<G4endl;
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

void RMC01AnalysisManager::Book()
{
  //----------------------
  //Creation of histograms
  //----------------------

  //Energy binning of the histograms : 60 log bins over [1keV-1GeV]

  G4double emin=1.*keV;
  G4double emax=1.*GeV;

  //file_name
  fFileName[0]="forward_sim";
  if (fAdjoint_sim_mode) fFileName[0]="adjoint_sim";

  //Histo manager
   G4AnalysisManager* theHistoManager = G4AnalysisManager::Instance();
   G4String extension = theHistoManager->GetFileType();
   fFileName[1] = fFileName[0] + "." + extension;
   theHistoManager->SetFirstHistoId(1);

   G4bool fileOpen = theHistoManager->OpenFile(fFileName[0]);
   if (!fileOpen) {
     G4cout << "\n---> RMC01AnalysisManager::Book(): cannot open "
             << fFileName[1]
              << G4endl;
     return;
   }

    // Create directories
   theHistoManager->SetHistoDirectoryName("histo");

  //Histograms for :
  //        1)the forward simulation results
  //        2)the Reverse MC  simulation results normalised to a user spectrum
  //------------------------------------------------------------------------

   G4int idHisto =
      theHistoManager->CreateH1(G4String("Edep_vs_prim_ekin"),
      G4String("edep vs e- primary energy"),60,emin,emax,
      "none","none",G4String("log"));
  fEdep_vs_prim_ekin = theHistoManager->GetH1(idHisto);

  idHisto = theHistoManager->CreateH1(G4String("elecron_current"),
        G4String("electron"),60,emin,emax,
        "none","none",G4String("log"));

  fElectron_current  =  theHistoManager->GetH1(idHisto);

  idHisto= theHistoManager->CreateH1(G4String("proton_current"),
        G4String("proton"),60,emin,emax,
        "none","none",G4String("log"));
  fProton_current=theHistoManager->GetH1(idHisto);

  idHisto= theHistoManager->CreateH1(G4String("gamma_current"),
          G4String("gamma"),60,emin,emax,
          "none","none",G4String("log"));
  fGamma_current=theHistoManager->GetH1(idHisto);

  //Response matrices for the adjoint simulation only
  //-----------------------------------------------
  if (fAdjoint_sim_mode){
  //Response matrices for external isotropic e- source
  //--------------------------------------------------

   idHisto =
    theHistoManager->CreateH1(G4String("Edep_rmatrix_vs_electron_prim_energy"),
    G4String("electron RM vs e- primary energy"),60,emin,emax,
    "none","none",G4String("log"));
   fEdep_rmatrix_vs_electron_prim_energy = theHistoManager->GetH1(idHisto);

   idHisto =
      theHistoManager->
         CreateH2(G4String("Electron_current_rmatrix_vs_electron_prim_energy"),
         G4String("electron current  RM vs e- primary energy"),
         60,emin,emax,60,emin,emax,
         "none","none","none","none",G4String("log"),G4String("log"));

   fElectron_current_rmatrix_vs_electron_prim_energy =
                                             theHistoManager->GetH2(idHisto);

   idHisto =
        theHistoManager->
           CreateH2(G4String("Gamma_current_rmatrix_vs_electron_prim_energy"),
           G4String("gamma current  RM vs e- primary energy"),
           60,emin,emax,60,emin,emax,
           "none","none","none","none",G4String("log"),G4String("log"));

   fGamma_current_rmatrix_vs_electron_prim_energy =
                                           theHistoManager->GetH2(idHisto);

  //Response matrices for external isotropic gamma source

   idHisto =
    theHistoManager->CreateH1(G4String("Edep_rmatrix_vs_gamma_prim_energy"),
         G4String("electron RM vs gamma primary energy"),60,emin,emax,
        "none","none",G4String("log"));
   fEdep_rmatrix_vs_gamma_prim_energy = theHistoManager->GetH1(idHisto);

   idHisto =
        theHistoManager->
          CreateH2(G4String("Electron_current_rmatrix_vs_gamma_prim_energy"),
          G4String("electron current  RM vs gamma primary energy"),
          60,emin,emax,60,emin,emax,
        "none","none","none","none",G4String("log"),G4String("log"));

   fElectron_current_rmatrix_vs_gamma_prim_energy =
                                               theHistoManager->GetH2(idHisto);

   idHisto =
          theHistoManager->
             CreateH2(G4String("Gamma_current_rmatrix_vs_gamma_prim_energy"),
             G4String("gamma current  RM vs gamma primary energy"),
             60,emin,emax,60,emin,emax,
             "none","none","none","none",G4String("log"),G4String("log"));

   fGamma_current_rmatrix_vs_gamma_prim_energy =
                                             theHistoManager->GetH2(idHisto);

  //Response matrices for external isotropic proton source
   idHisto =
     theHistoManager->CreateH1(G4String("Edep_rmatrix_vs_proton_prim_energy"),
         G4String("electron RM vs proton primary energy"),60,emin,emax,
         "none","none",G4String("log"));
   fEdep_rmatrix_vs_proton_prim_energy = theHistoManager->GetH1(idHisto);

   idHisto =
         theHistoManager->
           CreateH2(G4String("Electron_current_rmatrix_vs_proton_prim_energy"),
           G4String("electron current  RM vs proton primary energy"),
           60,emin,emax,60,emin,emax,
         "none","none","none","none",G4String("log"),G4String("log"));

    fElectron_current_rmatrix_vs_proton_prim_energy =
                                                theHistoManager->GetH2(idHisto);

   idHisto =
           theHistoManager->
              CreateH2(G4String("Gamma_current_rmatrix_vs_proton_prim_energy"),
              G4String("gamma current  RM vs proton primary energy"),
              60,emin,emax,60,emin,emax,
              "none","none","none","none",G4String("log"),G4String("log"));

   fGamma_current_rmatrix_vs_proton_prim_energy =
                                              theHistoManager->GetH2(idHisto);

   idHisto =
              theHistoManager->
              CreateH2(G4String("Proton_current_rmatrix_vs_proton_prim_energy"),
               G4String("proton current  RM vs proton primary energy"),
               60,emin,emax,60,emin,emax,
               "none","none","none","none",G4String("log"),G4String("log"));

      fProton_current_rmatrix_vs_proton_prim_energy =
                                                theHistoManager->GetH2(idHisto);
  }
  fFactoryOn = true;
  G4cout << "\n----> Histogram Tree is opened in " << fFileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManager::Save(G4double scaling_factor)
{ if (fFactoryOn) {
    G4AnalysisManager* theHistoManager = G4AnalysisManager::Instance();
    //scaling of results
    //-----------------

    for (G4int ind=1; ind<=theHistoManager->GetNofH1s();++ind){
       theHistoManager->SetH1Ascii(ind,true);
       theHistoManager->ScaleH1(ind,scaling_factor);
    }
    for (G4int ind=1; ind<=theHistoManager->GetNofH2s();++ind)
                        theHistoManager->ScaleH2(ind,scaling_factor);

    theHistoManager->Write();
    theHistoManager->CloseFile();
    G4cout << "\n----> Histogram Tree is saved in " << fFileName[1] << G4endl;

    delete G4AnalysisManager::Instance();
    fFactoryOn = false;
  }
}
