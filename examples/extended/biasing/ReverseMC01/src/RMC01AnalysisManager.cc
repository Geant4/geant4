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
// $Id: RMC01AnalysisManager.cc,v 1.7 2010-11-11 14:39:42 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	RMC01AnalysisManager
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
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
#include "RMC01AnalysisManagerMessenger.hh"
RMC01AnalysisManager* RMC01AnalysisManager::instance = 0;

////////////////////////////////////////////////////////////////////////////////
//
RMC01AnalysisManager::RMC01AnalysisManager()  
{ 
  
  theMsg = new RMC01AnalysisManagerMessenger(this);
  
  //----------------------
  //Creation of histograms
  //----------------------
  
  //Energy binning of the histograms : 60 log bins over [1keV-1GeV]
  
  G4double bins[61];
  size_t nbin=60;
  G4double Emin=1.*keV;
  G4double Emax=1.*GeV;
  for ( size_t i=0; i <= nbin; i++) {
   	G4double val_bin;
	val_bin=Emin * std::pow(10., i * std::log10(Emax/Emin)/nbin);
	G4double exp_10=4.-int(std::log10(val_bin));
	G4double factor =std::pow(10., exp_10);
	val_bin=int(factor*val_bin)/factor;
	bins[i] = val_bin;
  
  }
  
  
  //Histograms for :
  //	1)the forward simulation results 
  //	2)the Reverse MC  simulation results normalised to a user spectrum
  //-----------------------------------------------------------------------------------
   
  edep_vs_prim_ekin = new Histo1DVar("edep_vs_prim_ekin",bins,  nbin+1, LEFT);	
  electron_current  = new Histo1DVar("electron_current",bins,  nbin+1, LEFT);	
  proton_current= new Histo1DVar("proton_current",bins,  nbin+1, LEFT);	
  gamma_current= new Histo1DVar("gamma_current",bins,  nbin+1, LEFT);	
  
 
  
  //Answer matrices for the adjoint simulation only
  //-----------------------------------------------
  
  
  //answer matrices for external isotropic e- source
  
  edep_answer_matrix_vs_electron_prim_energy = new Histo1DVar("edep_answer_matrix_vs_electron_prim_energy", bins, nbin+1, LEFT); 
  	
  electron_current_answer_matrix_vs_electron_prim_energy = new Histo2DVar("electron_current_answer_matrix_vs_electron_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  gamma_current_answer_matrix_vs_electron_prim_energy = new Histo2DVar("gamma_current_answer_matrix_vs_electron_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  
  //answer matrices for external isotropic gamma source
  
  edep_answer_matrix_vs_gamma_prim_energy = new Histo1DVar("edep_answer_matrix_vs_gamma_prim_energy", bins, nbin+1, LEFT); 
  	
  electron_current_answer_matrix_vs_gamma_prim_energy = new Histo2DVar("electron_current_answer_matrix_vs_gamma_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  gamma_current_answer_matrix_vs_gamma_prim_energy = new Histo2DVar("gamma_current_answer_matrix_vs_gamma_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  
  
  //answer matrices for external isotropic proton source
  
  edep_answer_matrix_vs_proton_prim_energy = new Histo1DVar("edep_answer_matrix_vs_proton_prim_energy", bins, nbin+1, LEFT); 
  	
  electron_current_answer_matrix_vs_proton_prim_energy = new Histo2DVar("electron_current_answer_matrix_vs_proton_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  gamma_current_answer_matrix_vs_proton_prim_energy = new Histo2DVar("gamma_current_answer_matrix_vs_proton_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  proton_current_answer_matrix_vs_proton_prim_energy = new Histo2DVar("proton_current_answer_matrix_vs_proton_prim_energy", bins, nbin+1, LEFT, bins, nbin+1, LEFT); 
  	
  
  
  //--------------------------
  // Convergence test variables
  //--------------------------
  
  precision_to_reach =0.;
  stop_run_if_precision_reached =true;
  nb_evt_modulo_for_convergence_test =5000;
  
  
  
  //-------------
  //Timer for convergence vector
  //-------------
  
  theTimer = new G4Timer();
 
  
  
  
  

  
  //---------------------------------
  //Primary spectrum for normalisation of adjoint results
  //---------------------------------
  
  thePrimPDG_ID = G4Electron::Electron()->GetPDGEncoding();
  alpha_or_E0=.5*MeV;
  amplitude_prim_spectrum =1.;
  emin_prim_spectrum=1.*keV;
  emax_prim_spectrum = 20.*MeV;
  						
  
}
////////////////////////////////////////////////////////////////////////////////
//  
RMC01AnalysisManager::~RMC01AnalysisManager() 
{ 
  delete edep_vs_prim_ekin;
  delete electron_current;
  delete proton_current;
  delete gamma_current;
  
  delete edep_answer_matrix_vs_electron_prim_energy;
  delete electron_current_answer_matrix_vs_electron_prim_energy;
  delete gamma_current_answer_matrix_vs_electron_prim_energy;
  
  delete edep_answer_matrix_vs_gamma_prim_energy;
  delete electron_current_answer_matrix_vs_gamma_prim_energy;
  delete gamma_current_answer_matrix_vs_gamma_prim_energy;
  
  delete edep_answer_matrix_vs_proton_prim_energy;
  delete electron_current_answer_matrix_vs_proton_prim_energy;
  delete proton_current_answer_matrix_vs_proton_prim_energy;
  delete gamma_current_answer_matrix_vs_proton_prim_energy;
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
RMC01AnalysisManager* RMC01AnalysisManager::GetInstance()
{
  if (instance == 0) instance = new RMC01AnalysisManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//  
void RMC01AnalysisManager::BeginOfRun(const G4Run* aRun)
{  accumulated_edep =0.;
   accumulated_edep2 =0.;
   relative_error=1.;
   mean_edep=0.;
   error_mean_edep=0.;
   adjoint_sim_mode =G4AdjointSimManager::GetInstance()->GetAdjointSimMode();
   if (adjoint_sim_mode){
   	nb_evt_per_adj_evt=aRun->GetNumberOfEventToBeProcessed()/G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
	ConvergenceFileOutput.open("ConvergenceOfAdjointSimulationResults.txt", std::ios::out);
  	ConvergenceFileOutput<<"Normalised Edep[MeV]\terror[MeV]\tcomputing_time[s]"<<std::endl;
   }
   else {
   	ConvergenceFileOutput.open("ConvergenceOfForwardSimulationResults.txt", std::ios::out);
  	ConvergenceFileOutput<<"Edep per event [MeV]\terror[MeV]\tcomputing_time[s]"<<std::endl;
   }
   ConvergenceFileOutput.setf(std::ios::scientific);
   ConvergenceFileOutput.precision(6); 	
   ResetHistograms();
   theTimer->Start();
   elapsed_time=0.;
}
////////////////////////////////////////////////////////////////////////////////
//  
void RMC01AnalysisManager::EndOfRun(const G4Run* aRun)
{ theTimer->Stop();
 
  
  if (!adjoint_sim_mode){
        
	G4cout<<"Results of forward simulation!"<<std::endl;
	G4cout<<"edep per event [MeV] = "<<mean_edep<<std::endl;
	G4cout<<"error[MeV] = "<<error_mean_edep<<std::endl;
	G4int nb_evt=aRun->GetNumberOfEvent();
 	WriteHisto(edep_vs_prim_ekin,1./nb_evt, G4String("Fwd_Edep_vs_EkinPrim.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tEdep[MeV]\terr_Edep[MeV]\n"));
 	WriteHisto(electron_current,1./nb_evt, G4String("Fwd_ElectronCurrent.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering electron\terr\n"));
 	WriteHisto(proton_current,1./nb_evt, G4String("Fwd_ProtonCurrent.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering proton\terr[]\n"));
  	WriteHisto(gamma_current,1./nb_evt, G4String("Fwd_GammaCurrent.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering gamma\terr[]\n"));
  }
  
  else {
  	G4cout<<"Results of reverse/adjoint simulation!"<<std::endl;
	G4cout<<"normalised edep [MeV] = "<<mean_edep<<std::endl;
	G4cout<<"error[MeV] = "<<error_mean_edep<<std::endl;

	
	G4double factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun()*nb_evt_per_adj_evt/aRun->GetNumberOfEvent();
		
	
	
	WriteHisto(edep_vs_prim_ekin,factor, G4String("Adj_Edep_vs_EkinPrim.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tEdep[MeV]\terr_Edep[MeV]\n"));
 	WriteHisto(electron_current,factor, G4String("Adj_ElectronCurrent.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering electron\terr\n"));
 	WriteHisto(proton_current,factor, G4String("Adj_ProtonCurrent.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering proton\terr[]\n"));
  	WriteHisto(gamma_current,factor, G4String("Adj_GammaCurrent.txt"),G4String("E1[MeV]\t\tE2[MeV]\t\tnb entering gamma\terr[]\n"));
  
  	
	WriteHisto(edep_answer_matrix_vs_electron_prim_energy,factor,
					   G4String("Adj_Edep_vs_EkinPrimElectron_Answer.txt"),
					   G4String("E1[MeV]\t\tE2[MeV]\t\tEdep Efficiency[MeV*cm2*MeV*str]\terr_Edep[MeV*cm2*MeV*str]\n"));
 	
  	WriteHisto(electron_current_answer_matrix_vs_electron_prim_energy,factor,
					   G4String("Adj_ElectronCurrent_vs_EkinPrimElectron_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
 	
	WriteHisto(gamma_current_answer_matrix_vs_electron_prim_energy,factor,
					   G4String("Adj_GammaCurrent_vs_EkinPrimElectron_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
 	
  	
	
	WriteHisto(edep_answer_matrix_vs_gamma_prim_energy,factor,
					   G4String("Adj_Edep_vs_EkinPrimGamma_Answer.txt"),
					   G4String("E1[MeV]\t\tE2[MeV]\t\tEdep Efficiency[MeV*cm2*MeV*str]\terr_Edep[MeV*cm2*MeV*str]\n"));
 	
  	WriteHisto(electron_current_answer_matrix_vs_gamma_prim_energy,factor,
					   G4String("Adj_ElectronCurrent_vs_EkinPrimGamma_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
 	
	WriteHisto(gamma_current_answer_matrix_vs_gamma_prim_energy,factor,
					   G4String("Adj_GammaCurrent_vs_EkinPrimGamma_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
 	
	
	
	WriteHisto(edep_answer_matrix_vs_proton_prim_energy,factor,
					   G4String("Adj_Edep_vs_EkinPrimProton_Answer.txt"),
					   G4String("E1[MeV]\t\tE2[MeV]\t\tEdep Efficiency[MeV*cm2*MeV*str]\terr_Edep[MeV*cm2*MeV*str]\n"));
 	
  	WriteHisto(electron_current_answer_matrix_vs_proton_prim_energy,factor,
					   G4String("Adj_ElectronCurrent_vs_EkinPrimProton_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
 	
	WriteHisto(gamma_current_answer_matrix_vs_proton_prim_energy,factor,
					   G4String("Adj_GammaCurrent_vs_EkinPrimProton_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
 	
	WriteHisto(proton_current_answer_matrix_vs_proton_prim_energy,factor,
					   G4String("Adj_ProtonCurrent_vs_EkinPrimProton_Answer.txt"),
					   G4String("Eprim1[MeV]\t\tEprim2[MeV]\t\tEsec1[MeV]\t\tEsec2[MeV]\t Current Efficiency[cm2*MeV*str]\terr[cm2*MeV*str]\n"));
  }
  ConvergenceFileOutput.close();
}
////////////////////////////////////////////////////////////////////////////////
//  
void RMC01AnalysisManager::BeginOfEvent(const G4Event* )
{;   
}
////////////////////////////////////////////////////////////////////////////////
//  
void RMC01AnalysisManager::EndOfEvent(const G4Event* anEvent)
{  



   if (adjoint_sim_mode) EndOfEventForAdjointSimulation(anEvent);
   else EndOfEventForForwardSimulation(anEvent);
   
   
   //Test convergence. The error is already computed
   //--------------------------------------
   G4int nb_event=anEvent->GetEventID()+1;
   G4double factor=1.;
   if (adjoint_sim_mode) {
   	G4double  n_adj_evt= nb_event/nb_evt_per_adj_evt;
	// nb_event/nb_evt_per_adj_evt;
	if (n_adj_evt*nb_evt_per_adj_evt == nb_event) {
		nb_event =static_cast<G4int>(n_adj_evt);
		factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
	}	
	else nb_event=0;
	
   }
   
   
   
   if (nb_event>100 && stop_run_if_precision_reached && precision_to_reach >relative_error) {
		G4cout<<precision_to_reach*100.<<"%  Precision reached!"<<std::endl;
		theTimer->Stop();
		elapsed_time+=theTimer->GetRealElapsed();
		ConvergenceFileOutput<<mean_edep<<'\t'<<error_mean_edep<<'\t'<<elapsed_time<<std::endl;
		G4RunManager::GetRunManager()->AbortRun(true);
   }
   
   
   if (nb_event>0 && nb_event % nb_evt_modulo_for_convergence_test == 0) {
   	theTimer->Stop();
	elapsed_time+=theTimer->GetRealElapsed();
	theTimer->Start();
	ConvergenceFileOutput<<mean_edep<<'\t'<<error_mean_edep<<'\t'<<elapsed_time<<std::endl;
	
   }	
   	
  
   
  
  
   
}   
////////////////////////////////////////////////////////////////////////////////
// 
void  RMC01AnalysisManager::EndOfEventForForwardSimulation(const G4Event* anEvent)
{  
   
   G4SDManager* SDman = G4SDManager::GetSDMpointer();
   G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
   RMC01DoubleWithWeightHitsCollection* EdepCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("edep")));
   RMC01DoubleWithWeightHitsCollection* ElectronCurrentCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_electron")));
   RMC01DoubleWithWeightHitsCollection* ProtonCurrentCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_proton")));
   RMC01DoubleWithWeightHitsCollection* GammaCurrentCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_gamma")));
   
   //Total energy deposited in Event
   //-------------------------------
   G4double totEdep=0; 
   G4int i;
   for (i=0;i<EdepCollection->entries();i++) totEdep+=(*EdepCollection)[i]->GetValue()*(*EdepCollection)[i]->GetWeight();
   
   if (totEdep>0.){
   	accumulated_edep +=totEdep ;
   	accumulated_edep2 +=totEdep*totEdep;
        G4PrimaryParticle* thePrimary=anEvent->GetPrimaryVertex()->GetPrimary();
   	G4double E0= thePrimary->GetG4code()->GetPDGMass();
   	G4double P=thePrimary->GetMomentum().mag();
   	G4double prim_ekin =std::sqrt(E0*E0+P*P)-E0;
   	edep_vs_prim_ekin->fill(prim_ekin,totEdep);
   } 
   ComputeMeanEdepAndError(anEvent,mean_edep,error_mean_edep);
   if (error_mean_edep>0) relative_error= error_mean_edep/mean_edep;
   		
   
   
   //Particle current on sensitive cylinder
   //-------------------------------------
   
   for (i=0;i<ElectronCurrentCollection->entries();i++) {
   	G4double ekin =(*ElectronCurrentCollection)[i]->GetValue();
	G4double weight=(*ElectronCurrentCollection)[i]->GetWeight();
	electron_current->fill(ekin,weight);
   }
   
   for (i=0;i<ProtonCurrentCollection->entries();i++) {
   	G4double ekin =(*ProtonCurrentCollection)[i]->GetValue();
	G4double weight=(*ProtonCurrentCollection)[i]->GetWeight();
	proton_current->fill(ekin,weight);
   }	
   
   for (i=0;i<GammaCurrentCollection->entries();i++) {
   	G4double ekin =(*GammaCurrentCollection)[i]->GetValue();
	G4double weight=(*GammaCurrentCollection)[i]->GetWeight();
	gamma_current->fill(ekin,weight);
   }
   
   
}
////////////////////////////////////////////////////////////////////////////////
// 
void  RMC01AnalysisManager::EndOfEventForAdjointSimulation(const G4Event* anEvent)
{  
  //Output from Sensitive volume computed during the forward tracking phase
  //-----------------------------------------------------------------------
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
  RMC01DoubleWithWeightHitsCollection* EdepCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("edep")));
  RMC01DoubleWithWeightHitsCollection* ElectronCurrentCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_electron")));
  RMC01DoubleWithWeightHitsCollection* ProtonCurrentCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_proton")));
  RMC01DoubleWithWeightHitsCollection* GammaCurrentCollection =  (RMC01DoubleWithWeightHitsCollection*)(HCE->GetHC(SDman->GetCollectionID("current_gamma")));
  
  //Output from adjoint tracking phase
  //----------------------------------------------------------------------------
  
  G4AdjointSimManager* theAdjointSimManager = G4AdjointSimManager::GetInstance();
  G4int pdg_nb =theAdjointSimManager->GetFwdParticlePDGEncodingAtEndOfLastAdjointTrack();
  G4double prim_ekin=theAdjointSimManager->GetEkinAtEndOfLastAdjointTrack();
  G4double adj_weight=theAdjointSimManager->GetWeightAtEndOfLastAdjointTrack();
 
  
  //Factor of normalisation to user selected primary spectrum (power law or exponential) 
  //------------------------------------------------------------------------------------
  
  G4double normalised_weight = 0.;
  if (pdg_nb== thePrimPDG_ID && prim_ekin>= emin_prim_spectrum && prim_ekin<= emax_prim_spectrum) 
  			normalised_weight = adj_weight*PrimDiffAndDirectionalFluxForAdjointSim(prim_ekin);
  
  
  //Answer matrices
  //-------------
   Histo1DVar* edep_answer_matrix =0;
   Histo2DVar* electron_current_answer_matrix =0;
   Histo2DVar* gamma_current_answer_matrix =0;
   Histo2DVar* proton_current_answer_matrix =0;
	
   
   if (pdg_nb == G4Electron::Electron()->GetPDGEncoding()){ //electron answer matrices
	edep_answer_matrix = edep_answer_matrix_vs_electron_prim_energy;
	electron_current_answer_matrix = electron_current_answer_matrix_vs_electron_prim_energy;
	gamma_current_answer_matrix = gamma_current_answer_matrix_vs_electron_prim_energy;
   }
   else if (pdg_nb == G4Gamma::Gamma()->GetPDGEncoding()){ //electron answer matrices
	edep_answer_matrix = edep_answer_matrix_vs_gamma_prim_energy;
	electron_current_answer_matrix = electron_current_answer_matrix_vs_gamma_prim_energy;
	gamma_current_answer_matrix = gamma_current_answer_matrix_vs_gamma_prim_energy;
   }
   else if (pdg_nb == G4Proton::Proton()->GetPDGEncoding()){ //electron answer matrices
	edep_answer_matrix = edep_answer_matrix_vs_proton_prim_energy;
	electron_current_answer_matrix = electron_current_answer_matrix_vs_proton_prim_energy;
	gamma_current_answer_matrix = gamma_current_answer_matrix_vs_proton_prim_energy;
	proton_current_answer_matrix = proton_current_answer_matrix_vs_proton_prim_energy;
   }
  
  
  
  //Registering of total energy deposited in Event
  //-------------------------------
   G4double totEdep=0; 
   G4int i;
   for (i=0;i<EdepCollection->entries();i++) totEdep+=(*EdepCollection)[i]->GetValue()*(*EdepCollection)[i]->GetWeight();
   
   G4bool new_mean_computed=false;
   if (totEdep>0.){
   	if (normalised_weight>0.){
		G4double edep=totEdep* normalised_weight;
		
		//Check if the edep is not wrongly too high
		//----------------------------------------- 
		G4double new_mean , new_error;
		
		accumulated_edep +=edep;
   		accumulated_edep2 +=edep*edep;
		
		ComputeMeanEdepAndError(anEvent,new_mean,new_error);
		G4double new_relative_error = 1.;
		if ( new_error >0) new_relative_error = new_error/ new_mean;
		if (relative_error <0.10 && new_relative_error>1.5*relative_error) { // rejected
			G4cout<<"Potential wrong adjoint weight!"<<std::endl;
			G4cout<<"The results of this event will not be registered!"<<std::endl;
			G4cout<<"previous mean edep [MeV] "<< mean_edep<<std::endl;
			G4cout<<"previous relative error "<< relative_error<<std::endl;
			G4cout<<"new rejected mean edep [MeV] "<< new_mean<<std::endl;
			G4cout<<"new rejected relative error "<< new_relative_error<<std::endl;
			accumulated_edep -=edep;
			accumulated_edep2 -=edep*edep;
			return; 
		}
		else { //accepted
			mean_edep = new_mean;
			error_mean_edep = new_error;
			relative_error =new_relative_error;
			new_mean_computed=true;
		}	 
		edep_vs_prim_ekin->fill(prim_ekin,edep); 
	}
	
	// Registering answer matrix
	//---------------------------
	
	edep_answer_matrix->fill(prim_ekin,totEdep*adj_weight/cm2);
	

   }
   if (!new_mean_computed){
   	 ComputeMeanEdepAndError(anEvent,mean_edep,error_mean_edep);
   	 if (error_mean_edep>0) relative_error= error_mean_edep/mean_edep;
   }	 
    
   
  //Registering of current of particles on the sensitive volume
  //------------------------------------------------------------
   
   for (i=0;i<ElectronCurrentCollection->entries();i++) {
   	G4double ekin =(*ElectronCurrentCollection)[i]->GetValue();
	G4double weight=(*ElectronCurrentCollection)[i]->GetWeight();
	electron_current->fill(ekin,weight*normalised_weight);
	electron_current_answer_matrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
   }
   
   for (i=0;i<ProtonCurrentCollection->entries();i++) {
   	G4double ekin =(*ProtonCurrentCollection)[i]->GetValue();
	G4double weight=(*ProtonCurrentCollection)[i]->GetWeight();
	proton_current->fill(ekin,weight*normalised_weight);
	proton_current_answer_matrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);
  }	
   
   for (i=0;i<GammaCurrentCollection->entries();i++) {
   	G4double ekin =(*GammaCurrentCollection)[i]->GetValue();
	G4double weight=(*GammaCurrentCollection)[i]->GetWeight();
	gamma_current->fill(ekin,weight*normalised_weight);
	gamma_current_answer_matrix->fill(prim_ekin,ekin,weight*adj_weight/cm2);

   }
   		
  
}
////////////////////////////////////////////////////////////////////////////////
//
G4double RMC01AnalysisManager::PrimDiffAndDirectionalFluxForAdjointSim(G4double prim_energy)
{ 
  G4double flux=amplitude_prim_spectrum;
  if ( the_prim_spectrum_type ==EXPO)	flux*=std::exp(-prim_energy/alpha_or_E0);
  else flux*=std::pow(prim_energy, -alpha_or_E0);
  return flux;
}
////////////////////////////////////////////////////////////////////////////////
//
void  RMC01AnalysisManager::WriteHisto(Histo1DVar* anHisto, G4double scaling_factor, G4String fileName, G4String header_lines)
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
////////////////////////////////////////////////////////////////////////////////
//
void  RMC01AnalysisManager::WriteHisto(Histo2DVar* anHisto, G4double scaling_factor, G4String fileName, G4String header_lines)
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
		          <<'\t'<<anHisto->get_bin_error(int(i),int(j))*scaling_factor<<std::endl;
   	}
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01AnalysisManager::ResetHistograms()
{ edep_vs_prim_ekin->reset();
  electron_current->reset();
  proton_current->reset();
  gamma_current->reset();
  
  edep_answer_matrix_vs_electron_prim_energy->reset();
  electron_current_answer_matrix_vs_electron_prim_energy->reset();
  gamma_current_answer_matrix_vs_electron_prim_energy->reset();
  
  edep_answer_matrix_vs_gamma_prim_energy->reset();
  electron_current_answer_matrix_vs_gamma_prim_energy->reset();
  gamma_current_answer_matrix_vs_gamma_prim_energy->reset();
  
  edep_answer_matrix_vs_proton_prim_energy->reset();
  electron_current_answer_matrix_vs_proton_prim_energy->reset();
  proton_current_answer_matrix_vs_proton_prim_energy->reset();
  gamma_current_answer_matrix_vs_proton_prim_energy->reset();
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01AnalysisManager::ComputeMeanEdepAndError(const G4Event* anEvent,G4double& mean,G4double& error)
{  
   G4int nb_event=anEvent->GetEventID()+1;
   G4double factor=1.;
   if (adjoint_sim_mode) {
   	nb_event /=nb_evt_per_adj_evt;
	factor=1.*G4AdjointSimManager::GetInstance()->GetNbEvtOfLastRun();
   }
   
   //error computation
   if (nb_event>1) {
   	  mean = accumulated_edep/nb_event;
	  G4double mean_x2 =accumulated_edep2/nb_event;
  	  error = factor*std::sqrt(mean_x2-mean*mean)/std::sqrt(G4double(nb_event));
	  mean *=factor;
   } else {
          mean=0;
          error=0;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01AnalysisManager::SetPrimaryExponentialSpectrumForAdjointSim(const G4String& particle_name, G4double omni_fluence,G4double E0, G4double Emin,G4double Emax)
{ the_prim_spectrum_type = EXPO;
  if (particle_name == "e-" ) thePrimPDG_ID = G4Electron::Electron()->GetPDGEncoding();
  else if (particle_name == "gamma") thePrimPDG_ID = G4Gamma::Gamma()->GetPDGEncoding();
  else if (particle_name == "proton") thePrimPDG_ID = G4Proton::Proton()->GetPDGEncoding();
  else {
  	G4cout<<"The particle that you did select is not in the candidate list for primary [e-, gamma, proton]!"<<G4endl;
  	return;
  }	
  

  alpha_or_E0 = E0 ;
  amplitude_prim_spectrum = omni_fluence/E0/(std::exp(-Emin/E0)-std::exp(-Emax/E0))/4./pi;
  emin_prim_spectrum = Emin ;
  emax_prim_spectrum = Emax;
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01AnalysisManager::SetPrimaryPowerLawSpectrumForAdjointSim(const G4String& particle_name, G4double omni_fluence,G4double alpha, G4double Emin,G4double Emax)
{ the_prim_spectrum_type  =POWER;
  if (particle_name == "e-" ) thePrimPDG_ID = G4Electron::Electron()->GetPDGEncoding();
  else if (particle_name == "gamma") thePrimPDG_ID = G4Gamma::Gamma()->GetPDGEncoding();
  else if (particle_name == "proton") thePrimPDG_ID = G4Proton::Proton()->GetPDGEncoding();
  else {
  	G4cout<<"The particle that you did select is not in the candidate list for primary [e-, gamma, proton]!"<<G4endl;
  	return;
  }	
  

 if (alpha ==1.) {
 	amplitude_prim_spectrum = omni_fluence/std::log(Emax/Emin)/4./pi;
 }
 else {
 	G4double p=1.-alpha;
 	amplitude_prim_spectrum = omni_fluence/p/(std::pow(Emax,p)-std::pow(Emin,p))/4./pi;
 	
 }

  alpha_or_E0 = alpha ;
  emin_prim_spectrum = Emin ;
  emax_prim_spectrum = Emax;
}
